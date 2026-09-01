"""Interactive analysis helpers for GlobalLightFitColumns outputs."""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def _tree_frame(root_file, tree_name, branches=None):
    import ROOT
    source = ROOT.TFile.Open(str(root_file), "READ")
    if not source or source.IsZombie():
        raise OSError(f"Cannot open {root_file}")
    tree = source.Get(tree_name)
    if not tree:
        source.Close()
        raise KeyError(f"Missing tree {tree_name} in {root_file}")
    available = {branch.GetName() for branch in tree.GetListOfBranches()}
    names = sorted(available) if branches is None else [b for b in branches if b in available]
    arrays = ROOT.RDataFrame(tree).AsNumpy(names)
    source.Close()
    return pd.DataFrame({name: np.asarray(arrays[name]) for name in names})


def load_studies(files):
    """Load any number of labelled column-fit ROOT files."""
    studies = {}
    for label, filename in files.items():
        path = Path(filename).expanduser().resolve()
        if not path.exists():
            print(f"Skipping {label}: missing {path}")
            continue
        fit = _tree_frame(path, "global_fit")
        columns = _tree_frame(path, "global_fit_column_charges")
        input_hits = _tree_frame(path, "fiber_hits", ["event", "charge"])
        dbscan_hits = _tree_frame(path, "global_fit_dbscan_fibres", ["event"])
        studies[label] = {"path": path, "fit": fit, "columns": columns,
                          "input_hits": input_hits, "dbscan_hits": dbscan_hits}
    if not studies:
        raise FileNotFoundError("None of the configured column-fit files exists")
    return studies


def load_truth(reference_file):
    segment_branches = ["event", "segment", "detector", "primary_id", "primary_pdg", "contributors",
                        "contributor_track_ids", "contributor_pdgs",
                        "cube_x", "cube_y", "cube_z", "start_x", "start_y", "start_z",
                        "stop_x", "stop_y", "stop_z", "energy_deposit", "track_length"]
    track_branches = ["event", "track_id", "parent_id", "pdg", "point", "x", "y", "z"]
    return (_tree_frame(reference_file, "mc_virtual_segments", segment_branches),
            _tree_frame(reference_file, "mc_track_points", track_branches))


def primary_muon_ids(track_points):
    muons = track_points[np.abs(track_points.pdg) == 13]
    return muons.groupby("event").track_id.apply(lambda values: set(values.astype(int))).to_dict()


def select_truth_segments(segments, track_points, mode="all"):
    selected = segments[segments.detector == 0].copy()
    if mode == "all":
        return selected
    if mode != "primary_muon":
        raise ValueError("truth mode must be 'all' or 'primary_muon'")
    ids = primary_muon_ids(track_points)
    mask = [int(primary) in ids.get(int(event), set())
            for event, primary in zip(selected.event, selected.primary_id)]
    return selected[np.asarray(mask)]


def match_columns_to_truth(columns, segments, track_points, mode="all"):
    """Match MC energy deposits to each fitted physical column interval."""
    truth = select_truth_segments(segments, track_points, mode)
    midpoint = {
        0: 0.5 * (truth.start_x.to_numpy() + truth.stop_x.to_numpy()),
        1: 0.5 * (truth.start_y.to_numpy() + truth.stop_y.to_numpy()),
        2: 0.5 * (truth.start_z.to_numpy() + truth.stop_z.to_numpy()),
    }
    truth = truth.assign(_mx=midpoint[0], _my=midpoint[1], _mz=midpoint[2])
    coordinate = ["_mx", "_my", "_mz"]
    truth_by_event = {int(event): frame for event, frame in truth.groupby("event", sort=False)}
    rows = []
    for event, event_columns in columns.groupby("event", sort=False):
        event_truth = truth_by_event.get(int(event), truth.iloc[0:0])
        if event_columns.empty:
            continue
        axis = int(event_columns.axis.iloc[0])
        values = event_truth[coordinate[axis]].to_numpy()
        deposits = event_truth.energy_deposit.to_numpy()
        for row in event_columns.itertuples(index=False):
            inside = (values >= row.column_low) & (values < row.column_high)
            data = row._asdict()
            data["mc_energy"] = float(deposits[inside].sum())
            data["mc_segments"] = int(inside.sum())
            rows.append(data)
    return pd.DataFrame(rows)


def charge_calibration(matched):
    """Least-squares charge/energy scale through the origin."""
    energy = matched.mc_energy.to_numpy(float)
    charge = matched.fitted_selected_fibre_charge.to_numpy(float)
    denominator = np.dot(energy, energy)
    return np.dot(energy, charge) / denominator if denominator > 0 else np.nan


def add_calibrated_residuals(matched, scale=None):
    result = matched.copy()
    if scale is None:
        scale = charge_calibration(result)
    result["charge_per_energy"] = scale
    result["fitted_energy"] = result.fitted_selected_fibre_charge / scale
    result["energy_residual"] = result.fitted_energy - result.mc_energy
    result["energy_fractional_residual"] = np.where(
        result.mc_energy > 0, result.energy_residual / result.mc_energy, np.nan)
    return result, scale


def mip_prior_variants(studies, matched_by_study):
    """Return regularized/unregularized matched frames when MIP data exist.

    Files produced before the MIP prior are silently skipped.  Both variants
    use the stored independent COLUMN_CHARGE_PER_ENERGY scale when available;
    otherwise one scale fitted from the regularized columns is shared.
    """
    variants = {}
    for label, matched in matched_by_study.items():
        if "unregularized_fitted_selected_fibre_charge" not in matched:
            continue
        fit = studies[label]["fit"]
        enabled = ("column_mip_prior_enabled" in fit and
                   bool(np.any(fit.column_mip_prior_enabled.to_numpy())))
        if not enabled:
            continue
        scale = np.nan
        if "column_mip_prior_charge_per_energy" in fit:
            positive = fit.column_mip_prior_charge_per_energy.to_numpy(float)
            positive = positive[np.isfinite(positive) & (positive > 0)]
            if len(positive):
                scale = float(np.median(positive))
        if not np.isfinite(scale) or scale <= 0:
            scale = charge_calibration(matched)
        regularized, _ = add_calibrated_residuals(matched, scale)
        unregularized_source = matched.copy()
        unregularized_source["fitted_selected_fibre_charge"] = (
            unregularized_source.unregularized_fitted_selected_fibre_charge)
        unregularized, _ = add_calibrated_residuals(unregularized_source, scale)
        variants[label] = {"regularized": regularized,
                           "unregularized": unregularized, "scale": scale}
    return variants


def mip_prior_summary(studies, variants):
    """Event and column summary; returns an empty frame for legacy files."""
    rows = []
    for label, values in variants.items():
        fit = studies[label]["fit"]
        columns = values["regularized"]
        eligible = (columns.mip_prior_eligible.astype(bool)
                    if "mip_prior_eligible" in columns else
                    np.zeros(len(columns), bool))
        shift = (values["regularized"].fitted_energy.to_numpy() -
                 values["unregularized"].fitted_energy.to_numpy())
        rows.append({
            "lightmap": label, "events": len(fit),
            "prior_converged": int(fit.column_mip_prior_converged.sum()),
            "prior_stalled": int(fit.column_mip_prior_stalled.sum()),
            "prior_converged_fraction": float(
                fit.column_mip_prior_converged.sum()/max(len(fit), 1)),
            "eligible_columns": int(eligible.sum()),
            "charge_per_energy": values["scale"],
            "mean_regularization_shift": float(np.mean(shift)),
            "rms_regularization_shift": float(np.sqrt(np.mean(shift**2))),
            "median_prior_iterations": float(
                fit.column_mip_prior_iterations.median()),
            "median_prior_penalty": float(fit.column_mip_prior_penalty.median()),
        })
    return pd.DataFrame(rows).set_index("lightmap") if rows else pd.DataFrame()


def plot_mip_charge_comparison(variants, bins=80, energy_range=None,
                               fractional_range=(-2., 2.),
                               minimum_fractional_truth_energy=.5):
    """Plot truth, regularized, and retained unregularized column charges."""
    if not variants:
        print("No enabled MIP-prior branches found; skipping comparison plots")
        return None, pd.DataFrame()
    if energy_range is None:
        all_values = []
        for values in variants.values():
            for key in ("regularized", "unregularized"):
                all_values.extend(values[key].fitted_energy.replace(
                    [np.inf, -np.inf], np.nan).dropna())
            all_values.extend(values["regularized"].mc_energy.replace(
                [np.inf, -np.inf], np.nan).dropna())
        energy_range = (0., max(float(np.quantile(all_values, .995)), 1e-12))
    labels = list(variants)
    figure, axes = plt.subplots(len(labels), 4, figsize=(22, 4.2*len(labels)),
                                squeeze=False, constrained_layout=True)
    summary_rows = []
    for row_axes, label in zip(axes, labels):
        regularized = variants[label]["regularized"]
        unregularized = variants[label]["unregularized"]
        energy_edges = np.linspace(*energy_range, int(bins)+1)
        row_axes[0].hist(regularized.mc_energy, bins=energy_edges,
                         histtype="step", density=True, color="black",
                         linewidth=2, label="MC truth")
        row_axes[0].hist(unregularized.fitted_energy, bins=energy_edges,
                         histtype="step", density=True, color="#4c78a8",
                         label="without regularization")
        row_axes[0].hist(regularized.fitted_energy, bins=energy_edges,
                         histtype="step", density=True, color="#e45756",
                         label="with MIP regularization")
        row_axes[0].set(title=label, xlabel="column energy", ylabel="density")
        stable_fraction = regularized.mc_energy > minimum_fractional_truth_energy
        for frame, colour, curve_label in (
                (unregularized, "#4c78a8", "without regularization"),
                (regularized, "#e45756", "with MIP regularization")):
            absolute = frame.energy_residual.replace(
                [np.inf, -np.inf], np.nan).dropna()
            row_axes[1].hist(absolute, bins=bins, histtype="step",
                             density=True, color=colour, label=curve_label)
            residual = frame.loc[stable_fraction, "energy_fractional_residual"].replace(
                [np.inf, -np.inf], np.nan).dropna()
            row_axes[2].hist(residual, bins=bins, range=fractional_range,
                             histtype="step", density=True, color=colour,
                             label=curve_label)
            summary_rows.append({"lightmap": label, "solution": curve_label,
                                 "columns": len(residual),
                                 "absolute_residual_mean": absolute.mean(),
                                 "absolute_residual_std": absolute.std(),
                                 "absolute_residual_rms": np.sqrt(
                                     np.mean(absolute.to_numpy()**2)),
                                 "fractional_residual_mean": residual.mean(),
                                 "fractional_residual_std": residual.std(),
                                 "fractional_residual_rms": np.sqrt(
                                     np.mean(residual.to_numpy()**2))})
        row_axes[1].set(xlabel="reconstructed - MC column energy",
                        ylabel="density")
        row_axes[2].set(xlabel="(reconstructed - MC) / MC", ylabel="density",
                        xlim=fractional_range)
        shift = regularized.fitted_energy-unregularized.fitted_energy
        eligible = (regularized.mip_prior_eligible.astype(bool)
                    if "mip_prior_eligible" in regularized else
                    np.ones(len(regularized), bool))
        row_axes[3].hist(shift[~eligible], bins=bins, histtype="step",
                         density=True, color="0.5", label="not eligible")
        row_axes[3].hist(shift[eligible], bins=bins, histtype="step",
                         density=True, color="#e45756", label="prior eligible")
        row_axes[3].axvline(0., color="black", linestyle="--", linewidth=1)
        row_axes[3].set(xlabel="regularized - unregularized column energy",
                        ylabel="density")
        for axis in row_axes:
            axis.grid(alpha=.2); axis.legend(fontsize=8)
    return figure, pd.DataFrame(summary_rows)


def plot_mip_event_profiles(variants, event):
    """Overlay MC and both reconstructed solutions for one event."""
    if not variants:
        print("No enabled MIP-prior branches found; skipping event profiles")
        return None
    figure, axes = plt.subplots(len(variants), 1,
                                figsize=(12, 3.8*len(variants)),
                                squeeze=False, constrained_layout=True)
    for axis, (label, values) in zip(axes.flat, variants.items()):
        regularized = values["regularized"]
        regularized = regularized[regularized.event == event].sort_values(
            "column_index")
        unregularized = values["unregularized"]
        unregularized = unregularized[unregularized.event == event].sort_values(
            "column_index")
        if regularized.empty:
            axis.text(.5, .5, f"Event {event} unavailable", ha="center",
                      transform=axis.transAxes)
            continue
        coordinate = (.5*(regularized.column_low+regularized.column_high)).to_numpy(float)
        axis.step(coordinate, regularized.mc_energy.to_numpy(float), where="mid", color="black",
                  linewidth=2, label="MC truth")
        axis.step(coordinate, unregularized.fitted_energy.to_numpy(float), where="mid",
                  color="#4c78a8", label="without regularization")
        axis.step(coordinate, regularized.fitted_energy.to_numpy(float), where="mid",
                  color="#e45756", label="with MIP regularization")
        axis.set(title=f"{label}, event {event}",
                 xlabel="dominant-axis column centre [mm]",
                 ylabel="column energy")
        axis.grid(alpha=.2); axis.legend(fontsize=8)
    return figure


def convergence_table(studies):
    rows = []
    for label, study in studies.items():
        fit = study["fit"]
        available = fit.column_fit_available.astype(bool) if "column_fit_available" in fit else np.zeros(len(fit), bool)
        converged = fit.column_fit_converged.astype(bool) if "column_fit_converged" in fit else np.zeros(len(fit), bool)
        rows.append({
            "lightmap": label,
            "events": len(fit),
            "geometry_status_0": int((fit.status == 0).sum()),
            "geometry_status_3": int((fit.status == 3).sum()),
            "geometry_status_5": int((fit.status == 5).sum()),
            "column_results": int(available.sum()),
            "column_converged": int(converged.sum()),
            "column_converged_fraction": float(converged.sum() / max(available.sum(), 1)),
            "median_edm": float(fit.edm.median()) if "edm" in fit else np.nan,
            "median_column_iterations": float(fit.column_fit_iterations.median()) if "column_fit_iterations" in fit else np.nan,
        })
    return pd.DataFrame(rows).set_index("lightmap")


def plot_convergence(studies):
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for label, study in studies.items():
        fit = study["fit"]
        axes[0, 0].hist(fit.edm, bins=np.logspace(-7, 0, 60), histtype="step", label=label)
        if "column_fit_iterations" in fit:
            axes[0, 1].hist(fit.column_fit_iterations, bins=50, histtype="step", label=label)
        if "seed_nll" in fit:
            axes[1, 0].hist(fit.seed_nll - fit.nll, bins=60, histtype="step", label=label)
        if "column_nll" in fit:
            axes[1, 1].hist(fit.nll - fit.column_nll, bins=60, histtype="step", label=label)
    axes[0, 0].set(xscale="log", xlabel="geometry EDM", ylabel="events")
    axes[0, 1].set(xlabel="column EM iterations", ylabel="events")
    axes[1, 0].set(xlabel="seed NLL - geometry NLL", ylabel="events")
    axes[1, 1].set(xlabel="geometry NLL - column NLL", ylabel="events")
    for axis in axes.flat: axis.legend(fontsize=8)
    fig.tight_layout()
    return fig


def event_cutflow(studies,residuals=None,convergence="column"):
    """Count strictly sequential analysis losses plus parallel fit diagnostics."""
    if convergence not in (None,"none","column","geometry","both"):
        raise ValueError("convergence must be None, 'none', 'column', 'geometry', or 'both'")
    event_sets={};rows=[];sequences={}
    for label,study in studies.items():
        fit=study["fit"];hits=study["input_hits"]
        input_events=set(hits.event.astype(int))
        minimum_charge=float(fit.minimum_charge.iloc[0]) if len(fit) and "minimum_charge" in fit else 0.
        threshold_events=set(hits.loc[hits.charge>=minimum_charge,"event"].astype(int))
        dbscan_events=set(study["dbscan_hits"].event.astype(int))
        fit_events=set(fit.event.astype(int))
        geometry_events=set(fit.loc[fit.status==0,"event"].astype(int))
        available_events=set(fit.loc[fit.column_fit_available==1,"event"].astype(int))
        column_events=set(fit.loc[fit.column_fit_converged==1,"event"].astype(int))
        both_events=geometry_events&column_events
        if convergence in (None,"none"):selected=fit_events
        elif convergence=="geometry":selected=geometry_events
        elif convergence=="column":selected=column_events
        else:selected=both_events
        residual_events=(set(residuals[label].event.astype(int)) if residuals is not None and label in residuals else fit_events)
        analysis_events=selected&residual_events
        event_sets[label]=analysis_events
        sequences[label]=[("input events",input_events),(f"charge >= {minimum_charge:g}",threshold_events),
                          ("after DBSCAN",dbscan_events),("fit result written",fit_events),
                          (f"selected: {convergence or 'none'} convergence",selected),
                          ("selected with MC residual",analysis_events)]
        diagnostics=[("geometry status 0",geometry_events),("column fit available",available_events),
                     ("column fit converged",column_events),("geometry and column converged",both_events)]
        for stage,events in diagnostics:
            rows.append({"lightmap":label,"category":"diagnostic","stage":stage,"events":len(events),
                         "lost_from_previous":np.nan,"fraction_of_previous":np.nan,
                         "fraction_of_input":len(events)/max(len(input_events),1)})
    common=set.intersection(*event_sets.values()) if event_sets else set()
    for label,sequence in sequences.items():
        sequence=sequence+[("common across all lightmaps",common)];previous=None;input_count=len(sequence[0][1])
        for stage,events in sequence:
            count=len(events);lost=0 if previous is None else previous-count
            rows.append({"lightmap":label,"category":"stepwise","stage":stage,"events":count,
                         "lost_from_previous":lost,
                         "fraction_of_previous":1. if previous is None else count/max(previous,1),
                         "fraction_of_input":count/max(input_count,1)})
            previous=count
    return pd.DataFrame(rows),common


def plot_event_cutflow(cutflow):
    """Plot event counts for every cut-flow stage and light map."""
    cutflow=cutflow[cutflow.category=="stepwise"]
    table=cutflow.pivot(index="stage",columns="lightmap",values="events")
    order=list(dict.fromkeys(cutflow.stage))
    axis=table.reindex(order).plot(kind="bar",figsize=(15,6))
    axis.set(ylabel="events",xlabel="",title="GlobalFit event cut flow")
    axis.tick_params(axis="x",rotation=35);axis.legend(title="lightmap",fontsize=8)
    axis.figure.tight_layout();return axis.figure


def plot_column_resolution(matched_by_study, bins=100):
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    for label, matched in matched_by_study.items():
        calibrated, scale = add_calibrated_residuals(matched)
        positive = calibrated.mc_energy > 0
        axes[0].scatter(calibrated.loc[positive, "mc_energy"], calibrated.loc[positive, "fitted_energy"],
                        s=2, alpha=.15, label=f"{label} ({scale:.3g} q/energy)")
        axes[1].hist(calibrated.energy_residual, bins=bins, histtype="step", density=True, label=label)
        fractional = calibrated.loc[positive, "energy_fractional_residual"].replace([np.inf, -np.inf], np.nan).dropna()
        axes[2].hist(fractional.clip(-5, 5), bins=bins, range=(-5, 5), histtype="step", density=True, label=label)
    limit = max([axis.get_xlim()[1] for axis in axes[:1]] + [1])
    axes[0].plot([0, limit], [0, limit], "k--", lw=1, label="ideal")
    axes[0].set(xlabel="MC energy per column", ylabel="calibrated fitted energy per column")
    axes[1].set(xlabel="fitted - MC energy", ylabel="density")
    axes[2].set(xlabel="(fitted - MC) / MC", ylabel="density")
    for axis in axes: axis.legend(fontsize=8)
    fig.tight_layout()
    return fig


def plot_column_charge_distributions(matched_by_study, bins=80,
                                     histogram_range=None, density=True,
                                     include_zero_truth=True):
    """Overlay all-column q_c reco and truth distributions per light map.

    Reconstructed charge is divided by the sample-wide charge/MC-energy
    calibration, so q_c_reco and q_c_true=mc_energy use the same units.
    The same bin edges are used in every panel.
    """
    calibrated = {}
    combined = []
    rows = []
    for label, matched in matched_by_study.items():
        frame, scale = add_calibrated_residuals(matched)
        if not include_zero_truth:
            frame = frame[frame.mc_energy > 0]
        calibrated[label] = frame
        combined.extend(frame.fitted_energy.replace([np.inf, -np.inf], np.nan).dropna())
        combined.extend(frame.mc_energy.replace([np.inf, -np.inf], np.nan).dropna())
        rows.append({"lightmap": label, "columns": len(frame),
                     "charge_per_mc_energy": scale,
                     "qc_reco_mean": frame.fitted_energy.mean(),
                     "qc_reco_std": frame.fitted_energy.std(),
                     "qc_true_mean": frame.mc_energy.mean(),
                     "qc_true_std": frame.mc_energy.std()})
    if histogram_range is None:
        finite = np.asarray(combined, float)
        finite = finite[np.isfinite(finite)]
        upper = np.quantile(finite, .995) if len(finite) else 1.
        histogram_range = (0., max(float(upper), 1e-12))
    edges = np.linspace(histogram_range[0], histogram_range[1], int(bins)+1)
    labels = list(calibrated)
    columns = min(2, max(1, len(labels)))
    rows_of_axes = int(np.ceil(len(labels)/columns))
    figure, axes = plt.subplots(rows_of_axes, columns,
                                figsize=(7*columns, 4.5*rows_of_axes),
                                squeeze=False, sharex=True)
    for axis, label in zip(axes.flat, labels):
        frame = calibrated[label]
        axis.hist(frame.mc_energy, bins=edges, histtype="step", density=density,
                  color="black", linewidth=2, label=r"$q_c^{true}$")
        axis.hist(frame.fitted_energy, bins=edges, histtype="step", density=density,
                  linewidth=1.8, label=r"$q_c^{reco}$")
        axis.set(title=label, xlabel=r"column charge $q_c$ [MC-energy units]",
                 ylabel="density" if density else "columns")
        axis.legend()
    for axis in axes.flat[len(labels):]:
        axis.set_visible(False)
    figure.tight_layout()
    return pd.DataFrame(rows).set_index("lightmap"), figure


def gaussian_peak_fit(values,histogram_range=(-2,2),bins=160,window_sigma=1.5,
                      fixed_fit_range=None,minimum_peak_fraction=None):
    """Fit ROOT's ``gaus`` to a TH1D, restricted to the requested x range."""
    import ROOT
    import uuid
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    if len(values) < 10:
        return np.nan,np.nan,0,np.nan,np.nan,np.nan
    counts,edges=np.histogram(values,bins=bins,range=histogram_range)
    centres=.5*(edges[:-1]+edges[1:])
    if fixed_fit_range is None:
        peak=int(np.argmax(counts))
    else:
        candidate=np.flatnonzero((centres>=fixed_fit_range[0])&(centres<=fixed_fit_range[1]))
        if not len(candidate):raise ValueError("fixed_fit_range contains no histogram bins")
        peak=int(candidate[np.argmax(counts[candidate])])
    mean0=centres[peak]
    half=.5*counts[peak];left=peak;right=peak
    while left>0 and counts[left]>half:left-=1
    while right<len(counts)-1 and counts[right]>half:right+=1
    sigma0=max((centres[right]-centres[left])/2.355,(edges[1]-edges[0])*2)
    if fixed_fit_range is None:
        fit_low=mean0-window_sigma*sigma0;fit_high=mean0+window_sigma*sigma0
    else:
        fit_low,fit_high=map(float,fixed_fit_range)
        if not fit_low<fit_high:raise ValueError("fixed_fit_range must have low < high")
    if minimum_peak_fraction is not None:
        if not 0<minimum_peak_fraction<1:raise ValueError("minimum_peak_fraction must be between 0 and 1")
        threshold=minimum_peak_fraction*counts[peak];core_left=peak;core_right=peak
        while core_left>0 and counts[core_left-1]>=threshold:core_left-=1
        while core_right<len(counts)-1 and counts[core_right+1]>=threshold:core_right+=1
        fit_low=max(fit_low,edges[core_left]);fit_high=min(fit_high,edges[core_right+1])
    unique=uuid.uuid4().hex
    histogram=ROOT.TH1D(f"gaussian_data_{unique}","",int(bins),float(histogram_range[0]),float(histogram_range[1]))
    histogram.SetDirectory(0)
    for index,count in enumerate(counts,1):
        histogram.SetBinContent(index,float(count))
        histogram.SetBinError(index,float(np.sqrt(count)))
    histogram.SetEntries(float(len(values)))
    function=ROOT.TF1(f"gaussian_fit_{unique}","gaus",float(fit_low),float(fit_high))
    function.SetParameters(float(counts[peak]),float(mean0),float(sigma0))
    function.SetParLimits(0,0.,max(float(counts[peak])*10.,1.))
    function.SetParLimits(1,float(fit_low),float(fit_high))
    function.SetParLimits(2,float(edges[1]-edges[0])*.2,max(float(fit_high-fit_low),float(edges[1]-edges[0])))
    try:
        result=histogram.Fit(function,"RQ0S")
        if int(result)!=0:raise RuntimeError(f"ROOT Gaussian fit failed with status {int(result)}")
        amplitude,mean,sigma=(function.GetParameter(i) for i in range(3))
    except Exception:
        selected=values[(values>=fit_low)&(values<=fit_high)];mean=selected.mean();sigma=selected.std(ddof=1);amplitude=counts[peak]
    selected=(values>=fit_low)&(values<=fit_high)
    inside=((values>=histogram_range[0])&(values<=histogram_range[1])).sum()
    amplitude_density=amplitude/(max(inside,1)*(edges[1]-edges[0]))
    return float(mean),float(abs(sigma)),int(selected.sum()),float(fit_low),float(fit_high),float(amplitude_density)


def gaussian_core_fit(values, clip_sigma=1.5, iterations=0):
    """Backward-compatible local peak fit used by aggregation tables."""
    del iterations
    mean,sigma,entries,_,_,_=gaussian_peak_fit(values,window_sigma=clip_sigma)
    return mean,sigma,entries


def plot_relative_energy_resolution(matched_by_study, display_range=(-2, 2),
                                    gaussian_fit_range=(-0.4,0.4),
                                    bins=120):
    """Fit a pure Gaussian using exactly the requested x range."""
    figures = {}
    rows = []
    for label, matched in matched_by_study.items():
        figure,axis=plt.subplots(figsize=(9,6));figures[label]=figure
        calibrated, scale = add_calibrated_residuals(matched)
        values = calibrated.loc[calibrated.mc_energy > 0, "energy_fractional_residual"].to_numpy(float)
        values = values[np.isfinite(values)]
        axis.hist(values, bins=bins, range=display_range, histtype="step", density=True, label=f"{label} data")
        mean,sigma,entries,fit_low,fit_high,amplitude=gaussian_peak_fit(
            values,histogram_range=display_range,fixed_fit_range=gaussian_fit_range)
        if np.isfinite(sigma) and sigma > 0:
            x=np.linspace(fit_low,fit_high,400)
            axis.plot(x,amplitude*np.exp(-.5*((x-mean)/sigma)**2),"--",lw=2,
                      label=f"Gaussian fit: bias={mean:.3g}, sigma/E={sigma:.3g}")
            axis.axvspan(fit_low,fit_high,color="grey",alpha=.08,label=f"fit range [{fit_low:g}, {fit_high:g}]")
        rows.append({"lightmap": label, "charge_per_mc_energy": scale,
                     "gaussian_bias": mean, "sigma_over_E": sigma,
                     "fit_low":fit_low,"fit_high":fit_high,
                     "fit_entries": entries, "all_positive_mc_columns": len(values)})
        axis.set(xlabel="(calibrated fitted energy - MC energy) / MC energy",ylabel="density",
                 title=f"Relative column-energy resolution — {label}")
        axis.legend(fontsize=9);axis.set_xlim(display_range);figure.tight_layout()
    return pd.DataFrame(rows).set_index("lightmap"), figures


def truncated_dedx_table(matched_by_study, discard_fraction=0.30):
    """Event-level lower-fraction mean dE/dx for fitted and MC columns."""
    if not 0 <= discard_fraction < 1:
        raise ValueError("discard_fraction must be in [0, 1)")
    rows = []
    for label, matched in matched_by_study.items():
        calibrated, scale = add_calibrated_residuals(matched)
        for event, columns in calibrated.groupby("event"):
            valid = columns.path_length > 0
            values = columns.loc[valid].copy()
            values["fit_dedx"] = values.fitted_energy / values.path_length
            values["mc_dedx"] = values.mc_energy / values.path_length
            keep = max(1, int(np.floor((1.0 - discard_fraction) * len(values))))
            fit_kept = np.sort(values.fit_dedx.to_numpy(float))[:keep]
            mc_kept = np.sort(values.mc_dedx.to_numpy(float))[:keep]
            rows.append({"lightmap": label, "event": int(event), "columns": len(values),
                         "kept_columns": keep, "discard_fraction": discard_fraction,
                         "fit_truncated_dedx": float(fit_kept.mean()),
                         "mc_truncated_dedx": float(mc_kept.mean()),
                         "charge_per_mc_energy": scale})
    return pd.DataFrame(rows)


def plot_truncated_dedx(table):
    figure, axes = plt.subplots(1, 2, figsize=(14, 5))
    truth_drawn = False
    for label, frame in table.groupby("lightmap", sort=False):
        axes[0].hist(frame.fit_truncated_dedx, bins=80, histtype="step", density=True, label=f"{label} fit")
        if not truth_drawn:
            axes[0].hist(frame.mc_truncated_dedx, bins=80, histtype="step", density=True,
                         color="black", lw=2, label="MC truth")
            truth_drawn = True
        residual = (frame.fit_truncated_dedx-frame.mc_truncated_dedx) / frame.mc_truncated_dedx.replace(0, np.nan)
        axes[1].hist(residual.replace([np.inf,-np.inf],np.nan).dropna(), bins=100, range=(-2,2),
                     histtype="step", density=True, label=label)
    axes[0].set(xlabel="70% truncated mean dE/dx", ylabel="events")
    axes[1].set(xlabel="relative truncated-mean residual", ylabel="events")
    for axis in axes: axis.legend(fontsize=8)
    figure.tight_layout();return figure


def plot_truncated_dedx_resolution(table, relative_display_range=(-1,1), relative_fit_range=(-1,1), dedx_fit_range=None,
                                   relative_bins=50, dedx_bins=80):
    """Plot truncated dE/dx and its relative residual in separate figures."""
    if relative_bins < 5 or dedx_bins < 5:raise ValueError("Histogram bin counts must be at least 5")
    all_dedx=np.concatenate([table.fit_truncated_dedx.to_numpy(),table.mc_truncated_dedx.to_numpy()])
    all_dedx=all_dedx[np.isfinite(all_dedx)];dedx_range=tuple(np.quantile(all_dedx,[.005,.995]))
    if dedx_fit_range is None:dedx_fit_range=dedx_range
    dedx_figure,dedx_axis=plt.subplots(figsize=(10,6));residual_figure,residual_axis=plt.subplots(figsize=(10,6));rows=[]
    mc_drawn=False
    for label,frame in table.groupby("lightmap",sort=False):
        fit_values=frame.fit_truncated_dedx.to_numpy(float);mc_values=frame.mc_truncated_dedx.to_numpy(float)
        fit_mean,fit_sigma,fit_entries,fit_low,fit_high,fit_amplitude=gaussian_peak_fit(fit_values,histogram_range=dedx_range,fixed_fit_range=dedx_fit_range)
        mc_mean,mc_sigma,mc_entries,mc_low,mc_high,mc_amplitude=gaussian_peak_fit(mc_values,histogram_range=dedx_range,fixed_fit_range=dedx_fit_range)
        dedx_axis.hist(fit_values,bins=dedx_bins,range=dedx_range,histtype="step",density=True,label=f"{label} data")
        xd=np.linspace(fit_low,fit_high,400)
        dedx_axis.plot(xd,fit_amplitude*np.exp(-.5*((xd-fit_mean)/fit_sigma)**2),"--",lw=1.5,
                       label=f"{label} fit: sigma/mean={fit_sigma/fit_mean:.3g}")
        if not mc_drawn:
            dedx_axis.hist(mc_values,bins=dedx_bins,range=dedx_range,histtype="step",density=True,color="black",lw=1.5,label="MC truth data")
            xm=np.linspace(mc_low,mc_high,400)
            dedx_axis.plot(xm,mc_amplitude*np.exp(-.5*((xm-mc_mean)/mc_sigma)**2),"k--",lw=1.5,label=f"MC fit: sigma/mean={mc_sigma/mc_mean:.3g}")
            mc_drawn=True
        denominator=frame.mc_truncated_dedx.replace(0,np.nan)
        values=((frame.fit_truncated_dedx-frame.mc_truncated_dedx)/denominator).replace([np.inf,-np.inf],np.nan).dropna().to_numpy()
        residual_axis.hist(values,bins=relative_bins,range=relative_display_range,histtype="step",density=True,label=f"{label} data")
        mean,sigma,entries,fit_low,fit_high,residual_amplitude=gaussian_peak_fit(values,histogram_range=relative_display_range,bins=relative_bins,fixed_fit_range=relative_fit_range)
        x=np.linspace(fit_low,fit_high,400)
        residual_axis.plot(x,residual_amplitude*np.exp(-.5*((x-mean)/sigma)**2),"--",label=f"{label}: bias={mean:.3g}, sigma/E={sigma:.3g}")
        rows.append({"lightmap":label,"fit_peak_mean":fit_mean,"fit_peak_sigma":fit_sigma,
                     "fit_sigma_over_mean":fit_sigma/fit_mean,"mc_peak_mean":mc_mean,"mc_peak_sigma":mc_sigma,
                     "mc_sigma_over_mean":mc_sigma/mc_mean,"relative_residual_bias":mean,
                     "relative_residual_sigma":sigma,"residual_fit_low":fit_low,"residual_fit_high":fit_high,
                     "residual_fit_entries":entries,"events":len(values)})
    dedx_axis.set(xlabel="70% truncated mean dE/dx",ylabel="density",title="Truncated-mean dE/dx distributions and peak fits")
    residual_axis.set(xlabel="(fitted truncated dE/dx - MC truncated dE/dx) / MC truncated dE/dx",ylabel="density",title="Relative truncated-mean dE/dx resolution")
    residual_axis.set_xlim(relative_display_range)
    dedx_axis.legend(fontsize=8,ncol=2);residual_axis.legend(fontsize=8,ncol=2)
    dedx_figure.tight_layout();residual_figure.tight_layout()
    return pd.DataFrame(rows).set_index("lightmap"),(dedx_figure,residual_figure)


def neighbour_migration_table(matched_by_study):
    """Adjacent-column calibrated residual pairs for migration studies."""
    rows = []
    for label, matched in matched_by_study.items():
        calibrated, _ = add_calibrated_residuals(matched)
        for event, columns in calibrated.groupby("event"):
            columns = columns.sort_values("column_index")
            values = columns[["column_index", "mc_energy", "fitted_energy", "energy_residual"]].to_numpy(float)
            for left, right in zip(values[:-1], values[1:]):
                if int(right[0]) != int(left[0]) + 1: continue
                rows.append({"lightmap": label, "event": int(event), "column_index": int(left[0]),
                             "mc_energy_k": left[1], "residual_k": left[3],
                             "mc_energy_next": right[1], "residual_next": right[3],
                             "pair_mc_energy": left[1]+right[1],
                             "pair_fitted_energy": left[2]+right[2],
                             "pair_residual": left[3]+right[3]})
    return pd.DataFrame(rows)


def plot_neighbour_migration(table, limit_quantile=0.98):
    labels = list(table.lightmap.unique());figure, axes = plt.subplots(
        int(np.ceil(len(labels)/2)), 2, figsize=(13, 5.5*np.ceil(len(labels)/2)), squeeze=False)
    for axis, label in zip(axes.flat, labels):
        frame = table[table.lightmap == label]
        limit = np.quantile(np.abs(frame[["residual_k","residual_next"]].to_numpy()), limit_quantile)
        correlation = frame.residual_k.corr(frame.residual_next)
        image = axis.hexbin(frame.residual_k, frame.residual_next, gridsize=70,
                            extent=(-limit,limit,-limit,limit), bins="log", mincnt=1, cmap="viridis")
        axis.axhline(0,color="grey",lw=.8);axis.axvline(0,color="grey",lw=.8)
        axis.set(xlabel=r"$r_k=E_k^{fit}-E_k^{MC}$", ylabel=r"$r_{k+1}$",
                 title=f"{label}: adjacent residual correlation = {correlation:.3f}")
        figure.colorbar(image, ax=axis, label="log10(entries)")
    for axis in axes.flat[len(labels):]: axis.set_visible(False)
    figure.tight_layout();return figure


def aggregation_resolution_table(matched_by_study, widths=(1,2,3)):
    """Gaussian relative resolution for sliding sums of adjacent columns."""
    rows=[]
    for label, matched in matched_by_study.items():
        calibrated, _ = add_calibrated_residuals(matched)
        for width in widths:
            residuals=[]
            for _, columns in calibrated.groupby("event"):
                columns=columns.sort_values("column_index")
                index=columns.column_index.to_numpy(int);fit=columns.fitted_energy.to_numpy(float);truth=columns.mc_energy.to_numpy(float)
                for start in range(max(0,len(columns)-width+1)):
                    if index[start+width-1]-index[start] != width-1: continue
                    mc=truth[start:start+width].sum()
                    if mc>0: residuals.append((fit[start:start+width].sum()-mc)/mc)
            mean,sigma,entries=gaussian_core_fit(residuals)
            rows.append({"lightmap":label,"aggregation_columns":width,"aggregation_mm":10*width,
                         "gaussian_bias":mean,"sigma_over_E":sigma,"fit_entries":entries,"all_windows":len(residuals)})
    return pd.DataFrame(rows)


def plot_aggregation_resolution(table):
    figure,axes=plt.subplots(1,2,figsize=(12,5))
    for label,frame in table.groupby("lightmap",sort=False):
        axes[0].plot(frame.aggregation_mm,frame.sigma_over_E,"o-",label=label)
        axes[1].plot(frame.aggregation_mm,frame.gaussian_bias,"o-",label=label)
    axes[0].set(xlabel="longitudinal aggregation width [mm]",ylabel=r"Gaussian $\sigma/E$")
    axes[1].set(xlabel="longitudinal aggregation width [mm]",ylabel="Gaussian relative bias")
    for axis in axes: axis.legend(fontsize=8);axis.set_xticks(sorted(table.aggregation_mm.unique()))
    figure.tight_layout();return figure


def profile_shift_table(matched_by_study, maximum_lag=5):
    """Best signed fitted-to-MC profile lag and cumulative-shape distance."""
    rows=[]
    for label,matched in matched_by_study.items():
        calibrated,_=add_calibrated_residuals(matched)
        for event,columns in calibrated.groupby("event"):
            columns=columns.sort_values("column_index");lo=int(columns.column_index.min());hi=int(columns.column_index.max())
            index=np.arange(lo,hi+1);fit=np.zeros(len(index));truth=np.zeros(len(index))
            position=columns.column_index.to_numpy(int)-lo
            fit[position]=columns.fitted_energy.to_numpy(float);truth[position]=columns.mc_energy.to_numpy(float)
            if fit.sum()<=0 or truth.sum()<=0: continue
            fit/=fit.sum();truth/=truth.sum();scores=[]
            for lag in range(-maximum_lag,maximum_lag+1):
                shifted=np.zeros_like(fit)
                if lag>=0: shifted[lag:]=fit[:len(fit)-lag] if lag else fit
                else: shifted[:lag]=fit[-lag:]
                denominator=np.linalg.norm(shifted)*np.linalg.norm(truth)
                scores.append(np.dot(shifted,truth)/denominator if denominator else -np.inf)
            best_lag=int(np.arange(-maximum_lag,maximum_lag+1)[int(np.argmax(scores))])
            cumulative_distance=float(np.max(np.abs(np.cumsum(fit)-np.cumsum(truth))))
            rows.append({"lightmap":label,"event":int(event),"best_lag_columns":best_lag,
                         "best_lag_mm":10*best_lag,"correlation":max(scores),
                         "maximum_cumulative_difference":cumulative_distance})
    return pd.DataFrame(rows)


def plot_profile_shifts(table, maximum_lag=5):
    figure,axes=plt.subplots(1,2,figsize=(13,5))
    bins=np.arange(-maximum_lag-.5,maximum_lag+1.5)
    for label,frame in table.groupby("lightmap",sort=False):
        axes[0].hist(frame.best_lag_columns,bins=bins,histtype="step",density=True,label=label)
        axes[1].hist(frame.maximum_cumulative_difference,bins=70,histtype="step",density=True,label=label)
    axes[0].set(xlabel="best fitted-to-MC signed lag [columns]",ylabel="events",xticks=range(-maximum_lag,maximum_lag+1))
    axes[1].set(xlabel="max |cumulative fitted fraction - cumulative MC fraction|",ylabel="events")
    for axis in axes: axis.legend(fontsize=8)
    figure.tight_layout();return figure


def plot_event_cumulative_profiles(matched_by_study,event):
    figure,axis=plt.subplots(figsize=(12,5));truth_drawn=False
    for label,matched in matched_by_study.items():
        calibrated,_=add_calibrated_residuals(matched);frame=calibrated[calibrated.event==event].sort_values("column_index")
        if frame.empty: continue
        fit=frame.fitted_energy.to_numpy(float);truth=frame.mc_energy.to_numpy(float)
        fit=fit/fit.sum() if fit.sum() else fit;truth=truth/truth.sum() if truth.sum() else truth
        coordinate=.5*(frame.column_low+frame.column_high)
        axis.plot(coordinate,np.cumsum(fit),label=f"{label} fit")
        if not truth_drawn:axis.step(coordinate,np.cumsum(truth),where="mid",color="black",lw=2,label="MC truth");truth_drawn=True
    axis.set(xlabel="physical dominant-axis column centre [mm]",ylabel="cumulative profile fraction",title=f"Event {event}")
    axis.legend(fontsize=8,ncol=2);figure.tight_layout();return figure


def _line_residuals(points, point, direction):
    direction = direction / np.linalg.norm(direction)
    delta = points - point
    return delta - np.outer(delta @ direction, direction)


def muon_cube_start_positions(segments,track_points,maximum_distance_mm=2.0):
    """Recover one segment start per cube along the primary PDG-13 trajectory.

    Existing flat trees do not store contributor track IDs. Match every cube
    start to the dominant primary-muon trajectory polyline and retain the
    closest start in each cube when it is within maximum_distance_mm.
    """
    candidates=track_points[(np.abs(track_points.pdg)==13)&(track_points.parent_id==0)]
    rows=[]
    detector_segments=segments[segments.detector==0]
    if "contributor_pdgs" in detector_segments:
        has_muon=np.asarray([any(abs(int(pdg))==13 for pdg in pdgs) for pdgs in detector_segments.contributor_pdgs])
        detector_segments=detector_segments[has_muon]
    segment_groups={int(event):frame for event,frame in detector_segments.groupby("event",sort=False)}
    for event,event_tracks in candidates.groupby("event",sort=False):
        counts=event_tracks.groupby("track_id").size()
        if counts.empty or int(event) not in segment_groups:continue
        trajectory=event_tracks[event_tracks.track_id==counts.idxmax()].sort_values("point")[["x","y","z"]].to_numpy(float)
        if len(trajectory)<2:continue
        starts=segment_groups[int(event)].copy();points=starts[["start_x","start_y","start_z"]].to_numpy(float)
        begin=trajectory[:-1];direction=trajectory[1:]-begin;length2=np.sum(direction*direction,axis=1)
        distance=np.full(len(points),np.inf)
        for index,point in enumerate(points):
            projection=np.divide(np.sum((point-begin)*direction,axis=1),length2,
                                 out=np.zeros_like(length2),where=length2>0)
            projection=np.clip(projection,0,1);closest=begin+projection[:,None]*direction
            distance[index]=np.sqrt(np.min(np.sum((closest-point)**2,axis=1)))
        starts["muon_trajectory_distance"]=distance
        starts=starts[starts.muon_trajectory_distance<=maximum_distance_mm]
        starts=starts.sort_values("muon_trajectory_distance").drop_duplicates(["cube_x","cube_y","cube_z"])
        rows.append(starts)
    return pd.concat(rows,ignore_index=True) if rows else segments.iloc[0:0].copy()


def segment_line_residuals(studies,segments,track_points,maximum_muon_match_distance_mm=2.0):
    """PDG-13 cube-start residuals to fitted and pure-MC PCA lines."""
    primary=muon_cube_start_positions(segments,track_points,maximum_muon_match_distance_mm)
    output = {}
    mc_rows = []
    for event, values in primary.groupby("event"):
        points = values[["start_x", "start_y", "start_z"]].to_numpy()
        if len(points) < 3: continue
        centre = points.mean(axis=0)
        _, _, vh = np.linalg.svd(points - centre, full_matrices=False)
        residual = _line_residuals(points, centre, vh[0])
        for vector in residual: mc_rows.append((int(event), *vector, np.linalg.norm(vector)))
    output["MC straight-line reference"] = pd.DataFrame(mc_rows, columns=["event", "rx", "ry", "rz", "distance"])
    for label, study in studies.items():
        rows = []
        indexed = study["fit"].set_index("event")
        for event, values in primary.groupby("event"):
            if event not in indexed.index: continue
            fit = indexed.loc[event]
            points = values[["start_x", "start_y", "start_z"]].to_numpy()
            residual = _line_residuals(points,
                np.array([fit.fit_x, fit.fit_y, fit.fit_z]),
                np.array([fit.fit_dx, fit.fit_dy, fit.fit_dz]))
            for vector in residual: rows.append((int(event), *vector, np.linalg.norm(vector)))
        output[label] = pd.DataFrame(rows, columns=["event", "rx", "ry", "rz", "distance"])
    return output


def plot_line_residuals(residuals, distance_max=20, signed_fit_range=(-1,1), bins=100):
    """Plot residuals and return ROOT-Gaussian signed-component resolutions.

    The non-negative perpendicular distance is not Gaussian; its 68% quantile
    is reported instead of assigning it a misleading Gaussian sigma.
    """
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fields = [("rx", "X residual [mm]"), ("ry", "Y residual [mm]"),
              ("rz", "Z residual [mm]"), ("distance", "perpendicular distance [mm]")]
    rows=[]
    for label, frame in residuals.items():
        for axis, (field, title) in zip(axes.flat, fields):
            values = frame[field].replace([np.inf,-np.inf],np.nan).dropna().to_numpy(float)
            limits = (0, distance_max) if field == "distance" else (-distance_max, distance_max)
            colour="black" if label.startswith("MC ") else axis._get_lines.get_next_color()
            linewidth=2.0 if label.startswith("MC ") else 1.5
            histogram_label=label
            if field != "distance":
                mean,sigma,entries,fit_low,fit_high,amplitude=gaussian_peak_fit(
                    values,histogram_range=limits,bins=bins,fixed_fit_range=signed_fit_range)
                histogram_label=f"{label}: sigma={sigma:.3g} mm"
                if np.isfinite(sigma) and sigma>0:
                    x=np.linspace(fit_low,fit_high,400)
                    axis.plot(x,amplitude*np.exp(-.5*((x-mean)/sigma)**2),"--",lw=1.4,color=colour)
                rows.append({"sample":label,"component":field,"gaussian_mean_mm":mean,
                             "gaussian_sigma_mm":sigma,"fit_low_mm":fit_low,
                             "fit_high_mm":fit_high,"fit_entries":entries})
            axis.hist(values,bins=bins,range=limits,histtype="step",density=True,
                      label=histogram_label,color=colour,linewidth=linewidth)
            axis.set(xlabel=title, ylabel="density")
        distance=frame.distance.replace([np.inf,-np.inf],np.nan).dropna().to_numpy(float)
        component_rows=[row for row in rows if row["sample"]==label]
        sigmas=np.asarray([row["gaussian_sigma_mm"] for row in component_rows],float)
        rows.append({"sample":label,"component":"point_resolution_summary",
                     "gaussian_mean_mm":np.nan,
                     "gaussian_sigma_mm":np.sqrt(np.nansum(sigmas**2)/2.),
                     "fit_low_mm":np.nan,"fit_high_mm":np.nan,"fit_entries":len(distance),
                     "distance_68_mm":np.quantile(distance,.68) if len(distance) else np.nan})
    for axis in axes.flat: axis.legend(fontsize=8)
    fig.tight_layout()
    return pd.DataFrame(rows),fig


def filter_line_residuals(residuals, studies, event="ALL", convergence="column", common_events=True):
    """Select all/single-event residuals and optionally converged fits.

    With common_events=True, every light map and the MC reference use the same
    intersection of accepted events, which makes normalized overlays fair.
    """
    if convergence not in (None, "none", "column", "geometry", "both"):
        raise ValueError("convergence must be None, 'none', 'column', 'geometry', or 'both'")
    accepted = {}
    for label, study in studies.items():
        fit = study["fit"]
        mask = np.ones(len(fit), dtype=bool)
        if convergence in ("column", "both"):
            if "column_fit_converged" not in fit:
                raise KeyError(f"{label} has no column_fit_converged branch; regenerate this fit")
            mask &= fit.column_fit_converged.to_numpy() == 1
        if convergence in ("geometry", "both"):
            mask &= fit.status.to_numpy() == 0
        accepted[label] = set(fit.loc[mask, "event"].astype(int))
    if str(event).upper() != "ALL":
        wanted = {int(event)}
        accepted = {label: events & wanted for label, events in accepted.items()}
    if common_events:
        common = set.intersection(*accepted.values()) if accepted else set()
        selected = {label: common for label in accepted}
        mc_events = common
    else:
        selected = accepted
        mc_events = set.union(*accepted.values()) if accepted else set()
    output = {}
    for label, frame in residuals.items():
        events = mc_events if label == "MC straight-line reference" else selected.get(label, set())
        output[label] = frame[frame.event.isin(events)].copy()
    return output, accepted, mc_events


def plot_event_profile(studies, matched_by_study, event, labels=None, normalize=False):
    """Plot one event's fitted and MC column-charge profiles without widgets.

    This is also the widget-independent fallback for notebook frontends whose
    JavaScript widget manager is unavailable or stale.
    """
    available = list(studies)
    labels = available if labels is None else list(labels)
    unknown = [label for label in labels if label not in studies]
    if unknown:
        raise KeyError(f"Unknown light-map labels: {unknown}; available: {available}")
    fig, axis = plt.subplots(figsize=(13, 5))
    truth_drawn = False
    for label in labels:
        frame = matched_by_study[label]
        frame = frame[frame.event == int(event)].sort_values("column_index")
        if frame.empty:
            continue
        fitted = frame.fitted_selected_fibre_charge.to_numpy(float)
        truth = frame.mc_energy.to_numpy(float)
        if normalize:
            fitted = fitted / fitted.sum() if fitted.sum() else fitted
            truth = truth / truth.sum() if truth.sum() else truth
        else:
            fitted = fitted / charge_calibration(matched_by_study[label])
        coordinate = .5 * (frame.column_low + frame.column_high)
        axis.plot(coordinate, fitted, marker=".", label=f"{label} fit")
        if not truth_drawn:
            axis.step(coordinate, truth, where="mid", color="black", lw=2, label="MC truth")
            truth_drawn = True
    if not axis.lines:
        raise ValueError(f"Event {event} is absent from the selected light maps")
    axis.set(xlabel="physical dominant-axis column centre [mm]",
             ylabel="profile fraction" if normalize else "MC energy / calibrated fitted energy",
             title=f"Event {int(event)}: fitted qk versus MC truth")
    axis.legend(fontsize=8, ncol=2)
    axis.grid(alpha=.2)
    fig.tight_layout()
    return fig, axis


def interactive_event(studies, matched_by_study):
    import ipywidgets as widgets
    from IPython.display import display
    labels = list(studies)
    events = sorted(set.intersection(*[set(s["fit"].event.astype(int)) for s in studies.values()]))
    event_widget = widgets.IntSlider(value=events[0], min=min(events), max=max(events), step=1,
                                     description="Event", continuous_update=False,
                                     layout=widgets.Layout(width="650px"))
    maps_widget = widgets.SelectMultiple(options=labels, value=tuple(labels), description="Light maps",
                                          layout=widgets.Layout(width="450px", height="120px"))
    normalize_widget = widgets.Checkbox(value=False, description="Normalize profiles")
    output = widgets.Output()

    def draw(*_):
        with output:
            output.clear_output(wait=True)
            fig, _ = plot_event_profile(studies, matched_by_study, event_widget.value,
                                        labels=maps_widget.value,
                                        normalize=normalize_widget.value)
            plt.show()

    event_widget.observe(draw, names="value");maps_widget.observe(draw, names="value");normalize_widget.observe(draw, names="value")
    display(widgets.VBox([event_widget, maps_widget, normalize_widget]), output);draw()
    return output
