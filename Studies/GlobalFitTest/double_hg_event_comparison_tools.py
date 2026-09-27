"""Side-by-side legacy/double-HG GlobalFit displays."""
from pathlib import Path
import re
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

VIEWS={2:("XY (Z fibres)","x","y"),0:("XZ (Y fibres)","x","z"),1:("YZ (X fibres)","y","z")}

def _disable_capture():
    try:
        from IPython import get_ipython
        from JupyROOT.helpers import utils
        shell=get_ipython()
        for capture in list(utils.captures):
            for event in ("pre_execute","post_execute"):
                callback=getattr(capture,event,None)
                if shell and callback:
                    try: shell.events.unregister(event,callback)
                    except ValueError: pass
        utils.captures.clear()
    except (ImportError,AttributeError): pass

def tree_frame(filename,tree,branches):
    # ROOT 6.32 tries `%gui ROOT`, which newer IPython rejects because its
    # accepted GUI name list no longer contains uppercase ROOT.  This notebook
    # draws only with Matplotlib, so ROOT's GUI event thread is unnecessary.
    import ROOT
    ROOT.PyConfig.StartGuiThread=False
    _disable_capture()
    source=ROOT.TFile.Open(str(filename),"READ")
    if not source or source.IsZombie(): raise OSError(f"Cannot open {filename}")
    obj=source.Get(tree)
    if not obj: source.Close(); raise KeyError(f"Missing {tree} in {filename}")
    available={b.GetName() for b in obj.GetListOfBranches()}; selected=[b for b in branches if b in available]
    arrays=ROOT.RDataFrame(obj).AsNumpy(selected); source.Close()
    return pd.DataFrame({name:np.asarray(arrays[name]) for name in selected})

def discover_pairs(base):
    legacy,dhg={},{}
    for path in Path(base).expanduser().resolve().glob("*/global_fit.root"):
        name=path.parent.name
        if "double_hg" in name: match=re.search(r"double_hg_([0-9]+p[0-9]+)mm",name); target=dhg
        elif name.startswith("250514_"): match=re.search(r"_([0-9]+)mm_",name); target=legacy
        else: match=re.search(r"_([0-9]+p[0-9]+)mm_",name); target=legacy
        if match: target[float(match.group(1).replace("p","."))]=path
    return {s:{"Legacy":legacy[s],"Double-HG":dhg[s]} for s in sorted(set(legacy)&set(dhg))}

def load_pair(base,scatter_mm):
    pairs=discover_pairs(base); scatter_mm=float(scatter_mm)
    if scatter_mm not in pairs: raise KeyError(f"No pair for {scatter_mm:g} mm; available: {list(pairs)}")
    runs={}
    for label,path in pairs[scatter_mm].items():
        print(f"Loading {label}: {path}")
        runs[label]={"scatter":scatter_mm,"path":path,
          "hits":tree_frame(path,"fiber_hits",["event","projection","x","y","z","charge"]),
          "selected":tree_frame(path,"global_fit_fibres",["event","projection","x","y","z","charge"]),
          "fit":tree_frame(path,"global_fit",["event","status","edm","chi2","ndof","minimum_charge","seed_dx","seed_dy","seed_dz","fit_x","fit_y","fit_z","fit_dx","fit_dy","fit_dz","observations_before_clustering","observations_after_dbscan","observations"])}
    truth=tree_frame(runs["Legacy"]["path"],"mc_virtual_segments",["event","detector","primary_id","primary_pdg","segment","start_x","start_y","start_z","stop_x","stop_y","stop_z"])
    return runs,truth

def discover_reconstruction_pairs(base):
    """Runs containing both the previous and track-aware ROOT outputs."""
    found={}
    for directory in Path(base).expanduser().resolve().iterdir():
        old=directory/"global_fit.root"; new=directory/"global_fit_columns_track_aware.root"
        if not (directory.is_dir() and old.exists() and new.exists()): continue
        name=directory.name; model="Double-HG" if "double_hg" in name else "Legacy lightmap"
        pattern=r"double_hg_([0-9]+p[0-9]+)mm" if "double_hg" in name else r"_([0-9]+p[0-9]+)mm_"
        match=re.search(pattern,name)
        if match: scatter=float(match.group(1).replace("p","."))
        elif "_1mm_" in name: scatter=1.0
        else: continue
        # Prefer the consistently named non-dated legacy branch if both exist.
        key=(model,scatter)
        if key not in found or not name.startswith("250514_"): found[key]=directory
    return found

def load_reconstruction_pair(directory):
    directory=Path(directory); runs={}
    for label,filename in (("Previous fit","global_fit.root"),("Track-aware fit","global_fit_columns_track_aware.root")):
        path=directory/filename; print(f"Loading {label}: {path}")
        runs[label]={"scatter":np.nan,"path":path,
          "hits":tree_frame(path,"fiber_hits",["event","projection","x","y","z","charge"]),
          "selected":tree_frame(path,"global_fit_fibres",["event","projection","x","y","z","charge"]),
          "fit":tree_frame(path,"global_fit",["event","status","edm","chi2","ndof","minimum_charge","seed_dx","seed_dy","seed_dz","fit_x","fit_y","fit_z","fit_dx","fit_dy","fit_dz","observations_before_clustering","observations_after_dbscan","observations","track_clusters_found","track_clusters_merged","track_branches_rejected"]),
          "mc_line":tree_frame(path,"mc_segment_line_fit",["event","fit_x","fit_y","fit_z","fit_dx","fit_dy","fit_dz","rms_residual_mm"])}
    truth=tree_frame(runs["Previous fit"]["path"],"mc_virtual_segments",["event","detector","primary_id","primary_pdg","segment","start_x","start_y","start_z","stop_x","stop_y","stop_z"])
    return runs,truth

def direction_resolution_table(runs,accepted_statuses=None):
    rows=[]
    for label,run in runs.items():
        fit=run["fit"]
        if accepted_statuses is not None: fit=fit[fit.status.isin(accepted_statuses)]
        merged=fit.merge(run["mc_line"],on="event",suffixes=("_reco","_mc"))
        reco=merged[["fit_dx_reco","fit_dy_reco","fit_dz_reco"]].to_numpy(float)
        truth=merged[["fit_dx_mc","fit_dy_mc","fit_dz_mc"]].to_numpy(float)
        reco/=np.linalg.norm(reco,axis=1)[:,None]; truth/=np.linalg.norm(truth,axis=1)[:,None]
        merged["direction_error_deg"]=np.degrees(np.arccos(np.clip(np.abs(np.sum(reco*truth,axis=1)),0,1)))
        merged["model"]=label; rows.append(merged[["event","model","status","direction_error_deg"]])
    return pd.concat(rows,ignore_index=True)

def plot_direction_resolution(runs,accepted_statuses=None,maximum_degrees=10,bins=100):
    data=direction_resolution_table(runs,accepted_statuses); fig,ax=plt.subplots(figsize=(9,5))
    rows=[]
    for label,frame in data.groupby("model",sort=False):
        values=frame.direction_error_deg.to_numpy(float); ax.hist(values,bins=bins,range=(0,maximum_degrees),histtype="step",density=True,lw=1.6,label=label)
        rows.append({"method":label,"events":len(values),"median_deg":np.median(values),"68pct_deg":np.quantile(values,.68),"90pct_deg":np.quantile(values,.90)})
    ax.set(xlabel="absolute reconstructed–MC direction angle [deg]",ylabel="density",title="Track-direction resolution"); ax.legend(); fig.tight_layout()
    return data,pd.DataFrame(rows),fig

def _weighted_group(frame,value,weight="charge"):
    rows=[]
    for (event,column),group in frame.groupby(["event","column_center_x"],sort=False):
        w=group[weight].to_numpy(float); v=group[value].to_numpy(float)
        rows.append((int(event),float(column),len(group),float(np.dot(w,v)/w.sum()) if w.sum()>0 else np.nan,float(w.sum())))
    return pd.DataFrame(rows,columns=["event","column_center_x","multiplicity",f"coc_{value}","weight_sum"])

def _primary_crossings(track_points,column_centres):
    rows=[]; muons=track_points[(track_points.parent_id==0)&(np.abs(track_points.pdg)==13)]
    for event,event_tracks in muons.groupby("event",sort=False):
        counts=event_tracks.groupby("track_id").size(); trajectory=event_tracks[event_tracks.track_id==counts.idxmax()].sort_values("point")
        xyz=trajectory[["x","y","z"]].to_numpy(float)
        for xc in sorted(column_centres.get(int(event),())):
            candidates=[]
            for a,b in zip(xyz[:-1],xyz[1:]):
                if abs(b[0]-a[0])<1e-9: continue
                fraction=(xc-a[0])/(b[0]-a[0])
                if 0<=fraction<=1: candidates.append(a+fraction*(b-a))
            if candidates:
                point=min(candidates,key=lambda p:abs(p[0]-xc)); rows.append((int(event),xc,point[1],point[2]))
    return pd.DataFrame(rows,columns=["event","column_center_x","mc_cross_y","mc_cross_z"])

def transverse_column_positions(runs,segments,track_points,use_selected=True):
    """Charge centroids in x columns and MC energy-centroid reference.

    Projection 2 (Z fibres) measures y and projection 0 (Y fibres) measures z.
    Their x positions are shifted to their common physical cube centre.
    """
    output={}
    for label,run in runs.items():
        hits=run["selected" if use_selected else "hits"].copy(); hits=hits[hits.projection.isin([0,2])]
        hits["column_center_x"]=np.where(hits.projection==0,hits.x+2.5,hits.x-2.5)
        y=_weighted_group(hits[hits.projection==2],"y").rename(columns={"multiplicity":"multiplicity_y","coc_y":"reco_y","weight_sum":"charge_y"})
        z=_weighted_group(hits[hits.projection==0],"z").rename(columns={"multiplicity":"multiplicity_z","coc_z":"reco_z","weight_sum":"charge_z"})
        output[label]=y.merge(z,on=["event","column_center_x"],how="inner")
    detector=segments[segments.detector==0].copy(); detector["mid_x"]=(detector.start_x+detector.stop_x)/2; detector["mid_y"]=(detector.start_y+detector.stop_y)/2; detector["mid_z"]=(detector.start_z+detector.stop_z)/2
    detector["column_center_x"]=np.floor(detector.mid_x/10.0)*10.0+5.0
    truth_rows=[]
    for (event,column),group in detector.groupby(["event","column_center_x"],sort=False):
        w=group.energy_deposit.to_numpy(float); total=w.sum()
        if total>0: truth_rows.append((int(event),float(column),float(np.dot(w,group.mid_y)/total),float(np.dot(w,group.mid_z)/total),total,len(group)))
    truth_coc=pd.DataFrame(truth_rows,columns=["event","column_center_x","mc_coc_y","mc_coc_z","mc_energy","mc_multiplicity"])
    centres={int(e):set(g.column_center_x) for e,g in truth_coc.groupby("event")}; crossings=_primary_crossings(track_points,centres)
    truth_coc=truth_coc.merge(crossings,on=["event","column_center_x"],how="inner")
    truth_coc["mc_coc_residual_y"]=truth_coc.mc_coc_y-truth_coc.mc_cross_y; truth_coc["mc_coc_residual_z"]=truth_coc.mc_coc_z-truth_coc.mc_cross_z
    for label in output:
        output[label]=output[label].merge(truth_coc,on=["event","column_center_x"],how="inner")
        output[label]["reco_residual_y"]=output[label].reco_y-output[label].mc_cross_y; output[label]["reco_residual_z"]=output[label].reco_z-output[label].mc_cross_z
    return output,truth_coc

def plot_transverse_multiplicity(columns,bins=np.arange(.5,15.6,1)):
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True); rows=[]
    for label,frame in columns.items():
        for ax,axis in zip(axes,("y","z")):
            values=frame[f"multiplicity_{axis}"].to_numpy(); colour="tab:orange" if label=="Track-aware fit" else "tab:blue" if label=="Previous fit" else None; ax.hist(values,bins=bins,histtype="step",density=True,lw=1.6,label=label,color=colour)
            rows.append({"method":label,"coordinate":axis,"columns":len(values),"mean_multiplicity":np.mean(values),"median_multiplicity":np.median(values)})
    for ax,axis in zip(axes,("Y from XY/Z-fibre view","Z from XZ/Y-fibre view")): ax.set(xlabel=f"fibre multiplicity per x column: {axis}",ylabel="fraction of columns"); ax.legend()
    return pd.DataFrame(rows),fig

def scatter_length_multiplicity_scan(available,lightmap_model,fit_filename="global_fit_columns_track_aware.root"):
    """Average populated-column multiplicity versus scattering length."""
    rows=[]
    for (model,scatter),directory in sorted(available.items(),key=lambda item:item[0][1]):
        if model!=lightmap_model: continue
        path=Path(directory)/fit_filename
        if not path.exists(): continue
        print(f"Multiplicity {scatter:g} mm: {path}")
        hits=tree_frame(path,"global_fit_fibres",["event","projection","x","y","z","charge"])
        for projection,coordinate,shift in ((2,"y",-2.5),(0,"z",2.5)):
            view=hits[hits.projection==projection].copy(); view["column_center_x"]=view.x+shift
            counts=view.groupby(["event","column_center_x"]).size().to_numpy(float)
            rows.append({"scatter_length_mm":scatter,"coordinate":coordinate,"columns":len(counts),"mean_multiplicity":counts.mean(),"sem_multiplicity":counts.std(ddof=1)/np.sqrt(len(counts)),"median_multiplicity":np.median(counts)})
    data=pd.DataFrame(rows)
    if data.empty: raise FileNotFoundError(f"No {fit_filename} files for {lightmap_model}")
    combined=data.groupby("scatter_length_mm",as_index=False).agg(columns=("columns","sum"),mean_multiplicity=("mean_multiplicity","mean"),sem_multiplicity=("sem_multiplicity",lambda x:np.sqrt(np.sum(np.square(x)))/len(x)),median_multiplicity=("median_multiplicity","mean"))
    fig,ax=plt.subplots(figsize=(8.5,5),constrained_layout=True)
    combined=combined.sort_values("scatter_length_mm")
    ax.errorbar(combined.scatter_length_mm,combined.mean_multiplicity,yerr=combined.sem_multiplicity,marker="o",color="black",capsize=3)
    ax.set(xlabel="scattering length [mm]",ylabel="mean selected-fibre multiplicity per populated x column",title=f"{lightmap_model}: average transverse multiplicity"); ax.grid(alpha=.25)
    return combined,fig

def plot_transverse_column_resolution(columns,truth_coc,display_range=(-10,10),fit_range=(-3,3),bins=100):
    from column_fit_notebook_tools import gaussian_peak_fit
    fig,axes=plt.subplots(1,2,figsize=(13,5),constrained_layout=True); rows=[]
    samples={"MC energy center of charge":truth_coc}
    samples.update(columns)
    for ax,coordinate in zip(axes,("y","z")):
        for label,frame in samples.items():
            field=f"mc_coc_residual_{coordinate}" if label.startswith("MC ") else f"reco_residual_{coordinate}"; values=frame[field].replace([np.inf,-np.inf],np.nan).dropna().to_numpy(float)
            colour="black" if label.startswith("MC ") else "tab:orange" if label=="Track-aware fit" else "tab:blue" if label=="Previous fit" else ax._get_lines.get_next_color(); mean,sigma,entries,lo,hi,amplitude=gaussian_peak_fit(values,histogram_range=display_range,bins=bins,fixed_fit_range=fit_range)
            ax.hist(values,bins=bins,range=display_range,histtype="step",density=True,lw=1.6,color=colour,label=f"{label}: sigma={sigma:.3g} mm")
            if np.isfinite(sigma) and sigma>0:
                x=np.linspace(lo,hi,300); ax.plot(x,amplitude*np.exp(-.5*((x-mean)/sigma)**2),"--",color=colour,lw=1.2)
            rows.append({"sample":label,"coordinate":coordinate,"columns":len(values),"gaussian_mean_mm":mean,"gaussian_sigma_mm":sigma,"fit_low_mm":lo,"fit_high_mm":hi,"fit_entries":entries})
        ax.set(xlabel=f"{coordinate.upper()} center-of-charge minus MC crossing [mm]",ylabel="density",title=f"Per-column {coordinate.upper()} resolution"); ax.legend(fontsize=8)
    return pd.DataFrame(rows),fig

def common_events(runs): return sorted(set.intersection(*(set(r["hits"].event.astype(int)) for r in runs.values())))
def _edges(x):
    lo,hi=float(np.min(x))-5,float(np.max(x))+5
    return np.linspace(lo,hi,max(1,int(round((hi-lo)/10)))+1)

def draw_paired_event(runs,truth,event,threshold=None,log_charge=True,show_mc=True,show_fit=True,show_selected=True,shared_charge_scale=False,charge_min=1.0,charge_max=150.0):
    """Legacy left, double-HG right; rows are XY, XZ and YZ."""
    event=int(event); preferred=[x for x in ("Legacy","Double-HG","Previous fit","Track-aware fit") if x in runs]; labels=preferred or list(runs); data={}; allq=[]
    for label in labels:
        fit=runs[label]["fit"]; fit=fit[fit.event==event]
        stored=float(fit.minimum_charge.iloc[0]) if not fit.empty else 0.; cut=stored if threshold is None else float(threshold)
        hits=runs[label]["hits"]; hits=hits[(hits.event==event)&(hits.charge>=cut)]
        selected=runs[label]["selected"]; selected=selected[selected.event==event]
        data[label]=(hits,selected,fit,cut); allq.extend(hits.loc[hits.charge>0,"charge"])
    mc=truth[(truth.event==event)&(truth.detector==0)&(truth.primary_id==1)]
    shared=LogNorm(charge_min,charge_max) if log_charge else Normalize(charge_min,charge_max)
    fig,axes=plt.subplots(3,2,figsize=(14,14),constrained_layout=True,squeeze=False)
    for col,label in enumerate(labels):
        hits,selected,fit,cut=data[label]; q=hits.loc[hits.charge>0,"charge"]
        norm=shared
        image=None
        for row,(projection,(view_name,a,b)) in enumerate(VIEWS.items()):
            ax=axes[row,col]; view=hits[hits.projection==projection]
            if not view.empty:
                ae,be=_edges(view[a]),_edges(view[b]); charge,_,_=np.histogram2d(view[a],view[b],bins=(ae,be),weights=view.charge)
                image=ax.pcolormesh(ae,be,charge.T,cmap="viridis",norm=norm,shading="flat"); ax.set_xlim(ae[0],ae[-1]); ax.set_ylim(be[0],be[-1])
            if show_selected:
                accepted=selected[selected.projection==projection]; ax.scatter(accepted[a],accepted[b],facecolors="none",edgecolors="black",s=34,linewidths=.8,zorder=4)
            if show_mc:
                for s in mc.itertuples(index=False): ax.plot([getattr(s,"start_"+a),getattr(s,"stop_"+a)],[getattr(s,"start_"+b),getattr(s,"stop_"+b)],color="#159447",lw=1.6,zorder=5)
            if show_fit and not fit.empty:
                f=fit.iloc[0]; p=np.array([f.fit_x,f.fit_y,f.fit_z]); d=np.array([f.fit_dx,f.fit_dy,f.fit_dz]); ai,bi="xyz".index(a),"xyz".index(b); t=np.array([-2500.,2500.])
                ax.plot(p[ai]+t*d[ai],p[bi]+t*d[bi],color="red",lw=2,zorder=6)
            status=int(fit.status.iloc[0]) if not fit.empty else -1
            # ROOT's TH2 display fills the pad.  Keeping that convention also
            # magnifies the narrow transverse extent, which is the useful
            # feature of this track-width display.
            ax.set(title=f"{label} — {view_name} ({len(view)} fibres, status {status})",xlabel=f"{a.upper()} [mm]",ylabel=f"{b.upper()} [mm]"); ax.set_aspect("auto")
        if image is not None: fig.colorbar(image,ax=axes[:,col],shrink=.72,pad=.02,label="fibre charge [PE]")
    handles=[Line2D([0],[0],color="#159447",lw=2,label="MC truth"),Line2D([0],[0],color="red",lw=2,label="GlobalFit"),Line2D([0],[0],marker="o",markerfacecolor="none",markeredgecolor="black",linestyle="none",label="fit-selected fibre")]
    scatter=runs[labels[0]]["scatter"]; scatter_text=f" — scattering length {scatter:g} mm" if np.isfinite(scatter) else ""
    fig.legend(handles=handles,loc="upper center",ncol=3,bbox_to_anchor=(.5,1.015)); fig.suptitle(f"Event {event}{scatter_text} — colour range {charge_min:g}–{charge_max:g} PE",y=1.035,fontsize=15)
    return fig

def event_control_table(runs):
    out=[]
    for label,r in runs.items():
        x=r["hits"].groupby("event").charge.agg(hit_count="size",total_charge="sum",maximum_charge="max").reset_index(); x["model"]=label; out.append(x)
    return pd.concat(out,ignore_index=True)

def plot_input_controls(runs):
    data=event_control_table(runs); fig,axes=plt.subplots(1,3,figsize=(16,4.5),constrained_layout=True)
    for label,x in data.groupby("model",sort=False):
        for ax,key,name in zip(axes,("hit_count","total_charge","maximum_charge"),("fibre hits/event","total charge/event [PE]","maximum charge/event [PE]")): ax.hist(x[key],bins=40,histtype="step",density=True,lw=1.6,label=label); ax.set(xlabel=name,ylabel="density")
    axes[0].legend(); return data,fig

def plot_fit_controls(runs):
    out=[]
    for label,r in runs.items():
        x=r["fit"].copy(); x["model"]=label; a=x[["seed_dx","seed_dy","seed_dz"]].to_numpy(float); b=x[["fit_dx","fit_dy","fit_dz"]].to_numpy(float); a/=np.linalg.norm(a,axis=1)[:,None]; b/=np.linalg.norm(b,axis=1)[:,None]; x["seed_fit_angle_deg"]=np.degrees(np.arccos(np.clip(np.abs(np.sum(a*b,axis=1)),0,1))); x["chi2_ndof"]=np.where(x.ndof>0,x.chi2/x.ndof,np.nan); out.append(x)
    data=pd.concat(out,ignore_index=True); fig,axes=plt.subplots(2,2,figsize=(14,9),constrained_layout=True); pd.crosstab(data.model,data.status,normalize="index").plot.bar(stacked=True,ax=axes[0,0]); axes[0,0].set(ylabel="event fraction",title="GlobalFit status")
    for label,x in data.groupby("model",sort=False): axes[0,1].hist(x.edm[np.isfinite(x.edm)],bins=np.logspace(-6,1,45),histtype="step",density=True,label=label); axes[1,0].hist(x.seed_fit_angle_deg,bins=np.linspace(0,15,61),histtype="step",density=True,label=label); axes[1,1].hist(x.chi2_ndof.clip(upper=100),bins=50,histtype="step",density=True,label=label)
    axes[0,1].set(xscale="log",xlabel="EDM",ylabel="density"); axes[1,0].set(xlabel="seed-fit angle [deg]",ylabel="density"); axes[1,1].set(xlabel="chi2/ndof (clipped at 100)",ylabel="density"); axes[0,1].legend(); return data,fig

def interactive_pair(runs,truth):
    import ipywidgets as widgets
    from IPython.display import display
    events=common_events(runs); event=widgets.SelectionSlider(options=events,value=events[0],description="Event",continuous_update=False,layout=widgets.Layout(width="720px")); threshold=widgets.FloatText(value=-1,description="Threshold"); log=widgets.Checkbox(value=True,description="Log charge"); shared=widgets.Checkbox(value=False,description="Shared scale"); mc=widgets.Checkbox(value=True,description="MC truth"); fit=widgets.Checkbox(value=True,description="Fit"); selected=widgets.Checkbox(value=True,description="Selected"); output=widgets.Output()
    def redraw(*_):
        with output: output.clear_output(wait=True); draw_paired_event(runs,truth,event.value,None if threshold.value<0 else threshold.value,log.value,mc.value,fit.value,selected.value,shared.value); plt.show()
    for c in (event,threshold,log,shared,mc,fit,selected): c.observe(redraw,names="value")
    display(widgets.VBox([event,widgets.HBox([threshold,log,shared,mc,fit,selected])]),output); redraw(); return output


def draw_rl_event(run, truth, rl_columns, event, solution="regularized",
                  threshold=None, log_charge=True, charge_min=1., charge_max=150.):
    """Fibre-charge display with MC truth and only the movable RL trajectory."""
    event = int(event)
    hits = run["hits"]; hits = hits[hits.event == event].copy()
    if threshold is not None:
        hits = hits[hits.charge >= float(threshold)]
    columns = rl_columns[rl_columns.event == event].sort_values("column_index")
    if columns.empty:
        raise ValueError(f"Event {event} is absent from the RL output")
    prefix = "best_nll_" if solution == "best-nll" else ""
    coordinates = {
        axis: columns[prefix + f"centre_{axis}"].to_numpy(float)
        for axis in "xyz"
    }
    mc = truth[(truth.event == event) & (truth.detector == 0) &
               (truth.primary_id == 1)]
    norm = (LogNorm(charge_min, charge_max) if log_charge else
            Normalize(charge_min, charge_max))
    colour = "#d336c2" if solution == "regularized" else "#f28e2b"
    fig, axes = plt.subplots(3, 1, figsize=(8.5, 14), constrained_layout=True)
    image = None
    for axis, (projection, (view_name, a, b)) in zip(axes, VIEWS.items()):
        view = hits[hits.projection == projection]
        if not view.empty:
            ae, be = _edges(view[a]), _edges(view[b])
            charge, _, _ = np.histogram2d(view[a], view[b], bins=(ae, be),
                                           weights=view.charge)
            image = axis.pcolormesh(ae, be, charge.T, cmap="viridis", norm=norm,
                                    shading="flat")
            axis.set_xlim(ae[0], ae[-1]); axis.set_ylim(be[0], be[-1])
        for segment in mc.itertuples(index=False):
            axis.plot([getattr(segment, "start_" + a), getattr(segment, "stop_" + a)],
                      [getattr(segment, "start_" + b), getattr(segment, "stop_" + b)],
                      color="#159447", lw=1.8, zorder=5)
        axis.plot(coordinates[a], coordinates[b], "o-", ms=3, lw=1.8,
                  color=colour, zorder=6)
        axis.set(title=f"{view_name} ({len(view)} fibres)", xlabel=f"{a.upper()} [mm]",
                 ylabel=f"{b.upper()} [mm]")
    if image is not None:
        fig.colorbar(image, ax=axes, shrink=.72, pad=.02, label="fibre charge [PE]")
    handles = [Line2D([0], [0], color="#159447", lw=2, label="MC truth"),
               Line2D([0], [0], color=colour, marker="o", lw=2,
                      label=f"RL {solution}")]
    fig.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(.5, 1.01))
    fig.suptitle(f"Event {event}: RL-refined trajectory — colour range "
                 f"{charge_min:g}–{charge_max:g} PE", y=1.025, fontsize=15)
    return fig


def interactive_rl_event(run, truth, rl_columns, initial_event=None):
    """Interactive event selector for the fibre/MC/RL-only display."""
    import ipywidgets as widgets
    from IPython.display import display
    events = sorted(set(run["hits"].event.astype(int)) &
                    set(rl_columns.event.astype(int)))
    if not events:
        raise ValueError("No events are common to the fibre hits and RL output")
    initial = events[0] if initial_event not in events else int(initial_event)
    event = widgets.SelectionSlider(options=events, value=initial, description="Event",
                                    continuous_update=False,
                                    layout=widgets.Layout(width="720px"))
    solution = widgets.Dropdown(options=["regularized", "best-nll"],
                                value="regularized", description="RL solution")
    threshold = widgets.FloatText(value=-1, description="Threshold")
    log = widgets.Checkbox(value=True, description="Log charge")
    output = widgets.Output()
    def redraw(*_):
        with output:
            output.clear_output(wait=True)
            draw_rl_event(run, truth, rl_columns, event.value, solution.value,
                          None if threshold.value < 0 else threshold.value,
                          log.value, 1., 150.)
            plt.show()
    for control in (event, solution, threshold, log):
        control.observe(redraw, names="value")
    display(widgets.VBox([event, widgets.HBox([solution, threshold, log])]), output)
    redraw()
    return output
