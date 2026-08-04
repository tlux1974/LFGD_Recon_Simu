cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
find_package(nd280SoftwarePolicy 3.1 REQUIRED)
project(LFGDReconSimu)

ND280_PROJECT(LFGDReconSimu 0.1.0)
ND280_USE(ROOT *)
ND280_USE(oaEvent *)
ND280_USE(oaGeomInfo *)
ND280_USE(oaRuntimeParams *)
ND280_USE(hfgRecon *)
# This package contains only an executable.  This must be declared before
# ND280_EXECUTABLE so the ND280 policy does not link a nonexistent package
# library (-lLFGDReconSimu).
ND280_NO_LIBRARY(LFGDReconSimu)
ND280_EXECUTABLE(LFGDFLATTREE ../src/LFGDFlatTree.cxx)
ND280_END_PROJECT()
