# Mu2e tracker alignment track analyzer

This folder contains a focused workflow for studying wire alignment residuals with selected reconstructed tracks. It is a modification of the [Mu2e GitHub Tutorial](https://github.com/Mu2e/Tutorial), particularly its examples of building and running an `art` analyzer module. The tutorial’s [Module Writing Tutorial](https://mu2ewiki.fnal.gov/wiki/Module_Writing_Tutorial) and [Running Art Tutorial](https://mu2ewiki.fnal.gov/wiki/Running_Art_Tutorial) explain the framework, FHiCL configuration, and job execution.

## Contents

- `src/AlignmentTrackAnalyzer_module.cc` — the only C++ module in this package. It selects KalSeed tracks and produces track, hit, and wire-alignment diagnostic histograms and output suitable for residual studies.
- `fcl/` — the five retained configurations: `singlepanel.fcl`, `twopanel.fcl`, `fourpanelV.fcl`, `sixpanelZ.fcl`, and `10panelZ.fcl`. Each configures the analyzer’s track and panel-intersection selection.
- `drawPanelDifferences2.cpp` — reads one ROOT file produced by the analyzer, performs the residual fits, and produces diagnostic PDFs and CSV files.
- `plot_DOCA_XYV_folder_input.py` — reads the CSV outputs from `drawPanelDifferences2.cpp` and creates PDFs.
- `grid.sh` — submits reconstruction using the supplied grid workflow and calibration inputs.
- `gridMacro.sh` — submits the analyzer workflow using an embedded FHiCL configuration.
- `TrackerCalib_Palo/` — calibration, timing, panel map, and alignment text/FHiCL files used by reconstruction.
- `filelist.txt` — example input file list.

## Workflow overview

1. `grid.sh` runs reconstruction and produces the reconstructed art output. The reconstruction configuration refers to the calibration files in `TrackerCalib_Palo/`.
2. `gridMacro.sh` runs `AlignmentTrackAnalyzer_module.cc` on reconstructed events. The module performs track selection and histogramming and writes a ROOT file.
3. Run `drawPanelDifferences2.cpp` on that ROOT file to fit the wire-alignment residual distributions and save CSV/PDF results.
4. Run `plot_DOCA_XYV_folder_input.py` on the resulting CSV files to create the summary PDFs.

The scripts contain site-specific grid submission settings and the reconstruction input dataset name. Review those values for the target Mu2e environment and dataset before submitting. The selected FHiCL file can be changed in `gridMacro.sh` by editing the `--embed` argument.

## Build and run

Build this package within a compatible Mu2e Offline `art` environment, following the build and setup procedure in the Mu2e tutorials. Run the analyzer with one of the FHiCL files in `fcl/`, providing reconstructed art input and the appropriate tracker geometry/calibration conditions. Exact software releases, dataset names, and grid resource settings are environment-dependent.
