# Mu2e tracker alignment with reconstructed tracks

This repository contains a focused workflow for studying wire and panel alignment residuals with selected reconstructed tracks. The analyzer selects tracks and fills diagnostic histograms; ROOT-based post-processing fits the residuals and produces plots and tables used to evaluate alignment corrections.

This is the track-based stage of a broader geometry workflow. [MetrologyMain](https://github.com/dylancpalo/MetrologyMain) combines the camera, X-ray, tracker-frame, and mechanical measurements and computes offline values for alignment tables. Track-based studies provide an independent detector-response check and help validate or refine the geometry using reconstructed tracks.

## Contents

- `src/AlignmentTrackAnalyzer_module.cc` — the only C++ `art` analyzer module in this package; selects tracks and produces track, hit, and wire-alignment diagnostic histograms.
- `fcl/` — the retained FHiCL configurations for single-, two-, four-, six-, and ten-panel track selections.
- `drawPanelDifferences2.cpp` — reads an analyzer ROOT file, fits residual distributions, and produces diagnostic PDFs and CSV files.
- `plot_DOCA_XYV_folder_input.py` — reads CSV outputs from `drawPanelDifferences2.cpp` and creates summary PDFs.
- `grid.sh` — example grid workflow for reconstruction.
- `gridMacro.sh` — example grid workflow for running the analyzer.
- `TrackerCalib_Palo/` — calibration, timing, panel-map, and alignment text/FHiCL inputs used by the supplied workflows.
## Workflow

1. Run reconstruction with the desired calibration and alignment inputs.
2. Run `AlignmentTrackAnalyzer_module.cc` on the reconstructed events with a suitable FHiCL selection. The module performs track selection and histogramming and writes a ROOT file.
3. Run `drawPanelDifferences2.cpp` on the analyzer output to fit wire-alignment residuals and create CSV/PDF diagnostics.
4. Use `plot_DOCA_XYV_folder_input.py` to turn the CSV results into summary plots.
5. Compare the residual patterns with the metrology-derived geometry and use the combined evidence to assess alignment updates and weakly constrained modes.

The package follows the Mu2e `art` analyzer workflow described in the [Mu2e GitHub Tutorial](https://github.com/Mu2e/Tutorial), including its [Module Writing Tutorial](https://mu2ewiki.fnal.gov/wiki/Module_Writing_Tutorial) and [Running Art Tutorial](https://mu2ewiki.fnal.gov/wiki/Running_Art_Tutorial).

## Build and run

Build and run this package inside a compatible Mu2e Offline `art` environment. The scripts and FHiCL files contain environment-specific release, dataset, and grid settings; review those values and the input calibration files for the target dataset before submitting jobs.

## Related repositories

- [CameraMetrologyImageGrab_Processing](https://github.com/dylancpalo/CameraMetrologyImageGrab_Processing) captures and processes camera images and imports station-assembly measurements into PostgreSQL.
- [MetrologyMain](https://github.com/dylancpalo/MetrologyMain) combines camera, X-ray, tracker-frame, and survey data and produces offline fit values for alignment tables.
- [UltemAnalysis](https://github.com/dylancpalo/UltemAnalysis) supports Ultem quality control and machining decisions.
