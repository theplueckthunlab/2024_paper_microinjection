# Microinjection Raw Data

This repository contains raw and processed data from a series of experiments involving microinjection and live-cell fluorescence microscopy. The data underlies multiple scientific publications, including both published and under-review manuscripts.

## Related Publication

A major portion of this dataset has been used in the following peer-reviewed study:

**Vukovic, D., Winkelvoß, D., Kapp, J.N., Hänny, A.C., Bürgisser, H., Riermeier, L., Udovcic, A., Tiefenboeck, P. and Plückthun, A.**  
*Protein degradation kinetics measured by microinjection and live-cell fluorescence microscopy*.  
**Scientific Reports**, 14(1), p.27153 (2024).  
[https://www.nature.com/articles/s41598-024-76224-0](https://www.nature.com/articles/s41598-024-76224-0)

---

## File Structure

### Raw Data

- **`microinjection_raw_data.zip`**  
  Contains `.csv` files with extracted single-cell fluorescence intensity data from each experiment.

- **`conversion_analytes_v17.xlsx`**  
  Analyte harmonization table to match freezer samples across experiments.

- **`analytetranslation_paper02.csv`**  
  Adjusted analyte naming used for the second manuscript. Builds upon the master analyte conversion table.

### Image-to-Data Pipeline

- **`segmentation_pipeline.zip`**  
  A full example pipeline from raw microscope imagery to a final `.csv` result (`2023-04-13_output1.csv`) using multiple tools and languages.

  This includes:
  
  **Batch and Script Files**  
  | Filename | Purpose |
  |----------|---------|
  | `01_ff_folder_creation.bat` | Create folder structure for processing |
  | `02_sort_files_to_folders.bat` | Sort microscope images into position-based folders |
  | `03_creation_times.m` | Extract creation timestamps from `.tiff` metadata |
  | `04_copy_timepoints_CSV.bat` | Copy and prepare timestamp files for later steps |
  | `05_flatfielding.m` | Perform flatfield correction using dark/flat-field median images |
  | `06_rename_XX_numb.bat` | Standardize file naming for CellProfiler compatibility |
  | `07_CSV_for_cellprofiler_multich.R` | Generate file name tables for headless CellProfiler runs |
  | `08_segmentation_bgsub_min_bg80_triplech_v5.cppipe` | CellProfiler pipeline file for segmentation and background subtraction |
  | `09_cellprofiler_batch_v2_multich_tr1.bat` | Launch CellProfiler in parallel for batch processing |
  | `10_sort_CP_output_tr1.bat` | Relocate CellProfiler output files |
  | `11_copy_roifilev.bat` | Combine segmented ROIs for multiple cells in the same field |
  | `12_extract_single_cell_data_v7_4channels.R` | Extract and compile intensity data from segmented output |
  | `2023-04-13_output1.csv` | Final single-cell result for this example experiment |

### Data Processing

- **`02_single_cell_fit_data_v18.R`**  
  Performs model fitting of fluorescence decay curves and outlier filtering. Produces cleaned dataset for all figures.

### Data Visualization

- **`03_single_cell_fit_raw_graph_v5.R`**  
  Generates exploratory and publication-ready plots:
  - Raw and fitted data per cell
  - Analyte-level bar plots
  - Summary views for grouped analyte conditions

---

## Software Used

The pipeline is cross-platform and relies on the following software components:

- **R** (for data structuring, fitting, and plotting)
- **Batch scripting** (`.bat` files) for file handling on **Windows**
- **MATLAB R2021** for TIFF metadata extraction and flatfielding
- **CellProfiler 2.2.0 (rev ac0529e)** for cell segmentation and intensity extraction

---

## Data Description

Each `.csv` file in `microinjection_raw_data.zip` contains fluorescence intensities per cell per time point:

| Column Name       | Description |
|------------------|-------------|
| `exp.id`         | Unique experiment ID (e.g., `2019-11-01`, `2019-11-01-1`) |
| `pos.id`         | Microscope position ID |
| `cell.id`        | Unique cell ID per experiment |
| `tp.id`          | Time point index |
| `tp.posix`       | POSIX timestamp |
| `tp.postinj`     | Time after injection (±1 min accuracy) |
| `epi.405.int` to `epi.640.int` | Background-subtracted fluorescence intensities |
| `analyte`        | Freezer sample ID (see conversion files for harmonized names) |

---

## Notes

- Segmentation and analysis are performed in a fully documented, reproducible way.
- Each component in the pipeline can be reused or adapted to new datasets.
- Naming conventions were adapted between papers for clarity—see the analyte translation files for mappings.

---

If you use this repository or any part of the pipeline, **please cite** the publication above.  
For questions, suggestions, or collaborations, feel free to open an issue or contact the repository maintainer.
