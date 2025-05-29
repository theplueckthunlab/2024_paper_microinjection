# Microinjection Raw Data

This repository contains raw and processed data from a series of experiments involving microinjection and live-cell fluorescence microscopy. The data underlies multiple scientific publications, including both published and in-preparation manuscripts.

---

## Related Publications

### 1. Protein Degradation Kinetics by Microinjection

**Vukovic, D.**, Winkelvoß, D., Kapp, J.N., Hänny, A.C., Bürgisser, H., Riermeier, L., Udovcic, A., Tiefenboeck, P., and Plückthun, A.  
**Protein degradation kinetics measured by microinjection and live-cell fluorescence microscopy**  
*Scientific Reports*, 14(1), p.27153 (2024).  
[https://www.nature.com/articles/s41598-024-76224-0](https://www.nature.com/articles/s41598-024-76224-0)

---

### 2. Molecular Features Defining the Efficiency of bioPROTACs

Winkelvoß, D., **Vukovic, D.**, Hänny, A.-C., Riermeier, L., Udovcic, A., Honegger, A., Mittl, P., Michel, E., Hansen, S., and Kolibius, J.  
**Molecular features defining the efficiency of bioPROTACs**  
*Communications Biology* (in review).  
**Manuscript #: COMMSBIO-24-7628B**

---

### 3. Quantitative Comparison of bioPROTAC Degradation Strategies *(Manuscript in Preparation)*

**Vukovic, D.**, Winkelvoß, D., Udovcic, A., Riermeier, L., Jaime-Figueroa, S., Crews, C., and Plückthun, A.  
**Quantitative degradation rate assessment of bioPROTACs comparing peptide degrons, E3 domains, adapters and conjugated small molecules**  
*Manuscript in preparation.*

---

## File Structure

### Raw Data

- **`microinjection_raw_data.zip`**  
  Contains `.csv` files with extracted single-cell fluorescence intensity data from each experiment.

### Metadata and Analyte Harmonization

- **`conversion_analytes_v17.xlsx`**  
  Analyte harmonization table to match freezer samples across experiments.

- **`analytetranslation_paper02.csv`**  
  Adjusted analyte naming used specifically in paper #2 (bioPROTAC efficiency study), derived from the base conversion.

### Image-to-Data Pipeline

- **`segmentation_pipeline.zip`**  
  Complete example workflow (based on `2023-04-13_output1.csv`) that transforms raw microscopy images into the final `.csv` used in downstream analysis.

### Data Processing

- **`02_single_cell_fit_data_v18.R`**  
  Fits decay models to single-cell fluorescence trajectories, excluding outliers and producing clean data for figures.

### Data Visualization

- **`03_single_cell_fit_raw_graph_v5.R`**  
  Visualizes both raw and fitted data:
  - Per-cell fluorescence traces
  - Group-level bar graphs
  - Summary metrics by analyte and condition

---

## Segmentation Pipeline Contents

Contained in `segmentation_pipeline.zip`. Scripts automate image handling, flatfielding, segmentation, and intensity extraction:

| Script/File | Purpose |
|-------------|---------|
| `01_ff_folder_creation.bat` | Create initial folder structure |
| `02_sort_files_to_folders.bat` | Organize image files |
| `03_creation_times.m` | Extract image timestamps |
| `04_copy_timepoints_CSV.bat` | Copy and reformat time metadata |
| `05_flatfielding.m` | Apply flatfield corrections |
| `06_rename_XX_numb.bat` | Rename files for CellProfiler compatibility |
| `07_CSV_for_cellprofiler_multich.R` | Create input tables for headless runs |
| `08_segmentation_bgsub_min_bg80_triplech_v5.cppipe` | CellProfiler pipeline settings |
| `09_cellprofiler_batch_v2_multich_tr1.bat` | Run segmentation in parallel |
| `10_sort_CP_output_tr1.bat` | Reorganize segmentation outputs |
| `11_copy_roifilev.bat` | Consolidate POI/ROI data per cell |
| `12_extract_single_cell_data_v7_4channels.R` | Final extraction of time-series intensity data |
| `2023-04-13_output1.csv` | Final result from this pipeline |

---

## Software Used

- **R** – for data wrangling, fitting, and plotting
- **Windows Batch (.bat)** – for file management and automation
- **MATLAB R2021** – for TIFF metadata handling and flatfielding
- **CellProfiler 2.2.0 (rev ac0529e)** – for segmentation and feature quantification

---

## Data Description

| Column Name       | Description |
|------------------|-------------|
| `exp.id`         | Unique experiment identifier (e.g., `2019-11-01`, `2019-11-01-1`) |
| `pos.id`         | Microscope position |
| `cell.id`        | Unique per-cell identifier |
| `tp.id`          | Time point index |
| `tp.posix`       | Acquisition timestamp in POSIX format |
| `tp.postinj`     | Time since injection (in hours) |
| `epi.405.int` to `epi.640.int` | Fluorescence intensities for 4 channels (back
