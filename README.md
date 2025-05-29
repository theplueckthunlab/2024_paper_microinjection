# Microinjection Raw Data

This repository contains raw and processed data from a series of experiments involving microinjection and live-cell fluorescence microscopy. The data underlies multiple scientific publications, including both published and under-review manuscripts.

## Related Publication

A major portion of this dataset has been used in the following peer-reviewed study:

**Vukovic, D., Winkelvoß, D., Kapp, J.N., Hänny, A.C., Bürgisser, H., Riermeier, L., Udovcic, A., Tiefenboeck, P. and Plückthun, A.**  
*Protein degradation kinetics measured by microinjection and live-cell fluorescence microscopy*.  
**Scientific Reports**, 14(1), p.27153 (2024).  
[https://www.nature.com/articles/s41598-024-76224-0](https://www.nature.com/articles/s41598-024-76224-0)

## File Structure

### Raw Data

- **`microinjection_raw_data.zip`**  
  ZIP archive containing `.csv` files. Each file corresponds to one overnight experiment, with data extracted after image segmentation.

### Metadata and Analyte Harmonization

- **`conversion_analytes_v17.xlsx`**  
  Master lookup table for harmonizing analyte (sample) IDs across batches. Used for both data cleaning and combining equivalent experimental conditions.

- **`analytetranslation_paper02.csv`**  
  Custom analyte naming used specifically in the second manuscript. Builds upon `conversion_analytes_v17.xlsx`, adapting analyte names for clarity and consistency in publication.

### Data Processing

- **`02_single_cell_fit_data_v18.R`**  
  R script for model fitting of single-cell fluorescence data, with outlier exclusion. Generates the cleaned dataset used across figures.

### Data Visualization

- **`03_single_cell_fit_raw_graph_v5.R`**  
  Visualization script operating on the processed data from `02_single_cell_fit_data_v18.R`.  
  Includes:
  - Per-cell raw and fitted curves
  - Bar graphs of decay rates or intensities
  - Summarized views of all cells injected with a given analyte

## Data Description

Each `.csv` file in `microinjection_raw_data.zip` contains background-subtracted fluorescence data per cell. Column details:

| Column Name       | Description |
|------------------|-------------|
| `exp.id`         | Unique experiment ID (e.g., `2019-11-01`, or `2019-11-01-1` for multiple on same day). |
| `pos.id`         | Microscope field ID; multiple cells can be recorded in one view. |
| `cell.id`        | Unique per-cell ID within experiment (no reuse). |
| `tp.id`          | Time point index. |
| `tp.posix`       | POSIX timestamp of acquisition. |
| `tp.postinj`     | Time post-injection in hours (±1 min). |
| `epi.405.int`    | 405 nm fluorescence (background-subtracted). |
| `epi.488.int`    | 488 nm fluorescence (background-subtracted). |
| `epi.561.int`    | 561 nm fluorescence (background-subtracted). |
| `epi.640.int`    | 640 nm fluorescence (background-subtracted). |
| `analyte`        | Internal freezer sample ID. See conversion files for harmonized or publication-ready names. |

## Notes

- Image segmentation and data extraction were performed with **CellProfiler** to ensure consistency across fluorescence channels and time points.
- The full processing pipeline ensures clean, reproducible analysis—from raw intensities to statistical modeling and figure creation.
- For the second manuscript, `analytetranslation_paper02.csv` defines final analyte names as used in the figures and text.

---

If you use this repository or any of its contents, please cite the above publication.  
For questions, feedback, or collaborations, open an issue or contact the maintainer.
