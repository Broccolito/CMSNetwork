# CMSNetwork - Warp Initialization

## Project Overview
This project analyzes convergent molecular signatures (CMS) networks between Tibetan and Andean high-altitude populations, focusing on genetic adaptations and network overlaps.

## Directory Structure
```
CMSNetwork/
├── warp.md                              # This initialization file
├── .Rproj.user/                        # R project settings
├── annotated_network_overlap/           # Network overlap analysis results
├── csv_files/                           # CSV data files
├── manuscript/                          # Manuscript and publication files
├── network_overlap/                     # Network overlap computations
├── parquet_files/                       # Parquet format data files
├── run_v2g0.1/, run_v2g0.2/           # V2G analysis runs
├── andean_cms_annotated.csv            # Annotated Andean CMS data
├── andean_cms_components.csv           # Andean CMS network components
├── tibetan_cms_annotated.csv           # Annotated Tibetan CMS data
├── tibetan_cms_components.csv          # Tibetan CMS network components
├── v2g.csv.gz                          # Variant-to-gene mapping data
├── all_tibetan_andean_overlap.xlsx     # Combined overlap analysis
├── tibetan_andean_overlap_significant.xlsx # Significant overlaps
└── Various R scripts and analysis files
```

## Key Files & Their Purpose

### Data Files
- **v2g.csv.gz**: Large variant-to-gene mapping dataset (3.7GB)
- **tibetan_cms_components.csv**: Tibetan CMS network components (675MB)
- **andean_cms_components.csv**: Andean CMS network components (698MB)
- **tibetan_cms_annotated.csv**: Annotated Tibetan CMS data (5.2MB)
- **andean_cms_annotated.csv**: Annotated Andean CMS data (40.4MB)

### Analysis Scripts
- **determine_network_overlap.R**: Main script for network overlap analysis
- **determine_raw_gene_convergence.R**: Gene-level convergence analysis
- **determine_overlap.R**: General overlap determination
- **cleanup.R**: Data cleaning utilities
- **run_tkoi.R**: TKOI analysis runner
- **process_v2g.ipynb**: Jupyter notebook for V2G data processing

### Results & Outputs
- **network_convergence.png**: Network convergence visualization
- **gene_level_convergence.png**: Gene-level convergence plot
- **all_tibetan_andean_overlap.xlsx**: Complete overlap results
- **tibetan_andean_overlap_significant.xlsx**: Filtered significant results

## Quick Start Commands

### Environment Setup
```bash
# Load R environment
R

# Install required packages (if needed)
# install.packages(c("tidyverse", "igraph", "networkD3", "openxlsx"))
```

### Common Analysis Tasks
```bash
# Run network overlap analysis
Rscript determine_network_overlap.R

# Run gene convergence analysis
Rscript determine_raw_gene_convergence.R

# Process V2G data
jupyter notebook process_v2g.ipynb
```

### Data Exploration
```bash
# Quick look at CMS components
head -n 5 tibetan_cms_components.csv
head -n 5 andean_cms_components.csv

# Check file sizes
ls -lh *cms*.csv

# Explore parquet files
ls -la parquet_files/ | head -10
```

## Analysis Workflow

1. **Data Preparation**: Process V2G mappings and CMS components
2. **Network Construction**: Build networks from CMS data
3. **Overlap Analysis**: Determine overlaps between Tibetan and Andean networks
4. **Annotation**: Add functional annotations to significant overlaps
5. **Visualization**: Generate convergence plots and network visualizations
6. **Results Export**: Save results to Excel files for further analysis

## Key Research Questions
- What are the convergent molecular signatures between Tibetan and Andean populations?
- Which genes show evidence of convergent adaptation to high altitude?
- How do the network structures compare between the two populations?
- What functional pathways are enriched in the overlapping networks?

## Development Notes
- Project uses R for primary analysis with some Python (Jupyter) components
- Large datasets are stored in compressed formats (CSV.gz, parquet)
- Results are exported to Excel for collaboration and publication
- Network analysis uses igraph and custom algorithms
- V2G analysis incorporates multiple evidence sources

## Getting Help
- Check R script comments for detailed parameter explanations
- Review Jupyter notebook for V2G processing methodology
- Examine CSV file headers to understand data structure
- Use `str()` and `summary()` in R to explore data objects

---
*Last updated: August 20, 2025*
*Project: Convergent Molecular Signatures in High-Altitude Populations*
