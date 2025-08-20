# CMSNetwork: Genome-to-Network Evidence for Convergent Evolution in Tibetans and Andeans

## Project Overview

This project provides systematic network analysis of convergent genetic evolution signals between Tibetan and Andean high-altitude populations. Using an integrated population genetics and network-based framework, we characterize convergent signatures of adaptation to hypoxic stress across multiple biological layers.

**Full Title**: Genome-to-Network Evidence for Convergent Evolution in Tibetans and Andeans Under High-Altitude Stress

**Authors**: Wanjun Gu, Elijah Lawrence, Sergio Baranzini, Tatum Simonson  
**Affiliations**: 
- Department of Neurology, Weill Institute for Neurosciences, University of California, San Francisco, CA
- School of Medicine, Duke University, Durham, NC
- Division of Pulmonary, Critical Care, Sleep Medicine, and Physiology, Department of Medicine, University of California San Diego School of Medicine, La Jolla, CA

## Research Background

High-altitude adaptation represents a striking example of natural selection in humans. While Tibetans have inhabited the Tibetan Plateau for tens of thousands of years and Andean populations have lived at high altitude for millennia, both groups face comparable hypoxic stressors despite their distinct ancestry and demographic history. This suggests the possibility of convergent evolution beyond individual candidate genes like EPAS1.

## Methodology

### Population Genetic Analysis
- **Samples**: 27 Tibetan and 40 Andean individuals
- **Selection Detection**: Composite of Multiple Signals (CMS) test
- **Threshold**: CMS score ≥ 6 for strong positive selection candidates

### Variant-to-Gene Attribution
- **Framework**: Open Targets Genetics variant-to-gene mapping
- **Integration**: Pathogenicity, genomic proximity, GTEx expression data
- **Attribution Score**: ≥ 0.2 for likely effector genes

### Network Analysis
- **Knowledge Graph**: ~1 million nodes, 13 million edges, 18 entity types
- **Enrichment**: Evaluated against null simulations
- **Contextualization**: Literature review and automated retrieval with large-language model agent

## Key Findings

### Convergent Genes (n=16)
**Metabolic & Regulatory**:
- **EP300**, **ACO2**: Metabolic efficiency
- **TBC1D7**: mTOR regulation
- **TEF**: Circadian control

**Stress Response & Vascular**:
- **VNN2**: Oxidative stress defense
- **COL18A1**: Angiogenesis
- **CD86**: Immune regulation

### Anatomical Convergence (n=211 structures)
- **Respiratory**: 29% - Enhanced oxygen extraction
- **Lymphatic**: 21% - Immune surveillance
- **Digestive**: 15% - Nutrient absorption
- **Cerebellar Nervous System**: 12% - Motor coordination

### Cellular Convergence (n=26 cell types)
- **Immune Populations**: 57.7%
- **Vascular-Associated Cells**: 15.4%
- **Cerebellar Neurons**: Motor control
- **Other**: Colonocytes, mesenchymal stem cells, melanophages

## Project Structure

```
CMSNetwork/
├── README.md                           # This file
├── warp.md                            # Project initialization guide
├── .gitignore                         # Version control exclusions
├── manuscript/
│   └── CMS Network Abstract v0.1.docx # Research manuscript
├── data/
│   ├── v2g.csv.gz                     # Variant-to-gene mappings (3.7GB)
│   ├── tibetan_cms_components.csv     # Tibetan CMS network (674MB)
│   ├── andean_cms_components.csv      # Andean CMS network (698MB)
│   ├── tibetan_cms_annotated.csv      # Annotated Tibetan data
│   └── andean_cms_annotated.csv       # Annotated Andean data
├── analysis/
│   ├── determine_network_overlap.R    # Main overlap analysis
│   ├── determine_raw_gene_convergence.R # Gene convergence
│   ├── determine_overlap.R            # General overlap detection
│   ├── process_v2g.ipynb             # V2G data processing
│   └── run_tkoi.R                    # TKOI analysis
├── results/
│   ├── run_v2g0.1/                   # V2G analysis v0.1
│   ├── run_v2g0.2/                   # V2G analysis v0.2
│   ├── network_overlap/              # Network analysis results
│   ├── annotated_network_overlap/    # Annotated results
│   ├── all_tibetan_andean_overlap.xlsx
│   └── tibetan_andean_overlap_significant.xlsx
└── plots/
    ├── network_convergence.png
    └── gene_level_convergence.png
```

## Quick Start

### Prerequisites
```bash
# R packages
install.packages(c("tidyverse", "igraph", "networkD3", "openxlsx"))

# Python packages (for Jupyter notebook)
pip install pandas numpy jupyter
```

### Run Analysis Pipeline
```bash
# 1. Network overlap analysis
Rscript determine_network_overlap.R

# 2. Gene convergence analysis  
Rscript determine_raw_gene_convergence.R

# 3. Process V2G mappings
jupyter notebook process_v2g.ipynb
```

### Explore Results
```bash
# View significant overlaps
open tibetan_andean_overlap_significant.xlsx

# Check convergence plots
open network_convergence.png
open gene_level_convergence.png
```

## Data Files

**Large Files (>100MB - see .gitignore)**:
- `v2g.csv.gz` (3.7GB): Comprehensive variant-to-gene mappings
- `*cms_components.csv` (670-698MB): Population-specific CMS networks
- `*tkoi_result.rda` (~356MB each): TKOI analysis outputs

**Analysis Files**:
- CMS scores and annotations for both populations
- Network overlap calculations and significance tests
- Gene attribution and functional enrichment results

## Key Results Summary

This study provides **genome-wide evidence of convergent adaptation** between Tibetans and Andeans across multiple biological layers:

1. **Multi-systemic Adaptation**: Coordinated selection on metabolic, vascular, immune, and neural systems
2. **Pathway Convergence**: Oxygen utilization, vascular remodeling, metabolic timing, and stress resilience
3. **Cellular Specialization**: Enhanced immune surveillance, vascular adaptation, and neural coordination
4. **Network Architecture**: Integration of population genetics with gene attribution and network propagation

## Citation

```bibtex
@article{gu2024convergent,
  title={Genome-to-Network Evidence for Convergent Evolution in Tibetans and Andeans Under High-Altitude Stress},
  author={Gu, Wanjun and Lawrence, Elijah and Baranzini, Sergio and Simonson, Tatum},
  journal={In Preparation},
  year={2024}
}
```

## Contact

For questions about the analysis or data, please contact the corresponding authors.

---
*This project demonstrates the power of integrating population genetics with network biology to understand convergent human evolution.*
