# Molecular Dynamics Analysis of Protein Stability (2K6M)

## Project Overview
This project investigates the stability of the Human Supervillin Headpiece (PDB: 2K6M) using Molecular Dynamics (MD) simulations. The primary goal is to determine how **temperature influences protein stability** by tracking structural deviations across multiple simulations.

Inspired by [Motono et al. (2022)](https://doi.org/10.1093/bioinformatics/btab538), which proposed a contact frequency-based stability evaluation method for point mutations, this study explores whether the same method can be applied to temperature-dependent stability assessments. Instead of comparing two different mutants, we compare the same protein at **different temperatures** to identify conditions that induce significant motion **without causing unfolding**.

This approach allows us to **pre-screen protein designs computationally**, reducing laboratory screening costs and time.

---

## Simulation Details
- **Protein:** Human Supervillin Headpiece (PDB: 2K6M)  
- **Force Field:** CHARMM27  
- **Water Model:** TIP3P  
- **Ionic Conditions:** 0.15 M NaCl (Na⁺ and Cl⁻ ions)  
- **Temperature Range:** 300K - 370K (in 10K increments)  
- **Simulation Software:** GROMACS 2022  
- **Analysis Tools:** `mdtraj`, `Bio.PDB`, `matplotlib`, `pandas`  

The **trajectories** from different temperatures are analyzed to determine the optimal high-temperature condition where the protein undergoes substantial motion **without unfolding**.

---

## Project Organization
```
├── determine_stability.ipynb        <- Main notebook for analysis
├── environment.yml                  <- Conda environment file
├── README.md                        <- This readme
│
├── code
│   ├── mdcontactcom                 <- Modified scripts from https://gitlab.com/chiemotono/mdcontactcom
│   └── utils.py                     <- Helper functions for analysis
│
├── data
│   ├── params                       <- Simulation parameter files
│   ├── final                        <- Processed data outputs
│   ├── processed                    <- Intermediate transformed data
│   └── raw                          <- Raw input files (e.g., original PDB structure)
│
├── docs                             <- Referenced articles and project report
└── figures                          <- Figures produced in the notebook
```

---

## Reproducing the Analysis

The analysis was carride out on a Linux system. Running the notebook on a Windows system is not recommended.

### 1. Install Dependencies
The analysis requires **GROMACS 2022** and several Python libraries. The easiest way to install them is using **Conda**.

```bash
conda env create -f environment.yml 
conda activate gromacs
```

Alternatively, you can manually install dependencies:
```bash
pip install numpy matplotlib pandas mdtraj nglview
```

### 2. Download Input Structure
```bash
wget -O ./data/raw/2k6m.pdb "https://files.rcsb.org/download/2K6M.pdb"
```

### 3. Run the Notebook
```bash
cd code
jupyter-notebook
```


## References
- Motono, C., et al. (2022). "Contact frequency-based analysis of protein stability in molecular dynamics simulations." *Bioinformatics*. [DOI: 10.1093/bioinformatics/btab538](https://doi.org/10.1093/bioinformatics/btab538)
- GROMACS Documentation: [https://manual.gromacs.org/](https://manual.gromacs.org/)

---

