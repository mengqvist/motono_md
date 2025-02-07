
# Project info


In 2022 [Motono et al.](https://doi.org/10.1093/bioinformatics/btab538) published a method for evaluating protein stability of a point mutant by comparing the contact frequency of two MD trajectories. The method reqires two trajectories to compare and in the article these are trajectories from two different variants (differing by one mutation). 

For this project I want to see whether the method can be used to evaluate the stability of a protein. The trajectories I choose to compare are thus ones made with the same protien, but at different temperatures. However, arbitrarily picking a temperature is not ideal. I therefore go through and simulate several temperatures. I want to find a high temperature at which there is significant movement in the protein, without it unfolding. Applying computational methods for evaluating the fitness of protein designs can save considerable time and capital, by serving as a filter before variants are tested in the laboratory. 

In this project I carry out MD simulations of Human Supervillin Headpiece (2K6M) at different temperatures. Two of these will be used to evaluate its stability. The 2K6M structure is a solution NMR structure, which can be used to evaluate the result of the MD simulations.

For further details, plase refer to the document `./docs/final_report.pdf`.


# Project organization
    ├── determine_stability.ipynb        <- Main notebook for analysis
    ├── environment.yml                  <- Specefication of conda environment
    ├── README.md                        <- This readme
    │
    ├── code
    │   ├── mdcontactcom                 <- Folder containing heavily modified scrips originally from https://gitlab.com/chiemotono/mdcontactcom
    │   └── utils.py                     <- Helper functions for notebook
    │
    ├── data
    │   ├── params                       <- Folder holding parameters for the simulations
    │   ├── final                        <- Folder holding output files of processed data
    │   ├── processed                    <- Folder holding intermediate data that has been transformed
    │   └── raw                          <- Folder holding raw unmodified data
    │
    ├── docs                             <- Folder holding the referenced articles as well as the project report
    └── figures                          <- Folder holding figures produced in notebook    



# Reproducing the analysis

Requrements include GROMACS 2022, python modules: numpy, matplotlib, re, nglviewer, md_traj, pandas. For a full list of installed libraries, and their versions, refer to the printout at the beginning of the notebook in the `code` folder.

The most expedient route is to set up the environment using conda (https://docs.conda.io/en/latest/), and activate.
```bash
conda env create -f environment.yml 
conda activate gromacs
```

Download the input structure.
```bash
wget -O ./data/raw/2k6m.pdb "https://files.rcsb.org/download/2K6M.pdb"
```

Move to the code directory and fire up the notebook.
```bash
cd code
jupyter-notebook
```
