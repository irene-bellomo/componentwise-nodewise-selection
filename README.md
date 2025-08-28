# componentwise-nodewise-selection
This folder contains the implementation of the Componentwise Nodewise Selection procedure developed and applied in my MSc thesis in Statistics and Data Science at the University of Florence. 
The code, data, and results for the study of **Functional Graphical Models** on both real and simulated datasets, are based on the paper:

**Boxin Zhao, Percy S. Zhai, Y. Samuel Wang, Mladen Kolar. High-dimensional functional graphical model structure learning via neighborhood selection approach.**  
---

## Repository Structure

The repository is organized into five main folders:

### 1. `data`
Contains the real datasets used for the analyses:

- **fMRI ASD vs Control**: fMRI results for autistic individuals and controls.  
- **fMRI ADHD vs Control**: fMRI results for ADHD individuals and controls.  

> Both datasets are taken from Zhao et al.

---

### 2. `DGP`
Contains analyses performed on simulated data:

- Includes data generation and Monte Carlo estimation for the two DGPs presented in the thesis.  
- Provides a comparison of the results across simulations.

---

### 3. `Functions`
Contains all the functions used in the analyses:

- Functions are based on the methods described in Zhao et al.  
- Includes all utilities needed for generating the data, graph estimation, and analysis.

---

### 4. `Real Data Analysis`
Contains the scripts for analyzing the two real datasets (ADHD and ASD):

- Implements **ComponentWise Nodewise Selection** for graph estimation.  
- Compares the estimated graphs with the benchmark methods presented in the literature.  
- Produces the results used for method comparison and evaluation.
---

### 5. `Results`
Contains the results obtained from **FGLasso** and **Neighborhood Selection**:
