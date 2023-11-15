# Opening the Density-Functional Theory Black Box

<img align='left' src='https://github.com/tjz21/DFT_PIB_Code/blob/main/figures/graphical_abstract.png' width = "393" height = "452"> 

---

> 'It is nice to know that the computer understands the problem. But **I would like to understand it too**.' 
>  
> \- Eugene Wigner

---

This repository contains three Google Colab notebooks that are designed to facilitate understanding of Density-Functional Theory (DFT) through interactive visualizations. Our motivation for developing this software stems from the knowledge deficiency that is often produced from using DFT as a black box in commercial software. By applying DFT to the familiar particle in a box model system employing a real-space grid basis, we hope to have reduced DFT to its fundamental essence fit for pedagogy. Brief instructions for executing the code are provided at the beginning of each notebook and a [problem sheet](https://github.com/tjz21/DFT_PIB_Code/blob/main/DFT_worksheet.pdf) for getting started is attached. The notebooks can be accessed without any installation through Google Colab by simply clicking on the links and signing in with a Google account (offline alternative is provided [here](offline_jupyter/README.md)). Python programming knowledge is not required.
<br>
<br>
<br>



## Notebook 1&ndash;Particle in a 3D Box
<img align="right" src='https://github.com/tjz21/DFT_PIB_Code/blob/main/figures/NB1_wavefunction.png' width = "350" height = "286">
In this notebook, we’ll consider the particle in a three-dimensional box system treated in any undergraduate physical chemistry textbook. High-quality energy level diagrams and isosurface renderings of the wavefunction can be generated from user-specified box lengths. Depicted here is the 321 state of an anthracene-like box of dimensions 16 x 8 x 3 Bohr. 
<br />
<br>
<strong> Click here to open the notebook in Google Colab: </strong> 

<br>

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/tjz21/DFT_PIB_Code/blob/main/notebooks/NB1_3D_PIB.ipynb)

<br>
<br>
<br>
<br>
<br>


## Notebook 2&ndash;PAH Frontier Orbitals
<img align="right" src='https://github.com/tjz21/DFT_PIB_Code/blob/main/figures/NB2_anthracene.png' width = "300" height = "169">
Next, we’ll look at a real chemical system in the form of polycyclic aromatic hydrocarbons (PAHs). We can perform Hartree-Fock/STO-3G calculations to find the shapes and energies of their frontier molecular orbitals, which can make for interesting comparisons with the analogous results from Notebook 1.
<br />
<br>

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/tjz21/DFT_PIB_Code/blob/main/notebooks/NB2_PAH_HF.ipynb)

<br>
<br>

## Notebook 3&ndash;Density-Functional Theory
<img align="right" src='https://github.com/tjz21/DFT_PIB_Code/blob/main/figures/NB3_density.png' width = "295" height = "200">
Finally, we’ll reconsider the system from Notebook 1, but now we’ll turn on electron-electron interaction through the Kohn-Sham potential. We’ll consider each term of the single-particle Hamiltonian and put everything together into a self-consistent field (SCF) DFT calculation. We can then analyze the how the density and eigeneneriges change as a function of SCF iteration number. LDA and PBE are the available exchange-correlation functionals. <br>
<br>
Full theory notebook:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/tjz21/DFT_PIB_Code/blob/main/notebooks/NB3_DFT_PIB.ipynb)

<br>
Abbreviated notebook with the DFT calculator and analysis tools:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/tjz21/DFT_PIB_Code/blob/main/notebooks/NB3_DFT_PIB_calculator.ipynb)

<br>

---

### Citation &nbsp; &nbsp; [![DOI:<10.1021/acs.jchemed.3c00535>](http://img.shields.io/badge/DOI-10.1021/acs.jchemed.3c00535-blue.svg)](<http://dx.doi.org/10.1021/acs.jchemed.3c00535>)
If you would like to cite this work, please refer to the following publication:

> Hirschi, J. S.; Bashirova, D.; Zuehlsdorff, T. J.
> Opening the Density Functional Theory Black Box: A Collection of Pedagogic Jupyter Notebooks.
> *J. Chem. Educ.*
> **2023**,
> *100*, 4496-4503.
> https://doi.org/10.1021/acs.jchemed.3c00535
