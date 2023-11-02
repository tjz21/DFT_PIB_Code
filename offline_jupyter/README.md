# Offline Jupyter Notebooks
Instead of using Google Colab, offline versions of the notebooks have been made available in this directory and can be run through the conventional Jupyter IDE. Instructions for setting up the conda virtual environment are provided here and have been tested on wsl-Ubuntu 22.04 using the Firefox browser. Although NB1 works fine here, **this is not the recommended approach and most of the GUI features in NB2 and NB3 do not display correctly** (because of differences in how ipywidgets runs on Colab vs a local Jupyter Notebook).
## Conda Instructions
1. Create a conda virtual environment 
```
conda create --name DFT_code python=3.10
```
2. Install dependencies to this environment from the requirements.txt file. This takes ~4 minutes to complete.
```sh
conda activate DFT_code
```
```sh
pip install -r requirements.txt
```
3. Launch a Jupyter notebook instance.
```sh
jupyter notebook
```
4. Download the offline versions of the notebooks and open them through the Jupyter IDE. The visualizations should become active after executing the code cells.
