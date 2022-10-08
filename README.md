# CU-Heavy-Ions-Jet-Reco-ML
Exploration of machine learning for heavy ion jet pt reconstruction. Managed by the University of Colorado Boulder Heavy Ions Group led by Jamie Nagle and Dennis Perepelitsa. This project is currently managed by Jordan Lang.

# Order of Operations
### 1. Event_Generator_PYTHIA.C
_Generates p+p events in PYTHIA, generates thermal background from toy model, combines PYTHIA and thermal background, and runs FastJet to do jet clustering on combined and PYTHIA events._
* Uses: ROOT, PYTHIA 8, FastJet
* Edit `Makefile` paths to use your local installations of PYTHIA8 and FastJet
* Run `make` in Terminal

### 2. Event_Analyzer_ROOT.C (optional, for analysis)
_Outputs plots from Event Generator files for analysis._
* Uses: ROOT

### 3. Jet_ML_Prep_ROOT.C
_Prepares .root files of flat TTrees from clustered jets for use with machine learning._
* Uses: ROOT

### 4. ML_Estimators_PYTHON.ipynb
_Runs three machine learning regression estimators: Linear Regression, Random Forest Regression, and Multilayer Perceptron (Shallow Neural Network) Regression. Can adjust input features for training machine learning._
* Uses: Python 3.8, JuPyROOT (ROOT's version of Jupyter Notebook), SciKit-Learn
* Open JuPyRoot using `root --notebook` in Terminal.

### 5. Extract_Select_Variables_ROOT.C (optional, for analysis)
_Takes in .csv files exported from machine learning in Python and converts them to .root files for analysis._
* Uses: ROOT

### 6. Jet_ML_Plotter_ROOT.C (optional, for analysis)
_Outputs plots from machine learning for further analysis. Can show plots combining results of each estimator._
* Uses: ROOT
