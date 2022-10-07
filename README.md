# CU-Heavy-Ions-Jet-Reco-ML
Exploration of machine learning for heavy ion jet pt reconstruction. Managed by the University of Colorado Boulder Heavy Ions Group.

# Order of Operations
1. Event_Generator_Function_PYTHIA.C
2. Event_Analyzer_ROOT.C (optional, for analysis)
3. Jet_ML_Prep_Function_ROOT.C
4. ML_Estimators_PYTHON.ipynb
5. Jet_ML_Plotter_Function_ROOT.C (optional, for analysis)
6. Extract_Select_Variables_ROOT.C (optional, for analysis)

# Requirements
### Event_Generator_Function_PYTHIA.C
Uses ROOT, PYTHIA 8, FastJet. 
Must compile using Makefile and command `make` in Terminal.

### Event_Analyzer_ROOT.C
Uses ROOT only

### Jet_ML_Prep_Function_ROOT.C
Uses JuPyROOT (ROOT's version of Jupyter Notebook), SciKit-Learn.
Must open using `root --notebook` in Terminal.
