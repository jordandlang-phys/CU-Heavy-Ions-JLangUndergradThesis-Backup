# Add here if multiple files need to be compiled
all: MACRO_Jet_ML_Generator.cpp

# Compiles Event_Generator_PYTHIA.C regardless of changes
g++ MACRO_Jet_ML_Generator.cpp /Users/jordanlang/Applications/pythia8307/lib/libpythia8.a -o Jet_ML_Generator -I/Users/jordanlang/Applications/pythia8307/include -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC -Wl,-rpath,/Users/jordanlang/Applications/pythia8307/lib -ldl `~/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` -ldl `root-config --libs --cflags`

# NOTES:
# - The paths for 'libpythia8.a', 'pythia8307/include', 'pythia8307/lib', 'fastjet-config' need to be updated for local PYTHIA 8 and FastJet installations.
# - If this doesn't run, copy the g++ line to the terminal and run manually.
