# Add here if multiple files need to be compiled
# all: Event_Generator_PYTHIA.C

# Compiles Event_Generator_PYTHIA.C regardless of changes
g++ Event_Generator_PYTHIA.C /Users/jordanlang/Applications/pythia8307/lib/libpythia8.a -o Event_Generator_PYTHIA -I/Users/jordanlang/Applications/pythia8307/include -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC -Wl,-rpath,/Users/jordanlang/Applications/pythia8307/lib -ldl `~/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` -ldl `root-config --libs --cflags`

# NOTE: The above paths need to be updated for local PYTHIA 8 installation and FastJet installations
