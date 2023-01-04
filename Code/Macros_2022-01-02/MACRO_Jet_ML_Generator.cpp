#include "../Scripts_Cpp/Jet_ML_Generator_PYTHIA.cpp"
using namespace Jet_Generator;

/*
    NOTE: This will need to be compiled before running!
    If using command 'make' does not work, you can copy the 'g++ ...' line from the makefile and run manually in the Terminal.
 */
int main() {

    // Defines working directories for data and plots
    char dir_master[200];
    snprintf(dir_master, 200, "../../Files/Demo");
    
    char dir_data[200];
    snprintf(dir_data, 200, "%s/Data", dir_master);
    
    char dir_plots[200];
    snprintf(dir_plots, 200, "%s/Plots", dir_master);
    
    // Creates directories if they do not already exist
    std::__fs::filesystem::create_directories(dir_master);
    std::__fs::filesystem::create_directories(dir_data);
    std::__fs::filesystem::create_directories(dir_plots);
    
    // Generates jets for training
    Jet_Generator_Full(
        "Train",        // output_base_name (Typically just "Train" or "Test")
        dir_data,       // output_directory
        1000,           // event_count [int]
        4.,             // pt_bias_power (Generates events with bias of pT^N for event energies distribution)
        10.,            // jet_pt_min [GeV]
        90.,            // jet_pt_max [GeV]
        0.5,            // jet_eta_max (Reasonable value is detector_eta - fastjet_radius)
        7.5,            // pt_hat_min (Reasonable value is 0.75 * jet_pt_min)
        0.,             // pt_hat_max (If =0. then no upper limit)
        0.4,            // fastjet_radius
        5.,             // fastjet_pt_min (Minimum pT considered by FastJet for a jet)
        0.2,            // fastjet_match_radius (Reasonable value is less than fastjet_radius)
        0,              // softest_jet_index [int] (If =0, tries to match with all PYTHIA jets. If >0, matches only to top N PYTHIA jets)
        1,              // accept_jet_index [int] (If =0, accepts events with any jet in min/max. If >0, only accepts events with one of top N jets in min/max.)
        1800,           // thermal_mean
        200,            // thermal_sigma
        2760.,          // beam_power [GeV] (Simulated beam power in the detector)
        0.9             // detector_eta (Maximum rapidity of the detector)
        );
    
    // End of macro
    return 0;
}
