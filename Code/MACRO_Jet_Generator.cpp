#include "Jet_Generator_PYTHIA.cpp"
using namespace Jet_Generator;
#include "jet_ml_constants.h"
using namespace Jet_ML_Constants;

const bool print_out        = true;
const bool debug            = false;
const int  print_every_x    = 1000;

void Jet_Generator(
    char  output_base_name[100],    // Base/stem of output file name (no file extensions or prefixes).
    char  output_directory[400],    // Path to desired output directory.
    int   event_count,              // Number of collision events to generate. The number of jets will be greater than this.
    float pt_bias_power,            // Bias applied to the particle pT distribution.
    float jet_pt_min = 10.,         // Minimum jet pT to accept. At least 1 jet per event will be greater than this.
    float jet_pt_max = 90.,         // Maximum jet pT to accept. At least 1 jet per event will be less than this.
    float jet_eta_max = 0.5,        // Maximum rapidity to consider.
    float pt_hat_min = 7.5,         // [GeV] ptHatMin value.
    float pt_hat_max = 0.,          // [GeV] ptHatMax value. If 0. then no upper limit.
    float fastjet_radius = 0.4,     // Jet radius used by FastJet.
    float fastjet_pt_min = 5.0,     // [GeV] Minimum pT considered by FastJet for a
    float detector_eta = 0.9,       // [GeV] Simulated beam power in the dector. Default is 2.76
    float beam_power = 2760.        // [GeV] Simulated beam power in the dector. Default is 2.76 TeV for ALICE.
    ) {
    
//    char pythia_file_path[500];
//    sprintf(pythia_file_path, "%s/Pyth_%s_B%.0f_%.0f_%.0f_N%i.root", output_directory, output_base_name, pt_bias_power, jet_pt_min, jet_pt_max, event_count);
//
//    char thermal_file_path[500];
//    sprintf(thermal_file_path, "%s/Ther_%s_B%.0f_%.0f_%.0f_N%i.root", output_directory, output_base_name, pt_bias_power, jet_pt_min, jet_pt_max, event_count);
//
//    char combined_file_path[500];
//    sprintf(combined_file_path, "%s/Comb_%s_B%.0f_%.0f_%.0f_N%i.root", output_directory, output_base_name, pt_bias_power, jet_pt_min, jet_pt_max, event_count);
    
    std::cout << ">>> Generating PYTHIA Events <<<" << std::endl;
    char pythia_file_path[500];
    pythia_file_path = Pythia_Generator(
        output_base_name,
        output_directory,
        event_count,
        pt_bias_power,
        jet_pt_min,
        jet_pt_max,
        jet_eta_max,
        pt_hat_min,
        pt_hat_max,
        fastjet_radius,
        fastjet_pt_min
        );

    std::cout << ">>> Generating Thermal Events <<<" << std::endl;
    char thermal_file_path[500];
    thermal_file_path = Thermal_Generator(
        output_base_name,
        output_directory,
        event_count,
        pt_bias_power,
        jet_pt_min,
        jet_pt_max
        );
    
    std::cout << ">>> Combining PYTHIA and Thermal Events <<<" << std::endl;
    char combined_file_path[500];
    sprintf(combined_file_path, "%s", Combine_Events(
        output_base_name,
        output_directory,
        pythia_file_path,
        thermal_file_path
        ));
    
    std::cout << ">>> Clustering Combined Jets <<<" << std::endl;
    Jet_Clusterer(
        combined_file_path
        );
    
    std::cout << ">>> Preparing ML File <<<" << std::endl;
    Jet_ML_Prep(
        output_base_name,
        output_directory,
        combined_file_path,
        pythia_file_path
        );
        
    std::cout << ">>> Event Generation and Jet Clustering Complete! <<<" << std::endl;
}

int main() {
    // Establishes directories
    char dir_master[200];
    sprintf(dir_master, "../Files/Comparison_Trial2");
    
    char dir_data[200];
    sprintf(dir_data, "%s/Data", dir_master);
    
    char dir_plots[200];
    sprintf(dir_plots, "%s/Plots", dir_master);
    
    // Makes directories if they don't exist
    std::__fs::filesystem::create_directories(dir_master);
    std::__fs::filesystem::create_directories(dir_data);
    std::__fs::filesystem::create_directories(dir_plots);
    
    // Call Jet Generators below:
    
    Jet_Generator(
        "Train",    // output_base_name
        dir_data,   // output_directory
        500000,     // event_count
        8.,         // pt_bias_power
        10.,        // jet_pt_min
        90.         // jet_pt_max
        );
    
//    Jet_Generator(
//        "Test", // output_base_name (pt_bias, pt_min, pt_max, and event_cout will be added to this name in the functions)
//        dir_data, // output_directory
//        500000, // event_count
//        8., // pt_bias_power
//        10., // jet_pt_min
//        90. // jet_pt_max
//        );
    
    return 0;
}
