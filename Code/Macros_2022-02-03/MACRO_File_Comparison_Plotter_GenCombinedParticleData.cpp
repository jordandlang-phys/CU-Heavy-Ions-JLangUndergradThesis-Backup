#include "../Scripts_Cpp/File_Comparison_Plotter.cpp"
using namespace Comparison_Plotter;

void MACRO_File_Comparison_Plotter_GenCombinedParticleData() {
    
    const int num_bins = 5; // Number of branches/variables to plot comparisons for
    char output_directory[400];
    snprintf(output_directory, 400, "../../Files/Comparison_Test_4/Plots/Comparison_Joey_Jordan/GenCombinedParticleData");
    std::__fs::filesystem::create_directories(output_directory);
    
    char file_A_path[400];
    snprintf(file_A_path, 400, "../../Files/Comparison_Test_4/Data/Full_Train_B8_10_90_N500000.root");
//    snprintf(file_A_path, 400, "../../Files/Comparison_Test_4/Data/Full_Demo_B8_10_90_N10000.root");
    char tree_A_name[100];
    snprintf(tree_A_name, 100, "Combined_Particle_Tree");
    char file_A_label[100];
    snprintf(file_A_label, 100, "Jordan's Data");
    char tree_A_branches[num_bins][50] = {
        "particle_n",       "particle_pt",      "particle_eta",     "particle_phi",
        "particle_m"
        };
        
    char file_B_path[400];
//    snprintf(file_B_path, 400, "../../Files/Joey_Data/Data/Pyth_10_90_Train_Trees8.root");
    snprintf(file_B_path, 400, "../../../Jet_Reco_Joey/CU-Heavy-Ions-Jet-Reco-ML/Data/Pyth_10_90_Train_Trees8.root");
    char tree_B_name[100];
    snprintf(tree_B_name, 100, "Combined_Tree");
    char file_B_label[100];
    snprintf(file_B_label, 100, "Joey's Data");
    char tree_B_branches[num_bins][50] = {
        "particle_n",       "particle_pt",      "particle_eta",     "particle_phi",
        "particle_m"
        };
    
    int branch_bins[num_bins] = {
        50, 100, 51, 50,
        50
        };
    float branch_range[num_bins][2] = {
        {0., 3000.},     {0., 100.},   {-1., 1.},      {-3.2, 3.2},
        {0., 1.},
        };
    char branch_x_units[num_bins][50] = {
        "N_{particles}",        "p_{T} [GeV/c^{2}]",    "rapidity",             "azuimuthal angle",
        "mass [GeV]"
        };
    
    File_Comparison_TH1(
        file_A_path,        // file_A_path
        tree_A_name,        // tree_A_name
        file_A_label,       // file_A_label
        tree_A_branches,    // tree_A_branches
        file_B_path,        // file_B_path
        tree_B_name,        // tree_B_name
        file_B_label,       // file_B_label
        tree_B_branches,    // tree_B_branches
        output_directory,   // output_directory
        num_bins,           // num_bins
        branch_bins,        // branch_bins
        branch_range,       // branch_range
        branch_x_units      // branch_x_units
        );
}
