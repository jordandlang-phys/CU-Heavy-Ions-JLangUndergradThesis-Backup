#include "../Scripts_Cpp/File_Comparison_Plotter.cpp"
using namespace Comparison_Plotter;

void MACRO_File_Comparison_Plotter_MLPrepData() {
    
    const int num_bins = 21; // Number of branches/variables to plot comparisons for
    char output_directory[400];
    snprintf(output_directory, 400, "../../Files/Comparison_Test_4/Plots/Comparison_Joey_Jordan/MLPrepData");
    std::__fs::filesystem::create_directories(output_directory);
    
    char file_A_path[400];
    snprintf(file_A_path, 400, "../../Files/Comparison_Test_4/Data/Full_Train_B8_10_90_N500000_ML_Prep.root");
    char tree_A_name[100];
    snprintf(tree_A_name, 100, "ML_Train_B8_10_90_N500000");
    char file_A_label[100];
    snprintf(file_A_label, 100, "Jordan's Data");
    char tree_A_branches[num_bins][50] = {
        "jet_pt_raw",       "jet_pt_corr",      "jet_pt_true",      "jet_mass",
        "jet_area",         "jet_const_n",      "const_pt_mean",    "const_pt_median",
        "const_1_pt",       "const_2_pt",       "const_3_pt",       "const_4_pt",
        "const_5_pt",       "const_6_pt",       "const_7_pt",       "const_8_pt",
        "const_9_pt",       "const_10_pt",      "jet_y",            "jet_phi",
        "jet_rho"
        };
        
    char file_B_path[400];
    snprintf(file_B_path, 400, "../../Files/Joey_Data/Data/ML_Prep_10_90_Train8.root");
    char tree_B_name[100];
    snprintf(tree_B_name, 100, "Tree_Tree");
    char file_B_label[100];
    snprintf(file_B_label, 100, "Joey's Data");
    char tree_B_branches[num_bins][50] = {
        "jet_pt_raw",       "jet_pt_corr",      "jet_pt_true_pythia", "jet_mass",
        "jet_area",         "jet_const_n",      "const_pt_mean",    "const_pt_median",
        "const_1_pt",       "const_2_pt",       "const_3_pt",       "const_4_pt",
        "const_5_pt",       "const_6_pt",       "const_7_pt",       "const_8_pt",
        "const_9_pt",       "const_10_pt",      "jet_y",            "jet_phi",
        "jet_rho"
        };
    
    int branch_bins[num_bins] = {
        50, 50, 50, 50,
        50, 70, 50, 50,
        50, 50, 50, 50,
        50, 50, 50, 50,
        50, 50, 50, 50,
        50
        };
    float branch_range[num_bins][2] = {
        {0., 200.},     {-20., 180.},   {0., 100.},     {0., 60.},
        {0., 1.},       {20., 160.},    {0., 5.},       {0., 5.},
        {0., 50.},      {0., 35.},      {0., 20.},      {0., 16.},
        {0., 14.},      {0., 12.},      {0., 10.},      {0., 10.},
        {0., 10.},      {0., 10.},      {-1., 1.},      {0, 6.4},
        {40., 180.}
        };
    char branch_x_units[num_bins][50] = {
        "p_{T} [GeV/c^{2}]",    "p_{T} [GeV/c^{2}]",    "p_{T} [GeV/c^{2}]",    "mass [GeV]",
        "area [units^{2}]",     "N_{constituents}",     "p_{T} [GeV/c^{2}]",    "p_{T} [GeV/c^{2}]",
        "p_{T} [GeV/c^{2}]",    "p_{T} [GeV/c^{2}]",    "p_{T} [GeV/c^{2}]",    "p_{T} [GeV/c^{2}]",
        "p_{T} [GeV/c^{2}]",    "p_{T} [GeV/c^{2}]",    "p_{T} [GeV/c^{2}]",    "p_{T} [GeV/c^{2}]",
        "p_{T} [GeV/c^{2}]",    "p_{T} [GeV/c^{2}]",    "rapidity",             "azuimuthal angle",
        "p_{T}/area [GeV/c^{2}/units^{2}]"
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
