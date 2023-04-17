#include "../Scripts_Cpp/Jet_ML_Flattener_ROOT.cpp"
using namespace Jet_Flattener;

void MACRO_Jet_Flattener() {
//    Make_Flat_Distribution(
//        "../../Files/Thesis_Data/Full_Train_B8_10_90_N500000_ML_Prep.root", // mlprep_file_path
//        "ML_Train_B8_10_90_N500000", // input_tree_name
//        "jet_pt_true", // pt_true_label
//        80,     // n_bins
//        10.,    // jet_pt_min
//        90.,    // jet_pt_max
//        2500   // min_bin_n
//    );
    
    Make_Flat_Distribution(
        "../../Files/Thesis_Data/Full_Test_B8_10_90_N500000_ML_Prep.root", // mlprep_file_path
        "ML_Test_B8_10_90_N500000", // input_tree_name
        "jet_pt_true", // pt_true_label
        80,     // n_bins
        10.,    // jet_pt_min
        90.,    // jet_pt_max
        2500   // min_bin_n
    );

}
