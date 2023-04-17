#include "../Scripts_Cpp/Jet_ML_Plotter_ROOT.cpp"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace Jet_Plotter;



void Convert_CSV_to_ROOT(
    bool remake_root = TRUE
    ) {
    // ----- CHANGE THIS PART -----
    string dir_master   = "../../Files/Thesis_Data/";
    vector<string> label_trainbias_array    = {"B8_Flat"};
    vector<string> label_testbias_array     = {"B4", "B0"};
    vector<string> label_features_array     = {"F11"};
    vector<string> label_ptminmax           = {"10", "90"};
    vector<string> label_mltype_array       = {"LR", "RF", "MLP"};
    vector<string> label_coeffs_array       = {"LR_Coeffs", "RF_Coeffs"};
//    vector<string> label_trainbias_array    = {"B8_Flat_T2", "B8_Flat", "B8", "B4", "B0"};
//    vector<string> label_testbias_array     = {"B8_Flat", "B4", "B0"};
//    vector<string> label_features_array     = {"F12", "F11", "F3", "F1"};
//    vector<string> label_ptminmax           = {"10", "90"};
//    vector<string> label_mltype_array       = {"LR", "RF", "MLP"};
//    vector<string> label_coeffs_array       = {"LR_Coeffs", "RF_Coeffs"};
    vector<vector<vector<string>>> test_bin_array = {
        {{"Test_4GeV_Bins"}, {"18","22"}, {"28","32"}, {"38","42"}, {"48","52"}, {"58","62"}, {"68","72"}, {"78","82"}},
//        {{"Test_10GeV_Bins"}, {"10","20"}, {"20","30"}, {"30","40"}, {"40","50"}, {"50","60"}, {"60","70"}, {"70","80"}, {"80","90"}},
        {{"Test_Centered_Wide_Bins"}, {"10","90"}, {"20","80"}, {"30","70"}, {"40","60"}}
    };
    vector<vector<vector<string>>> traintest_bin_array = {
//        {{"Train_20GeV_Bins"}, {"10","30"}, {"20","40"}, {"30","50"}, {"40","60"}, {"50","70"}, {"60","80"}, {"70","90"}},
//        {{"Train_30GeV_Bins"}, {"10","40"}, {"20","50"}, {"30","60"}, {"40","70"}, {"50","80"}, {"60","90"}}
//        {{"Train_Centered_Test_40_60","40","60"}, {"10","40"}, {"20","50"}, {"30","60"}, {"40","70"}, {"50","80"}, {"60","90"}}
    };
    
    // ----- DO NOT CHANGE BELOW THIS!!! -----
    
    // Iterate over each training bias
    for ( int i=0 ; i < label_trainbias_array.size() ; i++ ){
        // Iterate through test biases
        for ( int p=0 ; p < label_testbias_array.size() ; p++ ){
            string dir_plot_csv = dir_master + "CSV_Plot_Info/";
            std::__fs::filesystem::create_directories(dir_plot_csv.c_str());
            std::ofstream plot_csv_data;
            string plot_csv_data_path = dir_plot_csv + "Train_" + label_trainbias_array[i] + "_Test_" + label_testbias_array[p] + "_Plot_Info.csv";
            plot_csv_data.open(plot_csv_data_path.c_str());
            plot_csv_data << "root_file_path, train_ptbias, test_ptbias, train_features, train_test_type, train_pt_min, train_pt_max, test_pt_min, test_pt_max, lr_coeffs_tree, rf_coeffs_tree, lr_results_tree, rf_results_tree, mlp_results_tree, plot_delta, plot_lr_coeffs, plot_comparison,\n";
                
            // Iterate through features
            for ( int j=0 ; j < label_features_array.size() ; j++ ){
                string dir_train = dir_master + "Train_" + label_trainbias_array[i] + "/Test_" + label_testbias_array[p] + "/";
                string dir_feature = dir_train + label_features_array[j] + "/";
                string dir_root = dir_train + "Root/";
                std::__fs::filesystem::create_directories(dir_root.c_str());
            
                // Iterate through test bin sets
                for ( int k=0 ; k < test_bin_array.size() ; k++ ){
                    string dir_csv   = dir_train + label_features_array[j] + "_" + test_bin_array[k][0][0] + "/";
                    string root_name = "Train_" + label_trainbias_array[i] + "_" + label_features_array[j] + "_" + label_ptminmax[0] + "_" + label_ptminmax[1];
                    string root_path = dir_root + root_name + "_" + test_bin_array[k][0][0] + ".root";
                    //std::cout << "ROOT File Path: " << root_path.c_str() << std::endl;
                    
                    // Make coefficient ROOT files and trees
                    string lr_coeffs_tree = "";
                    string rf_coeffs_tree = "";
                    for ( int l=0 ; l < label_coeffs_array.size() ; l++ ){
                        string csv_coeff_path   = dir_master + "Train_" + label_trainbias_array[i] + "/" + root_name + "_" + label_coeffs_array[l] + ".csv";
                        //std::cout << label_coeffs_array[l] << " File Path: " << csv_coeff_path.c_str() << std::endl;
                        string coeff_tree_name  = label_coeffs_array[l];
                        if (l==0) lr_coeffs_tree = coeff_tree_name;
                        if (l==1) rf_coeffs_tree = coeff_tree_name;
                        char input_file_path[400];
                        snprintf(input_file_path, 400, "%s", csv_coeff_path.c_str());
                        char output_file_path[400];
                        snprintf(output_file_path, 400, "%s", root_path.c_str());
                        char output_tree_name[100];
                        snprintf(output_tree_name, 100, "%s", coeff_tree_name.c_str());
                        if ( remake_root ) {
                            Build_Coeffs_TTree_FromCSV(
                                input_file_path,
                                output_file_path,
                                output_tree_name
                            );
                        }
                    }
                    // Make results ROOT files and trees, iterate through test bins
                    for ( int l=1 ; l < test_bin_array[k].size() ; l++ ){
                        string lr_results_tree  = "";
                        string rf_results_tree  = "";
                        string mlp_results_tree = "";
                        // Iterate through ML types
                        for (  int m=0 ; m < label_mltype_array.size() ; m++ ){
                            string tree_name = "Test_" + test_bin_array[k][l][0] + "_" + test_bin_array[k][l][1] + "_" + label_mltype_array[m];
                            if (m==0) lr_results_tree  = tree_name;
                            if (m==1) rf_results_tree  = tree_name;
                            if (m==2) mlp_results_tree = tree_name;
                            string csv_path  = dir_csv + "Train_" + label_trainbias_array[i] + "_" + label_features_array[j] + "_" + label_ptminmax[0] + "_" + label_ptminmax[1] + "_" + tree_name + ".csv";
                            //std::cout << label_features_array[j] << "_" << test_bin_array[k][0][0] << " File Path: " << csv_path.c_str() << std::endl;
                            char input_file_path[400];
                            snprintf(input_file_path, 400, "%s", csv_path.c_str());
                            char output_file_path[400];
                            snprintf(output_file_path, 400, "%s", root_path.c_str());
                            char output_tree_name[100];
                            snprintf(output_tree_name, 100, "%s", tree_name.c_str());
                            if ( remake_root ) {
                                Build_Results_TTree_FromCSV(
                                    input_file_path,
                                    output_file_path,
                                    output_tree_name
                                );
                            }
                        }
                        // Write info to CSV for each set of plots
                        plot_csv_data << (root_path + ",").c_str();
                        plot_csv_data << (label_trainbias_array[i] + ",").c_str();
                        plot_csv_data << (label_testbias_array[p] + ",").c_str();
                        plot_csv_data << (label_features_array[j] + ",").c_str();
                        plot_csv_data << (test_bin_array[k][0][0] + ",").c_str();
                        plot_csv_data << (label_ptminmax[0] + ",").c_str(); // train_pt_min
                        plot_csv_data << (label_ptminmax[1] + ",").c_str(); // train_pt_max
                        plot_csv_data << (test_bin_array[k][l][0] + ",").c_str(); // test_pt_min
                        plot_csv_data << (test_bin_array[k][l][1] + ",").c_str(); // test_pt_max
                        plot_csv_data << (lr_coeffs_tree + ",").c_str();
                        plot_csv_data << (rf_coeffs_tree + ",").c_str();
                        plot_csv_data << (lr_results_tree + ",").c_str();
                        plot_csv_data << (rf_results_tree + ",").c_str();
                        plot_csv_data << (mlp_results_tree + ",").c_str();
                        plot_csv_data << "true,"; // plot_delta
                        plot_csv_data << "false,"; // plot_lr_coeffs
                        plot_csv_data << "true,"; // plot_comparison
                        plot_csv_data << "\n";
                    }
                }
                for ( int k=0 ; k < traintest_bin_array.size() ; k++ ){
                    string dir_csv   = dir_train + label_features_array[j] + "_" + traintest_bin_array[k][0][0] + "/";
                    string root_path = dir_root + "Train_" + label_trainbias_array[i] + "_" + label_features_array[j] + "_" + label_ptminmax[0] + "_" + label_ptminmax[1] + "_" + traintest_bin_array[k][0][0] + ".root";
                    std::cout << "ROOT File Path: " << root_path.c_str() << std::endl;
                    for ( int l=1 ; l < traintest_bin_array[k].size() ; l++ ){
                        bool fixed_test = false;
                        if ( traintest_bin_array[k][0].size() == 3 ) fixed_test = true;
                        string train_pt_min = traintest_bin_array[k][l][0];
                        string train_pt_max = traintest_bin_array[k][l][1];
                        string test_pt_min  = traintest_bin_array[k][l][0];
                        string test_pt_max  = traintest_bin_array[k][l][1];
                        if ( fixed_test ) {
                            test_pt_min = traintest_bin_array[k][0][1];
                            test_pt_max = traintest_bin_array[k][0][2];
                        }
                        string csv_base = "Train_" + label_trainbias_array[i] + "_" + label_features_array[j] + "_" + traintest_bin_array[k][l][0] + "_" + traintest_bin_array[k][l][1];
                        string csv_results_path = dir_csv + csv_base + "_Test_" + test_pt_min + "_" + test_pt_max + ".csv";
                        string csv_coeff_path   = dir_csv + csv_base + "_LR_Coeffs.csv";
                        string results_tree = "Test_" + test_pt_min + "_" + test_pt_max + "_Results";
                        string coeff_tree   = "Test_" + test_pt_min + "_" + test_pt_max + "_LR_Coeffs";
                        if ( fixed_test ) {
                            string results_tree = "Train_" + train_pt_min + "_" + train_pt_max + "_Test_" + test_pt_min + "_" + test_pt_max + "_Results";
                            string coeff_tree   = "Train_" + train_pt_min + "_" + train_pt_max + "_Test_" + test_pt_min + "_" + test_pt_max + "_LR_Coeffs";
                        }
                        char output_file_path[400];
                        snprintf(output_file_path, 400, "%s", root_path.c_str());
                        char results_file_path[400];
                        snprintf(results_file_path, 400, "%s", csv_results_path.c_str());
                        char results_tree_name[100];
                        snprintf(results_tree_name, 100, "%s", results_tree.c_str());
                        char coeff_file_path[400];
                        snprintf(coeff_file_path, 400, "%s", csv_coeff_path.c_str());
                        char coeff_tree_name[100];
                        snprintf(coeff_tree_name, 100, "%s", coeff_tree.c_str());
                        if ( remake_root ) {
                            Build_Results_TTree_FromCSV(
                                results_file_path,
                                output_file_path,
                                results_tree_name
                            );
                            Build_Coeffs_TTree_FromCSV(
                                coeff_file_path,
                                output_file_path,
                                coeff_tree_name
                            );
                        }
                        plot_csv_data << (root_path + ",").c_str();
                        plot_csv_data << (label_trainbias_array[i] + ",").c_str();
                        plot_csv_data << (label_testbias_array[p] + ",").c_str();
                        plot_csv_data << (label_features_array[j] + ",").c_str();
                        plot_csv_data << (traintest_bin_array[k][0][0] + ",").c_str();
                        plot_csv_data << (train_pt_min + ",").c_str(); // train_pt_min
                        plot_csv_data << (train_pt_max + ",").c_str(); // train_pt_max
                        plot_csv_data << (test_pt_min + ",").c_str();  // test_pt_min
                        plot_csv_data << (test_pt_max + ",").c_str();  // test_pt_max
                        plot_csv_data << (coeff_tree + ",").c_str();
                        plot_csv_data << ",";
                        plot_csv_data << (results_tree + ",").c_str();
                        plot_csv_data << ",";
                        plot_csv_data << ",";
                        plot_csv_data << "true,"; // plot_delta
                        plot_csv_data << "true,"; // plot_lr_coeffs
                        plot_csv_data << "false,"; // plot_comparison
                        plot_csv_data << "\n";
                    }
                }
            }
            plot_csv_data.close();
        }
    }
}

void Plot_Delta(
    string input_file_path_str,
    string results_tree_name_str,
    string plot_file_dir_str,
    string plot_title_str,
    string test_pt_min_str,
    string test_pt_max_str
    ){
    
    string plot_dir_actual = plot_file_dir_str + "pT_Actual/";
    string plot_dir_delta  = plot_file_dir_str + "pT_Delta/";
    std::__fs::filesystem::create_directories(plot_dir_actual.c_str());
    std::__fs::filesystem::create_directories(plot_dir_delta.c_str());
    
    char  input_file_path[400];
    snprintf(input_file_path, 400, "%s", input_file_path_str.c_str());
    char  results_tree_name[200];
    snprintf(results_tree_name, 200, "%s", results_tree_name_str.c_str());
    char  plot_file_name_actual[400];
    snprintf(plot_file_name_actual, 400, "%s%s", plot_dir_actual.c_str(), plot_title_str.c_str());
    char  plot_file_name_delta[400];
    snprintf(plot_file_name_delta, 400, "%s%s", plot_dir_delta.c_str(), plot_title_str.c_str());
    char  plot_title[200];
    snprintf(plot_title, 200, "%s", plot_title_str.c_str());
    float test_pt_min = std::stof(test_pt_min_str.c_str());
    float test_pt_max = std::stof(test_pt_max_str.c_str());
    Plot_JetPt_ML_Corr_True_Difference(
        input_file_path,
        results_tree_name,
        plot_file_name_actual,
        plot_title,
        test_pt_min,
        test_pt_max,
        false // use_delta
    );
    Plot_JetPt_ML_Corr_True_Difference(
        input_file_path,
        results_tree_name,
        plot_file_name_delta,
        plot_title,
        test_pt_min,
        test_pt_max,
        true // use_delta
    );
}



void Plot_Comparison(
    string input_file_path_str,
    string plot_file_dir_str,
    string plot_title_str,
    string lr_tree_name_str,
    string rf_tree_name_str,
    string mlp_tree_name_str,
    string label_ptbias_str,
    string label_feature_str,
    string coeff_tree_name_str,
    string coeff_label_str,
    float train_pt_min,
    float train_pt_max,
    float test_pt_min,
    float test_pt_max
    ) {
    char  input_file_path[400];
    snprintf(input_file_path, 400, "%s", input_file_path_str.c_str());
    char  plot_file_dir[400];
    snprintf(plot_file_dir, 400, "%s", plot_file_dir_str.c_str());
    char  plot_title[100];
    snprintf(plot_title, 100, "%s", plot_title_str.c_str());
    char  lr_tree_name[100];
    snprintf(lr_tree_name, 100, "%s", lr_tree_name_str.c_str());
    char  rf_tree_name[100];
    snprintf(rf_tree_name, 100, "%s", rf_tree_name_str.c_str());
    char  mlp_tree_name[100];
    snprintf(mlp_tree_name, 100, "%s", mlp_tree_name_str.c_str());
    char  label_ptbias[100];
    snprintf(label_ptbias, 100, "%s", label_ptbias_str.c_str());
    char  label_feature[100];
    snprintf(label_feature, 100, "%s", label_feature_str.c_str());
    char  coeff_tree_name[100];
    snprintf(coeff_tree_name, 100, "%s", coeff_tree_name_str.c_str());
    char  coeff_label[100];
    snprintf(coeff_label, 100, "%s", coeff_label_str.c_str());
    
    Plot_JetPt_ML_Comparison(
        input_file_path,
        plot_file_dir,
        plot_title,
        lr_tree_name,
        rf_tree_name,
        mlp_tree_name,
        label_ptbias,
        label_feature,
        coeff_tree_name,
        coeff_label,
        train_pt_min,
        train_pt_max,
        test_pt_min,
        test_pt_max,
        true, // normalize
        true, // show_legend
        true, // show_widths
        false, // show_coeffs
        false // show_paper_plots
    );
}



void MACRO_Jet_ML_Plotter_ALL() {
    
    // Converts CSV files to ROOT files and exports a CSV of plots
    // Uncomment to run. Once run, this can be skipped
//    Convert_CSV_to_ROOT(
//        TRUE // remake_root
//        );
//    return;
    
    string dir_master   = "../../Files/Thesis_Data/";
    string dir_plots    = dir_master + "Plots/";
    string dir_plot_instr = dir_master + "CSV_Plot_Info/";
    std::__fs::filesystem::create_directories(dir_plots.c_str());
    // All training biases, only testing on B8_Flat
    vector<string> plot_instr_csv_array = {
        "Train_B0_Test_B8_Flat_Plot_Info.csv",
        "Train_B4_Test_B8_Flat_Plot_Info.csv",
        "Train_B8_Flat_Test_B8_Flat_Plot_Info.csv",
    };
    // Only training on B_Flat, testing on B0 and B4
//    vector<string> plot_instr_csv_array = {
//        "Train_B8_Flat_T2_Test_B8_Flat_Plot_Info.csv"
//        "Train_B8_Flat_Test_B0_Plot_Info.csv",
//        "Train_B8_Flat_Test_B4_Plot_Info.csv"
//    };

    
    for ( int i=0 ; i < plot_instr_csv_array.size() ; i++ ) {
        string plot_instr_path = dir_plot_instr + plot_instr_csv_array[i];
        ifstream plot_csv_file;
        plot_csv_file.open(plot_instr_path.c_str(), ios::in);
        bool skip_header = true;
        while (!plot_csv_file.eof()) {
            string plot_data_row;
            vector<string> plot_data;
            std::getline(plot_csv_file, plot_data_row);
            stringstream plot_data_row_ss(plot_data_row);
            if (skip_header) {
                skip_header = false;
                continue;
            }
            int  cell_counter = 0;
            bool skip_row = false;
            while (plot_data_row_ss.good()) {
                string cell;
                getline(plot_data_row_ss, cell, ',');
                if ( cell.empty() && cell_counter < 8 ) {
                    skip_row = true;
                    break;
                }
                if ( cell.empty() ) cell = "none";
                plot_data.push_back(cell);
                cell_counter++;
            }
            if ( skip_row ) continue;
            
            string root_path        = plot_data[0];
            string label_ptbias     = plot_data[1];
            string label_testbias   = plot_data[2];
            string label_feature    = plot_data[3];
            string label_test_type  = plot_data[4];
            string train_pt_min     = plot_data[5];
            string train_pt_max     = plot_data[6];
            string test_pt_min      = plot_data[7];
            string test_pt_max      = plot_data[8];
            string lr_coeffs_tree   = plot_data[9];
            string rf_coeffs_tree   = plot_data[10];
            string lr_results_tree  = plot_data[11];
            string rf_results_tree  = plot_data[12];
            string mlp_results_tree = plot_data[13];
            bool   plot_delta       = false;
            bool   plot_lr_coeffs   = false;
            bool   plot_comparison  = false;
            if ( plot_data[14].compare("true") == 0 ) plot_delta      = true;
            if ( plot_data[15].compare("true") == 0 ) plot_lr_coeffs  = true;
            if ( plot_data[16].compare("true") == 0 ) plot_comparison = true;
            
            std::cout << root_path +", "+ label_ptbias +", "+ label_feature +", "+ label_test_type +", "+ train_pt_min +", "+ train_pt_max +", "+ test_pt_min +", "+ test_pt_max +", "+ lr_coeffs_tree +", "+ rf_coeffs_tree +", "+ lr_results_tree +", "+ rf_results_tree +", "+ mlp_results_tree +", " << plot_delta << ", " << plot_lr_coeffs << ", " << plot_comparison << std::endl;
            
            string dir_plot_type = dir_plots + "Train_" + label_ptbias + "_" + label_feature + "_Test_" + label_testbias + "_" + label_test_type + "/";
            std::__fs::filesystem::create_directories(dir_plot_type.c_str());
            
            if (false){ //( plot_delta ) {
                string plot_delta_title_base = label_ptbias + "_" + label_feature + "_" + label_test_type + "_Train_" + train_pt_min + "_" + train_pt_max + "_Test_" + label_testbias + "_" + test_pt_min + "_" + test_pt_max;
                if (lr_results_tree.compare("none")) Plot_Delta(
                    root_path,                      // input_file_path_str
                    lr_results_tree,                // results_tree_name_str
                    dir_plot_type,                  // plot_file_dir_str,
                    "LR_" + plot_delta_title_base,  // plot_title_str
                    test_pt_min,                    // test_pt_min
                    test_pt_max                     // test_pt_max
                );
                if (rf_results_tree.compare("none")) Plot_Delta(
                    root_path,                      // input_file_path_str
                    rf_results_tree,                // results_tree_name_str
                    dir_plot_type,                  // plot_file_name_str,
                    "RF_" + plot_delta_title_base,  // plot_title_str
                    test_pt_min,                    // test_pt_min
                    test_pt_max                     // test_pt_max
                );
                if (mlp_results_tree.compare("none")) Plot_Delta(
                    root_path,                      // input_file_path_str
                    mlp_results_tree,               // results_tree_name_str
                    dir_plot_type,                  // plot_file_name_str,
                    "MLP_" + plot_delta_title_base, // plot_title_str
                    test_pt_min,                    // test_pt_min
                    test_pt_max                     // test_pt_max
                );
            }
            if ( plot_comparison  ) {
                string plot_comparison_title = label_ptbias + "_" + label_feature + "_" + label_test_type + "_Train_" + train_pt_min + "_" + train_pt_max + "_Test_" + label_testbias + "_" + test_pt_min + "_" + test_pt_max;
                Plot_Comparison(
                    root_path,
                    dir_plot_type,
                    ("Comparison_" + plot_comparison_title),
                    lr_results_tree,
                    rf_results_tree,
                    mlp_results_tree,
                    label_ptbias,
                    label_feature,
                    rf_coeffs_tree,
                    "Feature Importance",
                    std::stof(train_pt_min),
                    std::stof(train_pt_max),
                    std::stof(test_pt_min),
                    std::stof(test_pt_max)
                );
            }
        }
    }
}

