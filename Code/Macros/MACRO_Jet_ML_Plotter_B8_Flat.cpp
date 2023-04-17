#include "../Scripts_Cpp/Jet_ML_Plotter_ROOT.cpp"
using namespace Jet_Plotter;



void MACRO_Jet_ML_Plotter_B8_Flat() {
    
    string dir_master = "../../Files/Thesis_Data/";
    string dir_data   = dir_master + "LR_Coeff_Test/";
    string dir_plots  = dir_data + "Plots/";
    std::__fs::filesystem::create_directories(dir_plots.c_str());
    
    vector<string> dir_train     = {"Train_B8_Flat/"};
    vector<string> dir_test      = {"Test_B8_Flat/"};
    vector<string> label_bias    = {"B8_Flat"};
    vector<string> label_feature = {"F11"};
    
    vector<vector<vector<string>>> traintest_bins = {
        {{"Train_20GeV_Bins"}, {"10","30"}, {"20","40"}, {"30","50"}, {"40","60"}, {"50","70"}, {"60","80"}, {"70","90"}},
        {{"Train_30GeV_Bins"}, {"10","40"}, {"20","50"}, {"30","60"}, {"40","70"}, {"50","80"}, {"60","90"}},
    };
    
    bool build_root = FALSE;
    
    for ( int g=0 ; g < dir_train.size() ; g++ ) {
        for ( int h=0 ; h < dir_test.size() ; h++ ) {
            for ( int i=0 ; i < label_feature.size() ; i++ ) {
                for ( int j=0 ; j < label_bias.size() ; j++ ) {
                    for ( int k=0 ; k < traintest_bins.size() ; k++ ) {
                        string traintest_type = traintest_bins[k][0][0];
                        string root_file_path_str = dir_data + dir_train[g] + dir_test[h] + "Train_" + label_bias[j] + "_" + label_feature[i] + ".root";
                        float  test_min_max_array[20][2];
                        char   input_tree_names[10][100];
                        for ( int l=1 ; l < traintest_bins[k].size() ; l++ ) {
                            string traintest_min = traintest_bins[k][l][0];
                            string traintest_max = traintest_bins[k][l][1];
                            string csv_file_path_str = dir_data + dir_train[g] + dir_test[h] + label_feature[i] + "_" + traintest_type + "/Train_" + label_bias[j] + "_" + label_feature[i] + "_" + traintest_min + "_" + traintest_max + "_LR_Coeffs.csv";
                            string tree_name = "Train_" + label_bias[j] + "_" + label_feature[i] + "_" + traintest_min + "_" + traintest_max + "_LR_Coeffs";
                            snprintf(input_tree_names[l-1], 100, "%s", tree_name.c_str());
                            test_min_max_array[l-1][0] = std::stof(traintest_min);
                            test_min_max_array[l-1][1] = std::stof(traintest_max);
                            if ( build_root ) {
                                char csv_file_path[400];
                                snprintf(csv_file_path, 400, "%s", csv_file_path_str.c_str());
                                char root_file_path[400];
                                snprintf(root_file_path, 400, "%s", root_file_path_str.c_str());
                                std::cout << "\n" << "CSV Path: " << csv_file_path << "\nROOT Path: " << root_file_path << "\nTree Name: " << input_tree_names[l] << std::endl;
                                Build_Coeffs_TTree_FromCSV(
                                    csv_file_path,
                                    root_file_path,
                                    input_tree_names[l-1]
                                );
                            }
                        }
                        for ( int l=1 ; l < traintest_bins[k].size() ; l++ ) {
                            std::cout << input_tree_names[l-1];
                        }
                        for ( int l=1 ; l < traintest_bins[k].size() ; l++ ) {
                            string traintest_min = traintest_bins[k][l][0];
                            string traintest_max = traintest_bins[k][l][1];
                            string tree_name = "Train_" + label_bias[j] + "_" + label_feature[i] + "_" + traintest_min + "_" + traintest_max + "_LR_Coeffs";
                            string plot_file_name_str = dir_plots + "Train_" + label_bias[j] + "_" + label_feature[i] + "_" + traintest_type;
                            char  input_file_path[500];
                            snprintf(input_file_path, 500, "%s", root_file_path_str.c_str());
                            char  plot_file_name[200];
                            snprintf(plot_file_name, 200, "%s.pdf", plot_file_name_str.c_str());
                            char  plot_file_name_yint[100];
                            snprintf(plot_file_name_yint, 200, "%s_yint.pdf", plot_file_name_str.c_str());
                            int   bin_count = traintest_bins[k].size() - 1;
                            std::cout << "\n" << "ROOT Path: " << input_file_path << "\nPlot Name: " << plot_file_name << "\nBin Count: " << bin_count << std::endl;
                            Plot_LR_Coefficients(
                                input_file_path,
                                input_tree_names,
                                plot_file_name,
                                test_min_max_array,
                                bin_count,
                                false
                            );
                            Plot_LR_Coefficients(
                                input_file_path,
                                input_tree_names,
                                plot_file_name_yint,
                                test_min_max_array,
                                bin_count,
                                true
                            );
                            std::cout << "Coefficient Plot Completed" << std::endl;
                        }
                    }
                }
            }
        }
    }
    
    
//
//
//    char dir_master[200];
//    sprintf(dir_master, "../../Files/Thesis_Data/");
//
//    char dir_data[200];
//    sprintf(dir_data, "%sData/", dir_master);
//
//    char dir_plots[200];
//    sprintf(dir_plots, "%sPlots/", dir_master);
//
//    char output_file_name[100];
//    sprintf(output_file_name, "ML_Results_LR_Only/12_Feature/Test_4GeV_Bins/Train_B8_F12_10_90.root");
//
//    char ml_weights_array[1][2][100] = {
//        {"ML_Results_LR_Only/12_Feature/Test_4GeV_Bins/Train_B0_F12_10_90_LR_Coeffs.csv", "Train_B8_F12_10_90"}
//    };
//
//    // 20 GeV Train
//
//    char output_file_name_2[100];
//    sprintf(output_file_name_2, "ML_Results_LR_Only/12_Feature/Train_20GeV_Bins/Train_B8_F12_10_90.root");
//
//    char ml_weights_array_2[7][2][100] = {
//        {"Train_B8_Flat/F11_Train_20GeV_Bins/Train_B8_F12_10_30_LR_Coeffs.csv", "Train_B8_F11_10_30_Test_10_30"},
//        {"ML_Results_LR_Only/12_Feature/Train_20GeV_Bins/Train_B8_F12_20_40_LR_Coeffs.csv", "Train_B8_F12_20_40_Test_20_40"},
//        {"ML_Results_LR_Only/12_Feature/Train_20GeV_Bins/Train_B8_F12_30_50_LR_Coeffs.csv", "Train_B8_F12_30_50_Test_30_50"},
//        {"ML_Results_LR_Only/12_Feature/Train_20GeV_Bins/Train_B8_F12_40_60_LR_Coeffs.csv", "Train_B8_F12_40_60_Test_40_60"},
//        {"ML_Results_LR_Only/12_Feature/Train_20GeV_Bins/Train_B8_F12_50_70_LR_Coeffs.csv", "Train_B8_F12_50_70_Test_50_70"},
//        {"ML_Results_LR_Only/12_Feature/Train_20GeV_Bins/Train_B8_F12_60_80_LR_Coeffs.csv", "Train_B8_F12_60_80_Test_60_80"},
//        {"ML_Results_LR_Only/12_Feature/Train_20GeV_Bins/Train_B8_F12_70_90_LR_Coeffs.csv", "Train_B8_F12_70_90_Test_70_90"}
//    };
//
//    char ml_results_array_2[7][2][100] = {
//        {"ML_Results_LR_Only/12_Feature/Train_20GeV_Bins/Train_B8_F12_10_30_Test_10_30.csv", "Tree_Train_B8_F12_10_30_Test_10_30"},
//        {"ML_Results_LR_Only/12_Feature/Train_20GeV_Bins/Train_B8_F12_20_40_Test_20_40.csv", "Tree_Train_B8_F12_20_40_Test_20_40"},
//        {"ML_Results_LR_Only/12_Feature/Train_20GeV_Bins/Train_B8_F12_30_50_Test_30_50.csv", "Tree_Train_B8_F12_30_50_Test_30_50"},
//        {"ML_Results_LR_Only/12_Feature/Train_20GeV_Bins/Train_B8_F12_40_60_Test_40_60.csv", "Tree_Train_B8_F12_40_60_Test_40_60"},
//        {"ML_Results_LR_Only/12_Feature/Train_20GeV_Bins/Train_B8_F12_50_70_Test_50_70.csv", "Tree_Train_B8_F12_50_70_Test_50_70"},
//        {"ML_Results_LR_Only/12_Feature/Train_20GeV_Bins/Train_B8_F12_60_80_Test_60_80.csv", "Tree_Train_B8_F12_60_80_Test_60_80"},
//        {"ML_Results_LR_Only/12_Feature/Train_20GeV_Bins/Train_B8_F12_70_90_Test_70_90.csv", "Tree_Train_B8_F12_70_90_Test_70_90"}
//    };
//
//    // 30 GeV Train
//
//    char output_file_name_3[100];
//    sprintf(output_file_name_3, "ML_Results_LR_Only/12_Feature/Train_30GeV_Bins/Train_B8_F12_10_90.root");
//
//    char ml_weights_array_3[6][2][100] = {
//        {"ML_Results_LR_Only/12_Feature/Train_30GeV_Bins/Train_B8_F12_10_40_LR_Coeffs.csv", "Train_B8_F12_10_40_Test_10_40"},
//        {"ML_Results_LR_Only/12_Feature/Train_30GeV_Bins/Train_B8_F12_20_50_LR_Coeffs.csv", "Train_B8_F12_20_50_Test_20_50"},
//        {"ML_Results_LR_Only/12_Feature/Train_30GeV_Bins/Train_B8_F12_30_60_LR_Coeffs.csv", "Train_B8_F12_30_60_Test_30_60"},
//        {"ML_Results_LR_Only/12_Feature/Train_30GeV_Bins/Train_B8_F12_40_70_LR_Coeffs.csv", "Train_B8_F12_40_70_Test_40_70"},
//        {"ML_Results_LR_Only/12_Feature/Train_30GeV_Bins/Train_B8_F12_50_80_LR_Coeffs.csv", "Train_B8_F12_50_80_Test_50_80"},
//        {"ML_Results_LR_Only/12_Feature/Train_30GeV_Bins/Train_B8_F12_60_90_LR_Coeffs.csv", "Train_B8_F12_60_90_Test_60_90"}
//    };
//
//    char ml_results_array_3[6][2][100] = {
//        {"ML_Results_LR_Only/12_Feature/Train_30GeV_Bins/Train_B8_F12_10_40_Test_10_40.csv", "Tree_Train_B8_F12_10_40_Test_10_40"},
//        {"ML_Results_LR_Only/12_Feature/Train_30GeV_Bins/Train_B8_F12_20_50_Test_20_50.csv", "Tree_Train_B8_F12_20_50_Test_20_50"},
//        {"ML_Results_LR_Only/12_Feature/Train_30GeV_Bins/Train_B8_F12_30_60_Test_30_60.csv", "Tree_Train_B8_F12_30_60_Test_30_60"},
//        {"ML_Results_LR_Only/12_Feature/Train_30GeV_Bins/Train_B8_F12_40_70_Test_40_70.csv", "Tree_Train_B8_F12_40_70_Test_40_70"},
//        {"ML_Results_LR_Only/12_Feature/Train_30GeV_Bins/Train_B8_F12_50_80_Test_50_80.csv", "Tree_Train_B8_F12_50_80_Test_50_80"},
//        {"ML_Results_LR_Only/12_Feature/Train_30GeV_Bins/Train_B8_F12_60_90_Test_60_90.csv", "Tree_Train_B8_F12_60_90_Test_60_90"}
//    };
//
//    // 20 GeV Bins
//
//    char output_file_path_2[200];
//    sprintf(output_file_path_2, "%s%s", dir_data, output_file_name_2);
//
//    float test_min_max_array_2[7][2] = {
//        {10., 30.}, {20., 40.}, {30., 50.}, {40., 60.}, {50., 70.}, {60., 80.}, {70., 90.}
//    };
//
//    // Converts ML estimator feature weights into a TTree
//    for ( int i = 0 ; i < 7 ; i++ ) {
//        Build_Weights_TTree_FromCSV(
//            dir_data,
//            ml_weights_array_2[i][0],
//            output_file_name_2,
//            ml_weights_array_2[i][1]
//            );
//    }
//
//    // Converts ML results into a TTree
//    for ( int i = 0 ; i < 7 ; i++ ) {
//        Build_Results_TTree_FromCSV(
//            dir_data,
//            ml_results_array_2[i][0],
//            output_file_name_2,
//            ml_results_array_2[i][1]
//            );
//    }
//
//    char input_tree_names_2[7][100] {
//        "Train_B8_Flat_F12_10_30_Test_10_30",
//        "Train_B8_Flat_F12_20_40_Test_20_40",
//        "Train_B8_Flat_F12_30_50_Test_30_50",
//        "Train_B8_Flat_F12_40_60_Test_40_60",
//        "Train_B8_Flat_F12_50_70_Test_50_70",
//        "Train_B8_Flat_F12_60_80_Test_60_80",
//        "Train_B8_Flat_F12_70_90_Test_70_90"
//    };
//
//    char coeff_plot_file_path_2[200];
//    sprintf(coeff_plot_file_path_2, "%%Plots/Train_20GeV_Bins/ML_Train_B8_Flat_F12_10_90_Coefficients.pdf", dir_data);
//    char coeff_plot_file_path_2B[200];
//    sprintf(coeff_plot_file_path_2B, "%%Plots/Train_20GeV_Bins/ML_Train_B8_Flat_F12_10_90_Coefficients_ShowInt.pdf", dir_data);
//
//    // 30 GeV Bins
//
//    char output_file_path_3[200];
//    sprintf(output_file_path_3, "%s%s", dir_data, output_file_name_3);
//
//    float test_min_max_array_3[6][2] = {
//        {10., 40.}, {20., 50.}, {30., 60.}, {40., 70.}, {50., 80.}, {60., 90.}
//    };
//
//    // Converts ML estimator feature weights into a TTree
//    for ( int i = 0 ; i < 6 ; i++ ) {
//        Build_Weights_TTree_FromCSV(
//            dir_data,
//            ml_weights_array_3[i][0],
//            output_file_name_3,
//            ml_weights_array_3[i][1]
//            );
//    }
//
//    // Converts ML results into a TTree
//    for ( int i = 0 ; i < 6 ; i++ ) {
//        Build_Results_TTree_FromCSV(
//            dir_data,
//            ml_results_array_3[i][0],
//            output_file_name_3,
//            ml_results_array_3[i][1]
//            );
//    }
//
//    char input_tree_names_3[6][100] {
//        "Train_Flat_B8_F12_10_40_Test_10_40",
//        "Train_Flat_B8_F12_20_50_Test_20_50",
//        "Train_Flat_B8_F12_30_60_Test_30_60",
//        "Train_Flat_B8_F12_40_70_Test_40_70",
//        "Train_Flat_B8_F12_50_80_Test_50_80",
//        "Train_Flat_B8_F12_60_90_Test_60_90"
//    };
//
//    char coeff_plot_file_path_3[200];
//    sprintf(coeff_plot_file_path_3, "%Plots/Train_30GeV_Bins/ML_Train_B8_Flat_F12_10_90_Coefficients.pdf", dir_data);
//    char coeff_plot_file_path_3B[200];
//    sprintf(coeff_plot_file_path_3B, "%Plots/Train_30GeV_Bins/ML_Train_B8_Flat_F12_10_90_Coefficients_ShowInt.pdf", dir_data);
//
//    Plot_ML_LR_Coefficients_F12(
//        output_file_path_3,
//        input_tree_names_3,
//        coeff_plot_file_path_3,
//        test_min_max_array_3,
//        6,
//        false // show_intercept
//    );
//    Plot_ML_LR_Coefficients_F12(
//        output_file_path_3,
//        input_tree_names_3,
//        coeff_plot_file_path_3B,
//        test_min_max_array_3,
//        6,
//        true // show_intercept
//    );
//
//    std::cout << "Coefficient Plot Completed" << std::endl;

}

