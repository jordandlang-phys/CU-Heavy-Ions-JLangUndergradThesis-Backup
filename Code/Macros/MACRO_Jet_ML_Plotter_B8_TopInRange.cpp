#include "../Scripts_Cpp/Jet_ML_Plotter_ROOT.cpp"
using namespace Jet_Plotter;



void MACRO_Jet_ML_Plotter_B8_TopInRange() {
    
    char dir_master[200];
    sprintf(dir_master, "../../Files/Comparison_Test_4");
    
    char dir_data[200];
    sprintf(dir_data, "%s/Data", dir_master);
    
    char dir_plots[200];
    sprintf(dir_plots, "%s/Plots", dir_master);
    
    char output_file_name[100];
    sprintf(output_file_name, "ML_Results_TopInRange/Test_4GeV_Bins/Train_B8_F12_10_90.root");
    
    char output_file_path[200];
    sprintf(output_file_path, "%s/%s", dir_data, output_file_name);
    
    char ml_weights_array[1][2][100] = {
        {"ML_Results_TopInRange/Test_4GeV_Bins/Train_B8_F12_10_90_LR_Coeffs.csv", "Train_B8_F12_10_90"}
    };

    char ml_results_array[7][2][100] = {
        {"ML_Results_TopInRange/Test_4GeV_Bins/Train_B8_F12_10_90_Test_18_22.csv", "Tree_Train_B8_F12_10_90_Test_18_22"},
        {"ML_Results_TopInRange/Test_4GeV_Bins/Train_B8_F12_10_90_Test_28_32.csv", "Tree_Train_B8_F12_10_90_Test_28_32"},
        {"ML_Results_TopInRange/Test_4GeV_Bins/Train_B8_F12_10_90_Test_38_42.csv", "Tree_Train_B8_F12_10_90_Test_38_42"},
        {"ML_Results_TopInRange/Test_4GeV_Bins/Train_B8_F12_10_90_Test_48_52.csv", "Tree_Train_B8_F12_10_90_Test_48_52"},
        {"ML_Results_TopInRange/Test_4GeV_Bins/Train_B8_F12_10_90_Test_58_62.csv", "Tree_Train_B8_F12_10_90_Test_58_62"},
        {"ML_Results_TopInRange/Test_4GeV_Bins/Train_B8_F12_10_90_Test_68_72.csv", "Tree_Train_B8_F12_10_90_Test_68_72"},
        {"ML_Results_TopInRange/Test_4GeV_Bins/Train_B8_F12_10_90_Test_78_82.csv", "Tree_Train_B8_F12_10_90_Test_78_82"}
    };

    float test_min_max_array[7][2] = {
        {18., 22.}, {28., 32.}, {38., 42.}, {48., 52.}, {58., 62.}, {68., 72.}, {78., 82.}
    };
    
    // Converts ML estimator feature weights into a TTree
    for ( int i = 0 ; i < 7 ; i++ ) {
        Build_Weights_TTree_FromCSV(
            dir_data,
            ml_weights_array[i][0],
            output_file_name,
            ml_weights_array[i][1]
            );
    }

    // Converts ML results into a TTree
    for ( int i = 0 ; i < 7 ; i++ ) {
        Build_Results_TTree_FromCSV(
            dir_data,
            ml_results_array[i][0],
            output_file_name,
            ml_results_array[i][1]
            );
    }

    // Iterates through ML results to output plots
    for ( int i = 0 ; i < 7 ; i++ ) {

        char input_file_name[300];
        sprintf(
            input_file_name,
            "%s/%s",
            dir_data,
            output_file_name);

        char plot_file_name_actual[300];
        sprintf(
            plot_file_name_actual,
            "%s/ML_Results_TopInRange/Test_4GeV_Bins/Plots_Actual/Plot_Train_B8_F12_%.0f_%.0f_Test_%.0f_%.0f",
            dir_data,
            test_min_max_array[i][0],
            test_min_max_array[i][1],
            test_min_max_array[i][0],
            test_min_max_array[i][1]);

        char plot_file_name_delta[300];
        sprintf(
            plot_file_name_delta,
            "%s/ML_Results_TopInRange/Test_4GeV_Bins/Plots_Delta/Plot_Train_B8_F12_%.0f_%.0f_Test_%.0f_%.0f",
            dir_data,
            test_min_max_array[i][0],
            test_min_max_array[i][1],
            test_min_max_array[i][0],
            test_min_max_array[i][1]);

        char plot_title[100];
        sprintf(
            plot_title,
            "Test: p_{T}^{True} in %.0f-%.0f GeV [Train: p_{T}^{True} in 10-90 GeV, No p_{T} Bias, 12 Features]",
            test_min_max_array[i][0],
            test_min_max_array[i][1],
            test_min_max_array[i][0],
            test_min_max_array[i][1]);

        char ml_results[100];
        sprintf(ml_results, "%s", ml_results_array[i][1]);
        float test_pt_min = test_min_max_array[i][0];
        float test_pt_max = test_min_max_array[i][1];

        std::cout << input_file_name << std::endl;
        std::cout << ml_results << std::endl;
        std::cout << plot_file_name_actual << std::endl;
        std::cout << plot_title << std::endl;

        Plot_JetPt_True_Reco_Corr(
            input_file_name,
            ml_results,
            plot_file_name_actual,
            plot_title,
            test_pt_min,
            test_pt_max
        );

        Plot_JetPt_True_Reco_Corr(
            input_file_name,
            ml_results,
            plot_file_name_delta,
            plot_title,
            test_pt_min,
            test_pt_max,
            true
        );

        std::cout << "Distribution Plots Completed" << std::endl;
    
    }
    
    // 20 GeV Bins
    
    char output_file_name_2[100];
    sprintf(output_file_name_2, "ML_Results_TopInRange/Train_20GeV_Bins/Train_B0_F12_10_90.root");
    
    char output_file_path_2[200];
    sprintf(output_file_path_2, "%s/%s", dir_data, output_file_name_2);
    
    char ml_weights_array_2[7][2][100] = {
        {"ML_Results_TopInRange/Train_20GeV_Bins/Train_B8_F12_10_30_LR_Coeffs.csv", "Train_B8_F12_10_30_Test_10_30"},
        {"ML_Results_TopInRange/Train_20GeV_Bins/Train_B8_F12_20_40_LR_Coeffs.csv", "Train_B8_F12_20_40_Test_20_40"},
        {"ML_Results_TopInRange/Train_20GeV_Bins/Train_B8_F12_30_50_LR_Coeffs.csv", "Train_B8_F12_30_50_Test_30_50"},
        {"ML_Results_TopInRange/Train_20GeV_Bins/Train_B8_F12_40_60_LR_Coeffs.csv", "Train_B8_F12_40_60_Test_40_60"},
        {"ML_Results_TopInRange/Train_20GeV_Bins/Train_B8_F12_50_70_LR_Coeffs.csv", "Train_B8_F12_50_70_Test_50_70"},
        {"ML_Results_TopInRange/Train_20GeV_Bins/Train_B8_F12_60_80_LR_Coeffs.csv", "Train_B8_F12_60_80_Test_60_80"},
        {"ML_Results_TopInRange/Train_20GeV_Bins/Train_B8_F12_70_90_LR_Coeffs.csv", "Train_B8_F12_70_90_Test_70_90"}
    };

    char ml_results_array_2[7][2][100] = {
        {"ML_Results_TopInRange/Train_20GeV_Bins/Train_B8_F12_10_30_Test_10_30.csv", "Tree_Train_B8_F12_10_30_Test_10_30"},
        {"ML_Results_TopInRange/Train_20GeV_Bins/Train_B8_F12_20_40_Test_20_40.csv", "Tree_Train_B8_F12_20_40_Test_20_40"},
        {"ML_Results_TopInRange/Train_20GeV_Bins/Train_B8_F12_30_50_Test_30_50.csv", "Tree_Train_B8_F12_30_50_Test_30_50"},
        {"ML_Results_TopInRange/Train_20GeV_Bins/Train_B8_F12_40_60_Test_40_60.csv", "Tree_Train_B8_F12_40_60_Test_40_60"},
        {"ML_Results_TopInRange/Train_20GeV_Bins/Train_B8_F12_50_70_Test_50_70.csv", "Tree_Train_B8_F12_50_70_Test_50_70"},
        {"ML_Results_TopInRange/Train_20GeV_Bins/Train_B8_F12_60_80_Test_60_80.csv", "Tree_Train_B8_F12_60_80_Test_60_80"},
        {"ML_Results_TopInRange/Train_20GeV_Bins/Train_B8_F12_70_90_Test_70_90.csv", "Tree_Train_B8_F12_70_90_Test_70_90"}
    };

    float test_min_max_array_2[7][2] = {
        {10., 30.}, {20., 40.}, {30., 50.}, {40., 60.}, {50., 70.}, {60., 80.}, {70., 90.}
    };

    // Converts ML estimator feature weights into a TTree
    for ( int i = 0 ; i < 7 ; i++ ) {
        Build_Weights_TTree_FromCSV(
            dir_data,
            ml_weights_array_2[i][0],
            output_file_name_2,
            ml_weights_array_2[i][1]
            );
    }

    // Converts ML results into a TTree
    for ( int i = 0 ; i < 7 ; i++ ) {
        Build_Results_TTree_FromCSV(
            dir_data,
            ml_results_array_2[i][0],
            output_file_name_2,
            ml_results_array_2[i][1]
            );
    }

    // Iterates through ML results to output plots
    for ( int i = 0 ; i < 7 ; i++ ) {

        char input_file_name[300];
        sprintf(
            input_file_name,
            "%s/%s",
            dir_data,
            output_file_name_2);

        char plot_file_name_actual[300];
        sprintf(
            plot_file_name_actual,
            "%s/ML_Results_TopInRange/Train_20GeV_Bins/Plots_Actual/Plot_Train_B8_F12_%.0f_%.0f_Test_%.0f_%.0f",
            dir_data,
            test_min_max_array_2[i][0],
            test_min_max_array_2[i][1],
            test_min_max_array_2[i][0],
            test_min_max_array_2[i][1]);

        char plot_file_name_delta[300];
        sprintf(
            plot_file_name_delta,
            "%s/ML_Results_TopInRange/Train_20GeV_Bins/Plots_Delta/Plot_Train_B8_F12_%.0f_%.0f_Test_%.0f_%.0f",
            dir_data,
            test_min_max_array_2[i][0],
            test_min_max_array_2[i][1],
            test_min_max_array_2[i][0],
            test_min_max_array_2[i][1]);

        char plot_title[100];
        sprintf(
            plot_title,
            "Test: p_{T}^{True} in %.0f-%.0f GeV [Train: p_{T}^{True} in %.0f-%.0f GeV, No p_{T} Bias, 12 Features]",
            test_min_max_array[i][0],
            test_min_max_array[i][1],
            test_min_max_array[i][0],
            test_min_max_array[i][1]);

        char ml_results[100];
        sprintf(ml_results, "%s", ml_results_array_2[i][1]);
        float test_pt_min = test_min_max_array_2[i][0];
        float test_pt_max = test_min_max_array_2[i][1];

        std::cout << input_file_name << std::endl;
        std::cout << ml_results << std::endl;
        std::cout << plot_file_name_actual << std::endl;
        std::cout << plot_title << std::endl;

        Plot_JetPt_True_Reco_Corr(
            input_file_name,
            ml_results,
            plot_file_name_actual,
            plot_title,
            test_pt_min,
            test_pt_max
        );

        Plot_JetPt_True_Reco_Corr(
            input_file_name,
            ml_results,
            plot_file_name_delta,
            plot_title,
            test_pt_min,
            test_pt_max,
            true
        );

        std::cout << "Distribution Plots Completed" << std::endl;
    
    }
    
    char input_tree_names_2[7][100] {
        "Train_B8_F12_10_30_Test_10_30",
        "Train_B8_F12_20_40_Test_20_40",
        "Train_B8_F12_30_50_Test_30_50",
        "Train_B8_F12_40_60_Test_40_60",
        "Train_B8_F12_50_70_Test_50_70",
        "Train_B8_F12_60_80_Test_60_80",
        "Train_B8_F12_70_90_Test_70_90"
    };

    char coeff_plot_file_path_2[200];
    sprintf(coeff_plot_file_path_2, "%s/ML_Results_TopInRange/Train_20GeV_Bins/ML_Train_B8_F12_10_90_Coefficients.pdf", dir_data);
    char coeff_plot_file_path_2B[200];
    sprintf(coeff_plot_file_path_2B, "%s/ML_Results_TopInRange/Train_20GeV_Bins/ML_Train_B8_F12_10_90_Coefficients_ShowInt.pdf", dir_data);

    Plot_ML_LR_Coefficients_F12(
        output_file_path_2,
        input_tree_names_2,
        coeff_plot_file_path_2,
        test_min_max_array_2,
        7,
        false
    );
    Plot_ML_LR_Coefficients_F12(
        output_file_path_2,
        input_tree_names_2,
        coeff_plot_file_path_2B,
        test_min_max_array_2,
        7,
        true
    );

    std::cout << "Coefficient Plot Completed" << std::endl;
    
    // 30 GeV Bins
    
    char output_file_name_3[100];
    sprintf(output_file_name_3, "ML_Results_TopInRange/Train_30GeV_Bins/Train_B8_F12_10_90.root");
    
    char output_file_path_3[200];
    sprintf(output_file_path_3, "%s/%s", dir_data, output_file_name_3);
    
    char ml_weights_array_3[6][2][100] = {
        {"ML_Results_TopInRange/Train_30GeV_Bins/Train_B8_F12_10_40_LR_Coeffs.csv", "Train_B8_F12_10_40_Test_10_40"},
        {"ML_Results_TopInRange/Train_30GeV_Bins/Train_B8_F12_20_50_LR_Coeffs.csv", "Train_B8_F12_20_50_Test_20_50"},
        {"ML_Results_TopInRange/Train_30GeV_Bins/Train_B8_F12_30_60_LR_Coeffs.csv", "Train_B8_F12_30_60_Test_30_60"},
        {"ML_Results_TopInRange/Train_30GeV_Bins/Train_B8_F12_40_70_LR_Coeffs.csv", "Train_B8_F12_40_70_Test_40_70"},
        {"ML_Results_TopInRange/Train_30GeV_Bins/Train_B8_F12_50_80_LR_Coeffs.csv", "Train_B8_F12_50_80_Test_50_80"},
        {"ML_Results_TopInRange/Train_30GeV_Bins/Train_B8_F12_60_90_LR_Coeffs.csv", "Train_B8_F12_60_90_Test_60_90"}
    };

    char ml_results_array_3[6][2][100] = {
        {"ML_Results_TopInRange/Train_30GeV_Bins/Train_B8_F12_10_40_Test_10_40.csv", "Tree_Train_B8_F12_10_40_Test_10_40"},
        {"ML_Results_TopInRange/Train_30GeV_Bins/Train_B8_F12_20_50_Test_20_50.csv", "Tree_Train_B8_F12_20_50_Test_20_50"},
        {"ML_Results_TopInRange/Train_30GeV_Bins/Train_B8_F12_30_60_Test_30_60.csv", "Tree_Train_B8_F12_30_60_Test_30_60"},
        {"ML_Results_TopInRange/Train_30GeV_Bins/Train_B8_F12_40_70_Test_40_70.csv", "Tree_Train_B8_F12_40_70_Test_40_70"},
        {"ML_Results_TopInRange/Train_30GeV_Bins/Train_B8_F12_50_80_Test_50_80.csv", "Tree_Train_B8_F12_50_80_Test_50_80"},
        {"ML_Results_TopInRange/Train_30GeV_Bins/Train_B8_F12_60_90_Test_60_90.csv", "Tree_Train_B8_F12_60_90_Test_60_90"}
    };


    float test_min_max_array_3[6][2] = {
        {10., 40.}, {20., 50.}, {30., 60.}, {40., 70.}, {50., 80.}, {60., 90.}
    };
    
    // Converts ML estimator feature weights into a TTree
    for ( int i = 0 ; i < 6 ; i++ ) {
        Build_Weights_TTree_FromCSV(
            dir_data,
            ml_weights_array_3[i][0],
            output_file_name_3,
            ml_weights_array_3[i][1]
            );
    }

    // Converts ML results into a TTree
    for ( int i = 0 ; i < 6 ; i++ ) {
        Build_Results_TTree_FromCSV(
            dir_data,
            ml_results_array_3[i][0],
            output_file_name_3,
            ml_results_array_3[i][1]
            );
    }

    // Iterates through ML results to output plots
    for ( int i = 0 ; i < 6 ; i++ ) {

        char input_file_name[300];
        sprintf(
            input_file_name,
            "%s/%s",
            dir_data,
            output_file_name_3);

        char plot_file_name_actual[300];
        sprintf(
            plot_file_name_actual,
            "%s/ML_Results_TopInRange/Train_30GeV_Bins/Plots_Actual/Plot_Train_B8_F12_%.0f_%.0f_Test_%.0f_%.0f",
            dir_data,
            test_min_max_array_3[i][0],
            test_min_max_array_3[i][1],
            test_min_max_array_3[i][0],
            test_min_max_array_3[i][1]);

        char plot_file_name_delta[300];
        sprintf(
            plot_file_name_delta,
            "%s/ML_Results_TopInRange/Train_30GeV_Bins/Plots_Delta/Plot_Train_B8_F12_%.0f_%.0f_Test_%.0f_%.0f",
            dir_data,
            test_min_max_array_3[i][0],
            test_min_max_array_3[i][1],
            test_min_max_array_3[i][0],
            test_min_max_array_3[i][1]);

        char plot_title[100];
        sprintf(
            plot_title,
            "Test: p_{T}^{True} in %.0f-%.0f GeV [Train: p_{T}^{True} in %.0f-%.0f GeV, No p_{T} Bias, 12 Features]",
            test_min_max_array_3[i][0],
            test_min_max_array_3[i][1],
            test_min_max_array_3[i][0],
            test_min_max_array_3[i][1]);

        char ml_results[100];
        sprintf(ml_results, "%s", ml_results_array_3[i][1]);
        float test_pt_min = test_min_max_array_3[i][0];
        float test_pt_max = test_min_max_array_3[i][1];

        std::cout << input_file_name << std::endl;
        std::cout << ml_results << std::endl;
        std::cout << plot_file_name_actual << std::endl;
        std::cout << plot_title << std::endl;

        Plot_JetPt_True_Reco_Corr(
            input_file_name,
            ml_results,
            plot_file_name_actual,
            plot_title,
            test_pt_min,
            test_pt_max
        );

        Plot_JetPt_True_Reco_Corr(
            input_file_name,
            ml_results,
            plot_file_name_delta,
            plot_title,
            test_pt_min,
            test_pt_max,
            true
        );

        std::cout << "Distribution Plots Completed" << std::endl;

    }
    
    char input_tree_names_3[6][100] {
        "Train_B8_F12_10_40_Test_10_40",
        "Train_B8_F12_20_50_Test_20_50",
        "Train_B8_F12_30_60_Test_30_60",
        "Train_B8_F12_40_70_Test_40_70",
        "Train_B8_F12_50_80_Test_50_80",
        "Train_B8_F12_60_90_Test_60_90"
    };

    char coeff_plot_file_path_3[200];
    sprintf(coeff_plot_file_path_3, "%s/ML_Results_TopInRange/Train_30GeV_Bins/ML_Train_B8_F12_10_90_Coefficients.pdf", dir_data);
    char coeff_plot_file_path_3B[200];
    sprintf(coeff_plot_file_path_3B, "%s/ML_Results_TopInRange/Train_30GeV_Bins/ML_Train_B8_F12_10_90_Coefficients_ShowInt.pdf", dir_data);

    Plot_ML_LR_Coefficients_F12(
        output_file_path_3,
        input_tree_names_3,
        coeff_plot_file_path_3,
        test_min_max_array_3,
        6,
        false
    );
    Plot_ML_LR_Coefficients_F12(
        output_file_path_3,
        input_tree_names_3,
        coeff_plot_file_path_3B,
        test_min_max_array_3,
        6,
        true
    );

    std::cout << "Coefficient Plot Completed" << std::endl;

}

