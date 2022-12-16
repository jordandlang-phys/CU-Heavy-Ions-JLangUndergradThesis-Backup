#include "TFile.h"
#include "TTree.h"
#include "_header.h"
#include <fstream>
#include <vector>
#include <sstream>
using namespace Project_Constants;

void Build_Results_TTree_FromCSV(
    char file_dir[200],         // file directory path
    char input_file_name[200],  // input file name
    char output_file_name[200], // output file name
    char output_tree_name[200]  // name for TTree in output file
    ) {
    
    char input_file_path[200];
    sprintf(input_file_path, "%s/%s", file_dir, input_file_name);
    
    fstream csv_file;
    csv_file.open( input_file_path, ios::in );
    
    std::cout << "Accessed input file." << std::endl;
    
    char output_file_path[200];
    sprintf(output_file_path, "%s/%s", file_dir, output_file_name);
    
    TFile* output_file = new TFile(output_file_path, "UPDATE");
    
    TTree* output_tree = new TTree(output_tree_name, "TTree of results from machine learning testing.");
    std::cout << "Created output file." << std::endl;
    
    output_tree->ReadFile(input_file_path, "", ',');
    
    output_tree->Write("", TObject::kOverwrite);
    
    delete output_tree;
    
    output_file->Close();
    delete output_file;
    
    csv_file.close();
}



void Build_Weights_TTree_FromCSV(
    char* file_dir,
    char* input_file_name,
    char* output_file_name,
    char* output_tree_name
    ) {
    
    char input_file_path[200];
    sprintf(input_file_path, "%s/%s", file_dir, input_file_name);
    
    fstream csv_file;
    csv_file.open( input_file_path, ios::in );
    
    std::cout << "Accessed input file." << std::endl;
    
    char output_file_path[200];
    sprintf(output_file_path, "%s/%s", file_dir, output_file_name);
    
    TFile* output_file = new TFile(output_file_path, "UPDATE");
    
    TTree* output_tree = new TTree(output_tree_name, "TTree of feature weights for given ML method.");
    std::cout << "Created output file." << std::endl;
    
    output_tree->ReadFile(input_file_path, "", ',');
    
    output_tree->Write("", TObject::kOverwrite);
    
    delete output_tree;
    
    output_file->Close();
    delete output_file;
    
    csv_file.close();
}



void Plot_JetPt_True_Reco_Corr(
    char input_file_name[200],      // Full file name with path for input .root file
    char results_tree_name[200],    // Name of TTree to get ML data from
    char plot_file_name[200],       // Full output file name path WITH NO FILE EXTENSION
    char plot_title[200],           // Title to be printed on plot
    float test_pt_min,
    float test_pt_max,
    bool use_delta = false,
    bool use_normalized = true
    ) {
    
    TFile* input_file = new TFile(input_file_name, "READ");
    TTree* tree_results = (TTree*) input_file->Get(results_tree_name);
    std::cout << "Input file accessed." << std::endl;
    
    float jet_pt_true;
    float jet_pt_reco;
    float jet_pt_corr;
    
    tree_results->SetBranchAddress("jet_pt_true", &jet_pt_true);
    tree_results->SetBranchAddress("jet_pt_reco", &jet_pt_reco);
    tree_results->SetBranchAddress("jet_pt_corr", &jet_pt_corr);
    
    float x_min;
    float x_max;
    int   x_bins;
    char plot_file_name_ext[200];
    char plot_xaxis[50];
    char plot_yaxis[50];
    char plot_type[50];
    
    if (use_normalized) {
        sprintf(plot_yaxis, "N_{Jets}");
    }
    else {
        sprintf(plot_yaxis, "Probability Density");
    }
    
    if (use_delta) {
        x_min = -30.;
        x_max = 30.;
        x_bins = int((x_max - x_min) / 1);
        sprintf(plot_file_name_ext, "%s_delta.pdf", plot_file_name);
        sprintf(plot_xaxis, "p_{T}^{Reco} - p_{T}^{True} [GeV]");
        sprintf(plot_type, "_delta");
    }
    else {
        x_min = test_pt_min - 28.;
        x_max = test_pt_max + 28.;
        x_bins = int((x_max - x_min) / 1);
        sprintf(plot_file_name_ext, "%s.pdf", plot_file_name);
        sprintf(plot_xaxis, "p_{T, ch jet} [GeV]");
        sprintf(plot_type, "");
    }
    
    char plot_labels[200];
    sprintf(plot_labels, "%s; %s; %s", plot_title, plot_xaxis, plot_yaxis);
    TH1F* th1f_canvas_plot = new TH1F("th1f_canvas_plot", plot_labels, x_bins, x_min, x_max);
    th1f_canvas_plot->SetAxisRange(0., 0.2, "Y");
    
    char true_name[100];
    sprintf(true_name, "th1f_jet_pt_true_test_%.0f_%.0f%s", test_pt_min, test_pt_max, plot_type);
    TH1F* th1f_jet_pt_true = new TH1F(true_name, "True Jet pT; p_{T, ch jet} [GeV]; Probability Density", x_bins, x_min, x_max);
    
    char reco_name[100];
    sprintf(reco_name, "th1f_jet_pt_reco_test_%.0f_%.0f%s", test_pt_min, test_pt_max, plot_type);
    TH1F* th1f_jet_pt_reco = new TH1F(reco_name, "Jet pT Reconstructed by ML; p_{T, ch jet} [GeV]; Probability Density", x_bins, x_min, x_max);
    
    char corr_name[100];
    sprintf(corr_name, "th1f_jet_pt_corr_test_%.0f_%.0f%s", test_pt_min, test_pt_max, plot_type);
    TH1F* th1f_jet_pt_corr = new TH1F(corr_name, "Jet pT from Simple Correction; p_{T, ch jet} [GeV]; Probability Density", x_bins, x_min, x_max);
    
    th1f_jet_pt_true->Sumw2();
    th1f_jet_pt_reco->Sumw2();
    th1f_jet_pt_corr->Sumw2();
    
    for ( int j=0 ; j < tree_results->GetEntries() ; j++ ) {
        tree_results->GetEntry(j);
        th1f_jet_pt_true->Fill(jet_pt_true);
        
        if (use_delta) {
            th1f_jet_pt_reco->Fill(jet_pt_reco - jet_pt_true);
            th1f_jet_pt_corr->Fill(jet_pt_corr - jet_pt_true);
        }
        else {
            th1f_jet_pt_reco->Fill(jet_pt_reco);
            th1f_jet_pt_corr->Fill(jet_pt_corr);
        }
        
//        std::cout << jet_pt_true << ", " << jet_pt_reco << ", " << jet_pt_corr << std::endl;
    }
    
    // Uncomment to plot probability density
    if (!use_delta) th1f_jet_pt_true->Scale( 0.5 / th1f_jet_pt_true->Integral(),"WIDTH");
    th1f_jet_pt_reco->Scale( 1. / th1f_jet_pt_reco->Integral(),"WIDTH");
    th1f_jet_pt_corr->Scale( 1. / th1f_jet_pt_corr->Integral(),"WIDTH");
    
    std::cout << "Made it here" << std::endl;
    // TCanvas("name", "title", width (px), height (px))
    TCanvas* canvas = new TCanvas("canvas", "", 1000, 600);

    // Turns on ticks, turns off stats
    gPad->SetTicks();
    gStyle->SetOptStat(0);

    // Turns on Latex formatting relative to canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(kTRUE);
    
    // Set plot colors and styles
    th1f_jet_pt_true->SetLineColor(jet_blk_line);
    th1f_jet_pt_true->SetMarkerColor(jet_blk_mark);
    th1f_jet_pt_true->SetMarkerStyle(mark_circ_open[0]);
    th1f_jet_pt_true->SetMarkerSize(mark_circ_open[1]);
    
    th1f_jet_pt_reco->SetLineColor(jet_blu_line);
    th1f_jet_pt_reco->SetMarkerColor(jet_blu_mark);
    th1f_jet_pt_reco->SetMarkerStyle(mark_star_open[0]);
    th1f_jet_pt_reco->SetMarkerSize(mark_star_open[1]);
    
    th1f_jet_pt_corr->SetLineColor(jet_red_line);
    th1f_jet_pt_corr->SetMarkerColor(jet_red_mark);
    th1f_jet_pt_corr->SetMarkerStyle(mark_squa_open[0]);
    th1f_jet_pt_corr->SetMarkerSize(mark_squa_open[1]);
    
    TLegend* legend = new TLegend(0.15, 0.65, 0.35, 0.85);
    char reco_stats[100];
    float reco_mean   = th1f_jet_pt_reco->GetMean();
    float reco_stddev = th1f_jet_pt_reco->GetStdDev();
    sprintf(reco_stats, "Mean: %.2f, StdDev: %.2f", reco_mean, reco_stddev);
    
    char corr_stats[100];
    float corr_mean   = th1f_jet_pt_corr->GetMean();
    float corr_stddev = th1f_jet_pt_corr->GetStdDev();
    sprintf(corr_stats, "Mean: %.2f, StdDev: %.2f", corr_mean, corr_stddev);
    
    if (!use_delta) legend->AddEntry(th1f_jet_pt_true,  "p_{T}^{True} (0.5 scale)", "lp");
    legend->AddEntry(th1f_jet_pt_reco,  "Linear Regression", "lp");
    legend->AddEntry((TObject*)0, reco_stats, "");
    legend->AddEntry(th1f_jet_pt_corr,  "Area Correction", "lp");
    legend->AddEntry((TObject*)0, corr_stats, "");
    legend->SetLineWidth(0);
    legend->SetFillStyle(0);
    
    th1f_canvas_plot->Draw();
    
    if (!use_delta) th1f_jet_pt_true->Draw("same");
    th1f_jet_pt_corr->Draw("same");
    th1f_jet_pt_reco->Draw("same");
    
    legend->Draw("same");
    
    canvas->Print(plot_file_name_ext);
    std::cout << "Plotted " << plot_file_name_ext << std::endl;
    
    delete legend;
    delete latex;
    delete canvas;
    delete th1f_jet_pt_true;
    delete th1f_jet_pt_corr;
    delete th1f_jet_pt_reco;
    delete input_file;
    std::cout << "Made it to the end" << std::endl;
}



void Jet_ML_Plotter_ROOT() {
    
    char* file_dir = dir_data;
    char output_file_name[100];
    sprintf(output_file_name, "ML_Results/Train_B8_10_90.root");
    
    char ml_weights_array[3][2][100] = {
        {"ML_Results/Train_B8_F1_10_90_LR_Coeffs.csv", "Tree_Train_B8_F1_10_90_LR_Coeffs"},
        {"ML_Results/Train_B8_F3_10_90_LR_Coeffs.csv", "Tree_Train_B8_F3_10_90_LR_Coeffs"},
        {"ML_Results/Train_B8_F12_10_90_LR_Coeffs.csv", "Tree_Train_B8_F12_10_90_LR_Coeffs"}
    };
    
    char ml_results_array[21][2][100] = {
        {"ML_Results/Train_B8_F1_10_90_Test_18_22.csv", "Tree_Train_B8_F1_10_90_Test_18_22"},
        {"ML_Results/Train_B8_F1_10_90_Test_28_32.csv", "Tree_Train_B8_F1_10_90_Test_28_32"},
        {"ML_Results/Train_B8_F1_10_90_Test_38_42.csv", "Tree_Train_B8_F1_10_90_Test_38_42"},
        {"ML_Results/Train_B8_F1_10_90_Test_48_52.csv", "Tree_Train_B8_F1_10_90_Test_48_52"},
        {"ML_Results/Train_B8_F1_10_90_Test_58_62.csv", "Tree_Train_B8_F1_10_90_Test_58_62"},
        {"ML_Results/Train_B8_F1_10_90_Test_68_72.csv", "Tree_Train_B8_F1_10_90_Test_68_72"},
        {"ML_Results/Train_B8_F1_10_90_Test_78_82.csv", "Tree_Train_B8_F1_10_90_Test_78_82"},
        {"ML_Results/Train_B8_F3_10_90_Test_18_22.csv", "Tree_Train_B8_F3_10_90_Test_18_22"},
        {"ML_Results/Train_B8_F3_10_90_Test_28_32.csv", "Tree_Train_B8_F3_10_90_Test_28_32"},
        {"ML_Results/Train_B8_F3_10_90_Test_38_42.csv", "Tree_Train_B8_F3_10_90_Test_38_42"},
        {"ML_Results/Train_B8_F3_10_90_Test_48_52.csv", "Tree_Train_B8_F3_10_90_Test_48_52"},
        {"ML_Results/Train_B8_F3_10_90_Test_58_62.csv", "Tree_Train_B8_F3_10_90_Test_58_62"},
        {"ML_Results/Train_B8_F3_10_90_Test_68_72.csv", "Tree_Train_B8_F3_10_90_Test_68_72"},
        {"ML_Results/Train_B8_F3_10_90_Test_78_82.csv", "Tree_Train_B8_F3_10_90_Test_78_82"},
        {"ML_Results/Train_B8_F12_10_90_Test_18_22.csv", "Tree_Train_B8_F12_10_90_Test_18_22"},
        {"ML_Results/Train_B8_F12_10_90_Test_28_32.csv", "Tree_Train_B8_F12_10_90_Test_28_32"},
        {"ML_Results/Train_B8_F12_10_90_Test_38_42.csv", "Tree_Train_B8_F12_10_90_Test_38_42"},
        {"ML_Results/Train_B8_F12_10_90_Test_48_52.csv", "Tree_Train_B8_F12_10_90_Test_48_52"},
        {"ML_Results/Train_B8_F12_10_90_Test_58_62.csv", "Tree_Train_B8_F12_10_90_Test_58_62"},
        {"ML_Results/Train_B8_F12_10_90_Test_68_72.csv", "Tree_Train_B8_F12_10_90_Test_68_72"},
        {"ML_Results/Train_B8_F12_10_90_Test_78_82.csv", "Tree_Train_B8_F12_10_90_Test_78_82"},
    };
    
    float test_min_max_array[7][2] = {
        {18., 22.}, {28., 32.}, {38., 42.}, {48., 52.}, {58., 62.}, {68., 72.}, {78., 82.}
    };
    
//    // Converts ML estimator feature weights into a TTree
//    for ( int i=0 ; i < 3 ; i++ ) {
//        Build_Weights_TTree_FromCSV(
//            file_dir,
//            ml_weights_array[i][0],
//            output_file_name,
//            ml_weights_array[i][1]
//            );
//    }
//
//    // Converts ML results into a TTree
//    for ( int i=0 ; i < 21 ; i++ ) {
//        Build_Results_TTree_FromCSV(
//            file_dir,
//            ml_results_array[i][0],
//            output_file_name,
//            ml_results_array[i][1]
//            );
//    }

    // Iterates through ML results to output plots
    for ( int i=0 ; i<7 ; i++ ) {

        char input_file_name[300];
        sprintf(
            input_file_name,
            "%s/%s",
            file_dir,
            output_file_name);

        char plot_file_name_actual[300];
        sprintf(
            plot_file_name_actual,
            "%s/ML_Results/Plots_Actual/Plot_Train_B8_F12_10_90_Test_%.0f_%.0f",
            file_dir,
            test_min_max_array[i][0],
            test_min_max_array[i][1]);
            
        char plot_file_name_delta[300];
        sprintf(
            plot_file_name_delta,
            "%s/ML_Results/Plots_Delta/Plot_Train_B8_F12_10_90_Test_%.0f_%.0f",
            file_dir,
            test_min_max_array[i][0],
            test_min_max_array[i][1]);

        char plot_title[100];
        sprintf(
            plot_title,
            "Test: p_{T}^{True} in %.0f-%.0f GeV [Train: p_{T}^{True} in 10-90 GeV, p_{T}^{8} Bias, 12 Features]",
            test_min_max_array[i][0],
            test_min_max_array[i][1]);

        char ml_results[100];
        sprintf(ml_results, "%s", ml_results_array[14+i][1]);
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

        std::cout << "Plot Completed" << std::endl;
    }
    
}

