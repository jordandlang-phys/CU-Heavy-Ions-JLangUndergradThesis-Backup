#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include <cmath>
#include <iostream>
#include <filesystem>
#include <algorithm>

#include "_header.h"
using namespace Project_Constants;

const bool printOut = true;
const bool debug = true;

TFile* input_file;
const int max_jets = 100;
const int max_parts = 400;

double Gaussian_Func(double *x,double *par) {
    double arg = 0;
    if (par[2]!=0) arg = (x[0] - par[1])/par[2];
    double fitval = par[0]*TMath::Exp(-0.5*arg*arg);
    return fitval;
}

TF1* tf1_paper_simple_correction_fit;
TF1* tf1_paper_neural_network_fit;
TF1* tf1_paper_random_forest_fit;
TF1* tf1_paper_linear_regression_fit;

void Fit_Paper_Data() {
    // Create Graphs to Fit Paper Data
    double tgraph_paper_sc_x[25] = {-21., -19., -17., -15., -13., -11., -9., -7., -5., -3., -1., 1., 3., 5., 7., 9., 11., 12., 15., 17., 19., 21., 23., 25., 27.};
    double tgraph_paper_sc_y[25] = {-0.0003636363636363440, 0.0007272727272727430, 0.0018181818181818300, 0.0040000000000000000, 0.009090909090909090, 0.015272727272727300, 0.024363636363636400, 0.034181818181818200, 0.04000000000000000, 0.048727272727272700, 0.049090909090909100, 0.05127272727272730, 0.048000000000000000, 0.04436363636363640, 0.03672727272727280, 0.029454545454545500, 0.02145454545454550, 0.015636363636363600, 0.010545454545454600, 0.00545454545454549, 0.0032727272727272900, 0.0029090909090909200, 0.00036363636363637200, 0.00036363636363637200, 0.0};
    
    TGraph* tgraph_paper_simple_correction = new TGraph(25, tgraph_paper_sc_x, tgraph_paper_sc_y);
    tf1_paper_simple_correction_fit = new TF1("tf1_paper_simple_correction_fit", Gaussian_Func, -40., 40., 3);
    tf1_paper_simple_correction_fit->SetParameters(.1, 0., 5.);
    tgraph_paper_simple_correction->Fit(tf1_paper_simple_correction_fit,"RNL");
    
    double tgraph_paper_nn_x[16] = {-15., -13., -11., -9., -7., -5., -3., -1., 1., 3., 5., 7., 9., 11., 12., 15.};
    double tgraph_paper_nn_y[16] = {0.0, 0.0, 0.0007272727272727430, 0.003636363636363630, 0.015272727272727300, 0.03563636363636370, 0.07272727272727270, 0.10581818181818200, 0.11054545454545500, 0.07854545454545460, 0.04290909090909090, 0.02036363636363640, 0.00872727272727275, 0.0029090909090909200, 0.00036363636363637200, 0.0};
    
    TGraph* tgraph_paper_neural_network = new TGraph(16, tgraph_paper_nn_x, tgraph_paper_nn_y);
    tf1_paper_neural_network_fit = new TF1("tf1_paper_neural_network_fit", Gaussian_Func, -40., 40., 3);
    tf1_paper_neural_network_fit->SetParameters(.1, 0., 5.);
    tgraph_paper_neural_network->Fit(tf1_paper_neural_network_fit,"RNL");
    
    double tgraph_paper_rf_x[16] = {-15., -13., -11., -9., -7., -5., -3., -1., 1., 3., 5., 7., 9., 11., 12., 15.};
    double tgraph_paper_rf_y[16] = {0.0, 0.0, 0.0014545454545454600, 0.0058181818181818300, 0.02036363636363640, 0.04618181818181820, 0.07927272727272730, 0.09854545454545460, 0.09600000000000000, 0.0701818181818182, 0.042181818181818200, 0.023272727272727300, 0.009454545454545470, 0.003636363636363630, 0.0007272727272727430, 0.0};
    
    TGraph* tgraph_paper_random_forest = new TGraph(16, tgraph_paper_rf_x, tgraph_paper_rf_y);
    tf1_paper_random_forest_fit = new TF1("tf1_paper_random_forest_fit", Gaussian_Func, -40., 40., 3);
    tf1_paper_random_forest_fit->SetParameters(.1, 0., 5.);
    tgraph_paper_random_forest->Fit(tf1_paper_random_forest_fit,"RNL");
    
    double tgraph_paper_lr_x[16] = {-15., -13., -11., -9., -7., -5., -3., -1., 1., 3., 5., 7., 9., 11., 12., 15.};
    double tgraph_paper_lr_y[16] = {0.0, 0.00036363636363637200, 0.002545454545454570, 0.009818181818181840, 0.02872727272727270, 0.052727272727272700, 0.08145454545454550, 0.09454545454545460, 0.08836363636363640, 0.06618181818181820, 0.03927272727272730, 0.021090909090909100, 0.009818181818181840, 0.003636363636363630, 0.0007272727272727430, 0.0};
    
    TGraph* tgraph_paper_linear_regression = new TGraph(16, tgraph_paper_lr_x, tgraph_paper_lr_y);
    tf1_paper_linear_regression_fit = new TF1("tf1_paper_linear_regression_fit", Gaussian_Func, -40., 40., 3);
    tf1_paper_linear_regression_fit->SetParameters(.1, 0., 5.);
    tgraph_paper_linear_regression->Fit(tf1_paper_linear_regression_fit,"RNL");
}

void Plot_ML_pT_Comparison(char* plot_filename, char* plot_directory, char* plot_title_xlabel_ylabel, TH1D* th1d_simple_correction, TH1D* th1d_linear_regression, TH1D* th1d_random_forest, TH1D* th1d_neural_network, TH1D* th1d_feature_list, bool show_feature_list, bool show_paper_plots) {
    
    double x_min  = th1d_simple_correction->GetXaxis()->GetXmin();
    double x_max  = th1d_simple_correction->GetXaxis()->GetXmax();
    double x_bins = th1d_simple_correction->GetNbinsX();
    double max_arr[4] = {th1d_simple_correction->GetMaximum(), th1d_neural_network->GetMaximum(), th1d_random_forest->GetMaximum(), th1d_linear_regression->GetMaximum()};
    double y_max = 0.0;
    for ( int i = 0 ; i < 4 ; i++ ) {
        if ( y_max < max_arr[i] ) y_max = max_arr[i];
    }
    if (show_paper_plots) y_max = 0.20;
    
    TH1D* th1d_underlying_plot = new TH1D("th1d_underlying_plot", plot_title_xlabel_ylabel, x_bins, x_min, x_max);
    th1d_underlying_plot->SetAxisRange(0, 1.2 * y_max, "Y");
    
    // Fit to Each Plot
    TF1* tf1_simple_correction_fit  = new TF1("tf1_simple_correction_fit",  Gaussian_Func, x_min, x_max, 3);
    TF1* tf1_linear_regression_fit  = new TF1("tf1_linear_regression_fit",  Gaussian_Func, x_min, x_max, 3);
    TF1* tf1_random_forest_fit      = new TF1("tf1_random_forest_fit",      Gaussian_Func, x_min, x_max, 3);
    TF1* tf1_neural_network_fit     = new TF1("tf1_neural_network_fit",     Gaussian_Func, x_min, x_max, 3);
    if (show_paper_plots) {
        tf1_simple_correction_fit   ->SetParameters(0.8, 0., 5);
        tf1_linear_regression_fit   ->SetParameters(0.8, 0., 5);
        tf1_random_forest_fit       ->SetParameters(0.8, 0., 5);
        tf1_neural_network_fit      ->SetParameters(0.8, 0., 5);
    }
    else {
        tf1_simple_correction_fit   ->SetParameters(y_max/2., 0., 0.1);
        tf1_linear_regression_fit   ->SetParameters(y_max/2., 0., 0.1);
        tf1_random_forest_fit       ->SetParameters(y_max/2., 0., 0.1);
        tf1_neural_network_fit      ->SetParameters(y_max/2., 0., 0.1);
    }
    th1d_simple_correction  ->Fit(tf1_simple_correction_fit,"RNL");
    th1d_linear_regression  ->Fit(tf1_linear_regression_fit,"RNL");
    th1d_random_forest      ->Fit(tf1_random_forest_fit,"RNL");
    th1d_neural_network     ->Fit(tf1_neural_network_fit,"RNL");
    
    // TCanvas("name", "title", width (px), height (px))
    TCanvas* canvas = new TCanvas("canvas", "", 1000, 600);

    // Turns on ticks, turns off stats
    gPad->SetTicks();
    gStyle->SetOptStat(0);

    // Turns on Latex formatting relative to canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(kTRUE);
    
    // Build histogram plot
    th1d_simple_correction->SetLineColor(jet_red_line);
    th1d_simple_correction->SetMarkerColor(jet_red_mark);
    th1d_simple_correction->SetMarkerStyle(mark_circ_open[0]);
    th1d_simple_correction->SetMarkerSize(mark_circ_open[1]);
    tf1_simple_correction_fit->SetLineColor(jet_red_line);
    tf1_simple_correction_fit->SetLineWidth(1);
    
    th1d_linear_regression->SetLineColor(jet_vio_line);
    th1d_linear_regression->SetMarkerColor(jet_vio_mark);
    th1d_linear_regression->SetMarkerStyle(mark_squa_open[0]);
    th1d_linear_regression->SetMarkerSize(mark_squa_open[1]);
    tf1_linear_regression_fit->SetLineColor(jet_vio_line);
    tf1_linear_regression_fit->SetLineWidth(1);
    
    th1d_random_forest->SetLineColor(jet_tea_line);
    th1d_random_forest->SetMarkerColor(jet_tea_mark);
    th1d_random_forest->SetMarkerStyle(mark_diam_open[0]);
    th1d_random_forest->SetMarkerSize(mark_diam_open[1]);
    tf1_random_forest_fit->SetLineColor(jet_tea_line);
    tf1_random_forest_fit->SetLineWidth(1);

    th1d_neural_network->SetLineColor(jet_blu_line);
    th1d_neural_network->SetMarkerColor(jet_blu_mark);
    th1d_neural_network->SetMarkerStyle(mark_star_open[0]);
    th1d_neural_network->SetMarkerSize(mark_star_open[1]);
    tf1_neural_network_fit->SetLineColor(jet_blu_line);
    tf1_neural_network_fit->SetLineWidth(1);
    
    // Draws the histograms and lines
    th1d_underlying_plot        ->Draw();
    if (show_paper_plots) {
        Fit_Paper_Data();
        
        tf1_paper_neural_network_fit->SetLineColor(jet_blu_line);
        tf1_paper_neural_network_fit->SetLineWidth(1);
        tf1_paper_neural_network_fit->SetLineStyle(3);
    
        tf1_paper_random_forest_fit->SetLineColor(jet_tea_line);
        tf1_paper_random_forest_fit->SetLineWidth(1);
        tf1_paper_random_forest_fit->SetLineStyle(3);
    
        tf1_paper_linear_regression_fit->SetLineColor(jet_vio_line);
        tf1_paper_linear_regression_fit->SetLineWidth(1);
        tf1_paper_linear_regression_fit->SetLineStyle(3);
    
        tf1_paper_simple_correction_fit->SetLineColor(jet_red_line);
        tf1_paper_simple_correction_fit->SetLineWidth(1);
        tf1_paper_simple_correction_fit->SetLineStyle(3);
        
        tf1_paper_simple_correction_fit ->Draw("same");
        tf1_paper_neural_network_fit    ->Draw("same");
        tf1_paper_random_forest_fit     ->Draw("same");
        tf1_paper_linear_regression_fit ->Draw("same");
    }
    tf1_simple_correction_fit   ->Draw("same");
    tf1_neural_network_fit      ->Draw("same");
    tf1_random_forest_fit       ->Draw("same");
    tf1_linear_regression_fit   ->Draw("same");
    th1d_simple_correction      ->Draw("same");
    th1d_neural_network         ->Draw("same");
    th1d_random_forest          ->Draw("same");
    th1d_linear_regression      ->Draw("same");

    // Draws Legend as TLegend(xmin, ymin, xmax, ymax)
    TLegend* legend;
    if (show_paper_plots) legend = new TLegend(0.15, 0.45, 0.40, 0.85);
    else legend = new TLegend(0.15, 0.65, 0.40, 0.85);
    legend->AddEntry(th1d_random_forest,        "Random Forest", "lp");
    legend->AddEntry(th1d_neural_network,       "Neural Network", "lp");
    legend->AddEntry(th1d_linear_regression,    "Linear Regression", "lp");
    legend->AddEntry(th1d_simple_correction,    "Simple Correction", "lp");
    if (show_paper_plots) {
        legend->AddEntry(tf1_paper_random_forest_fit,       "RF (Paper)", "l");
        legend->AddEntry(tf1_paper_neural_network_fit,      "NN (Paper)", "l");
        legend->AddEntry(tf1_paper_linear_regression_fit,   "LR (Paper)", "l");
        legend->AddEntry(tf1_paper_simple_correction_fit,   "SC (Paper)", "l");
    }
    legend->SetLineWidth(0);
    legend->SetFillStyle(0);
    legend->Draw();
    
    // Draws Feature List
    int  feature_count = th1d_feature_list->GetNbinsX();
    char feature_list_8[9][100];
    char feature_vals_8[9][100];
    char feature_list_10[11][100];
    char feature_vals_10[11][100];
    if (show_feature_list && feature_count == 8) {
        sprintf(feature_list_8[0], "#scale[0.6]{Feature Importance}");
        sprintf(feature_list_8[1], "#scale[0.6]{#bf{Jet p_{T, raw}:}}");
        sprintf(feature_list_8[2], "#scale[0.6]{#bf{Jet p_{T, corr}:}}");
        sprintf(feature_list_8[3], "#scale[0.6]{#bf{Mean Const. p_{T}:}}");
        sprintf(feature_list_8[4], "#scale[0.6]{#bf{Median Const. p_{T}:}}");
        sprintf(feature_list_8[5], "#scale[0.6]{#bf{p_{T, Const.}^{1}:}}");
        sprintf(feature_list_8[6], "#scale[0.6]{#bf{p_{T, Const.}^{2}:}}");
        sprintf(feature_list_8[7], "#scale[0.6]{#bf{p_{T, Const.}^{3}:}}");
        sprintf(feature_list_8[8], "#scale[0.6]{#bf{p_{T, Const.}^{4}:}}");
        
        sprintf(feature_vals_8[0], "#scale[0.6]{}");
        sprintf(feature_vals_8[1], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(1));
        sprintf(feature_vals_8[2], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(2));
        sprintf(feature_vals_8[3], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(3));
        sprintf(feature_vals_8[4], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(4));
        sprintf(feature_vals_8[5], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(5));
        sprintf(feature_vals_8[6], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(6));
        sprintf(feature_vals_8[7], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(7));
        sprintf(feature_vals_8[8], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(8));
        
        for ( int i = 0 ; i <= 8 ; i++ ) {
            latex->DrawLatex(0.63, (.813 - 0.0494 * i), feature_list_8[i]);
        }
        for ( int i = 0 ; i <= 8 ; i++ ) {
            latex->DrawLatex(0.80, (.813 - 0.0494 * i), feature_vals_8[i]);
        }
    }
    if (show_feature_list && feature_count == 10) {
        sprintf(feature_list_10[0], "#scale[0.6]{Feature Importance}");
        sprintf(feature_list_10[1], "#scale[0.6]{#bf{Jet p_{T, corr}:}}");
        sprintf(feature_list_10[2], "#scale[0.6]{#bf{Jet Mass:}}");
        sprintf(feature_list_10[3], "#scale[0.6]{#bf{Jet Area:}}");
        sprintf(feature_list_10[4], "#scale[0.6]{#bf{N_{const}:}}");
        sprintf(feature_list_10[5], "#scale[0.6]{#bf{Mean Const. p_{T}:}}");
        sprintf(feature_list_10[6], "#scale[0.6]{#bf{p_{T, Const.}^{1}:}}");
        sprintf(feature_list_10[7], "#scale[0.6]{#bf{p_{T, Const.}^{2}:}}");
        sprintf(feature_list_10[8], "#scale[0.6]{#bf{p_{T, Const.}^{3}:}}");
        sprintf(feature_list_10[9], "#scale[0.6]{#bf{Jet y:}}");
        sprintf(feature_list_10[10], "#scale[0.6]{#bf{Jet #rho:}}");
        
        sprintf(feature_vals_10[0], "#scale[0.6]{}");
        sprintf(feature_vals_10[1], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(1));
        sprintf(feature_vals_10[2], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(2));
        sprintf(feature_vals_10[3], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(3));
        sprintf(feature_vals_10[4], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(4));
        sprintf(feature_vals_10[5], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(5));
        sprintf(feature_vals_10[6], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(6));
        sprintf(feature_vals_10[7], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(7));
        sprintf(feature_vals_10[8], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(8));
        sprintf(feature_vals_10[9], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(7));
        sprintf(feature_vals_10[10], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(8));
        
        for ( int i = 0 ; i <= 10 ; i++ ) {
            latex->DrawLatex(0.63, (.813 - 0.0494 * i), feature_list_10[i]);
        }
        for ( int i = 0 ; i <= 10 ; i++ ) {
            latex->DrawLatex(0.80, (.813 - 0.0494 * i), feature_vals_10[i]);
        }
    }
    
    // Prints out the plot
    char plot_output[300];
    sprintf(plot_output, "%s/%s", plot_directory, plot_filename);
    canvas->Print(plot_output);
    std::cout << "Plotted " << plot_filename << std::endl;
    gPad->SetLogy(0);
    
    return 0;
}

// ##############################
// #                            #
// #    MACRO STARTS HERE!!!    #
// #                            #
// ##############################

const char* input_file_name = "Jet_ML_Prep.root";
char* pt_train_bin = "40_60_pt3_Bias";
//char* pt_test_bin_arr[11] = {
//    "30_35", "35_40", "40_45", "45_50", "50_55", "55_60", "60_65", "65_70",
//    "40_60", "40_50", "40_41"};
//char* pt_test_bin_title_arr[11] = {
//    "30-35 GeV", "35-40 GeV", "40-45 GeV", "45-50 GeV", "50-55 GeV", "55-60 GeV", "60-65 GeV", "65-70 GeV",
//    "40-60 GeV", "40-50 GeV", "40-41 GeV"};
char* pt_test_bin_arr[4] = {
    "40_45", "45_50", "50_55", "55_60"};
char* pt_test_bin_title_arr[4] = {
    "40-45 GeV", "45-50 GeV", "50-55 GeV", "55-60 GeV"};
char* th1d_name_arr[40] = {
    "8feature_ptTrueA_feature_importance",
    "8feature_ptTrueA_simple_correction",   "8feature_ptTrueA_linear_regression",
    "8feature_ptTrueA_random_forest",       "8feature_ptTrueA_neural_network",
    "8feature_ptTrueA_compare_feature_importance",
    "8feature_ptTrueA_compare_simple_correction",   "8feature_ptTrueA_compare_linear_regression",
    "8feature_ptTrueA_compare_random_forest",       "8feature_ptTrueA_compare_neural_network",
    "8feature_ptTrueB_feature_importance",
    "8feature_ptTrueB_simple_correction",   "8feature_ptTrueB_linear_regression",
    "8feature_ptTrueB_random_forest",       "8feature_ptTrueB_neural_network",
    "8feature_ptTrueB_compare_feature_importance",
    "8feature_ptTrueB_compare_simple_correction",   "8feature_ptTrueB_compare_linear_regression",
    "8feature_ptTrueB_compare_random_forest",       "8feature_ptTrueB_compare_neural_network",
    "only_ptRaw_ptTrueA_compare_feature_importance",
    "only_ptRaw_ptTrueA_compare_simple_correction", "only_ptRaw_ptTrueA_compare_linear_regression",
    "only_ptRaw_ptTrueA_compare_random_forest",     "only_ptRaw_ptTrueA_compare_neural_network",
    "only_ptCorr_ptTrueA_compare_feature_importance",
    "only_ptCorr_ptTrueA_compare_simple_correction", "only_ptCorr_ptTrueA_compare_linear_regression",
    "only_ptCorr_ptTrueA_compare_random_forest",     "only_ptCorr_ptTrueA_compare_neural_network",
    "10feature_ptTrueA_compare_feature_importance",
    "10feature_ptTrueA_compare_simple_correction",  "10feature_ptTrueA_compare_linear_regression",
    "10feature_ptTrueA_compare_random_forest",      "10feature_ptTrueA_compare_neural_network",
    "10feature_ptTrueB_compare_feature_importance",
    "10feature_ptTrueB_compare_simple_correction",  "10feature_ptTrueB_compare_linear_regression",
    "10feature_ptTrueB_compare_random_forest",      "10feature_ptTrueB_compare_neural_network"
};

void Jet_Plotter_Function_ROOT() {
    
    // Opens and reads the Root output file
    char input_file_path[200];
    sprintf(input_file_path, "%s/%s", dir_data, input_file_name);
    input_file = new TFile(input_file_path, "READ");
    std::cout << "Input file read." << std::endl;
    
    char subdir_plots [200];
    sprintf(subdir_plots, "%s/MachineLearning", dir_plots);
    std::__fs::filesystem::create_directories(subdir_plots);
    
    // MACHINE LEARNING PLOTS
    
    
    // EXAMINING TOP 8 FEATURES
    for ( int pt_test_bin_index = 0 ; pt_test_bin_index < 11 ; pt_test_bin_index++ ) {
        char  plot_file_name[400];
        char  plot_labels[400];
        char  th1d_feature_importance[200];
        TH1D* th1d_feature_importance_list;
        char  th1d_simple_correction[200];
        TH1D* th1d_simple_correction_plot;
        char  th1d_linear_regression[200];
        TH1D* th1d_linear_regression_plot;
        char  th1d_random_forest[200];
        TH1D* th1d_random_forest_plot;
        char  th1d_neural_network[200];
        TH1D* th1d_neural_network_plot;
        
        
//        // Plot 8 input features, pt_true method A, no comparison
//        sprintf(plot_file_name, "th1d_%s_jet_pt_delta_8feature_ptTrueA.pdf", pt_test_bin_arr[pt_test_bin_index]);
//        std::cout << "----- Generating: " << plot_file_name << " -----" << std::endl;
//
//        sprintf(plot_labels, "Jet p_{T} Delta for %s #left[p_{T}^{true} = p_{T, ch jet}^{PYTHIA} #right]; (p_{T, ch jet}^{reco} - p_{T, ch jet}^{true})/p_{T, ch jet}^{true}; N_{ch jets}", pt_test_bin_title_arr[pt_test_bin_index]);
//
//        sprintf(th1d_feature_importance, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[0]);
//        th1d_feature_importance_list    = (TH1D*) input_file->Get(th1d_feature_importance);
//
//        sprintf(th1d_simple_correction, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[1]);
//        th1d_simple_correction_plot     = (TH1D*) input_file->Get(th1d_simple_correction);
//
//        sprintf(th1d_linear_regression, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[2]);
//        th1d_linear_regression_plot     = (TH1D*) input_file->Get(th1d_linear_regression);
//
//        sprintf(th1d_random_forest, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[3]);
//        th1d_random_forest_plot         = (TH1D*) input_file->Get(th1d_random_forest);
//
//        sprintf(th1d_neural_network, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[4]);
//        th1d_neural_network_plot        = (TH1D*) input_file->Get(th1d_neural_network);
//
//        Plot_ML_pT_Comparison(
//            plot_file_name,
//            subdir_plots,
//            plot_labels,
//            th1d_simple_correction_plot,
//            th1d_linear_regression_plot,
//            th1d_random_forest_plot,
//            th1d_neural_network_plot,
//            th1d_feature_importance_list,
//            true, // show feature list
//            false); // compare to paper plots
//
//
//        // Plot 8 input features, pt_true method B, no comparison
//        sprintf(plot_file_name, "th1d_%s_jet_pt_delta_8feature_ptTrueB.pdf", pt_test_bin_arr[pt_test_bin_index]);
//        std::cout << "----- Generating: " << plot_file_name << " -----" << std::endl;
//
//        sprintf(plot_labels, "Jet p_{T} Delta for %s #left[p_{T}^{true} = p_{T}^{raw} #Sigma p_{T, const. i}^{PYTHIA} / #Sigma p_{T, const. i} #right]; p_{T, ch jet}^{reco} - p_{T, ch jet}^{true} [GeV]; Probability Density", pt_test_bin_title_arr[pt_test_bin_index]);
//
//        sprintf(th1d_feature_importance, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[10]);
//        th1d_feature_importance_list    = (TH1D*) input_file->Get(th1d_feature_importance);
//
//        sprintf(th1d_simple_correction, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[11]);
//        th1d_simple_correction_plot     = (TH1D*) input_file->Get(th1d_simple_correction);
//
//        sprintf(th1d_linear_regression, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[12]);
//        th1d_linear_regression_plot     = (TH1D*) input_file->Get(th1d_linear_regression);
//
//        sprintf(th1d_random_forest, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[13]);
//        th1d_random_forest_plot         = (TH1D*) input_file->Get(th1d_random_forest);
//
//        sprintf(th1d_neural_network, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[14]);
//        th1d_neural_network_plot        = (TH1D*) input_file->Get(th1d_neural_network);
//
//        Plot_ML_pT_Comparison(
//            plot_file_name,
//            subdir_plots,
//            plot_labels,
//            th1d_simple_correction_plot,
//            th1d_linear_regression_plot,
//            th1d_random_forest_plot,
//            th1d_neural_network_plot,
//            th1d_feature_importance_list,
//            true, // show feature list
//            false); // compare to paper plots
        
        
        // Plot 8 input features, pt_true method A, with comparison to paper
        sprintf(plot_file_name, "th1d_%s_jet_pt_delta_8feature_ptTrueA_compare.pdf", pt_test_bin_arr[pt_test_bin_index]);
        std::cout << "----- Generating: " << plot_file_name << " -----" << std::endl;
        
        sprintf(plot_labels, "Jet p_{T} Delta for %s #left[p_{T}^{true} = p_{T, ch jet}^{PYTHIA} #right]; p_{T, ch jet}^{reco} - p_{T, ch jet}^{true} [GeV]; Probability Density", pt_test_bin_title_arr[pt_test_bin_index]);
        
        sprintf(th1d_feature_importance, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[5]);
        th1d_feature_importance_list    = (TH1D*) input_file->Get(th1d_feature_importance);
        
        sprintf(th1d_simple_correction, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[6]);
        th1d_simple_correction_plot     = (TH1D*) input_file->Get(th1d_simple_correction);
        
        sprintf(th1d_linear_regression, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[7]);
        th1d_linear_regression_plot     = (TH1D*) input_file->Get(th1d_linear_regression);
        
        sprintf(th1d_random_forest, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[8]);
        th1d_random_forest_plot         = (TH1D*) input_file->Get(th1d_random_forest);
        
        sprintf(th1d_neural_network, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[9]);
        th1d_neural_network_plot        = (TH1D*) input_file->Get(th1d_neural_network);
        
        th1d_simple_correction_plot     ->Scale( 1. / th1d_simple_correction_plot->Integral(),"WIDTH");
        th1d_linear_regression_plot     ->Scale( 1. / th1d_linear_regression_plot->Integral(),"WIDTH");
        th1d_random_forest_plot         ->Scale( 1. / th1d_random_forest_plot->Integral(),"WIDTH");
        th1d_neural_network_plot        ->Scale( 1. / th1d_neural_network_plot->Integral(),"WIDTH");
        
        Plot_ML_pT_Comparison(
            plot_file_name,
            subdir_plots,
            plot_labels,
            th1d_simple_correction_plot,
            th1d_linear_regression_plot,
            th1d_random_forest_plot,
            th1d_neural_network_plot,
            th1d_feature_importance_list,
            true, // show feature list
            true); // compare to paper plots
        
        
//        // Plot 8 input features, pt_true method B, with comparison to paper
//        sprintf(plot_file_name, "th1d_%s_jet_pt_delta_8feature_ptTrueB_compare.pdf", pt_test_bin_arr[pt_test_bin_index]);
//        std::cout << "----- Generating: " << plot_file_name << " -----" << std::endl;
//
//        sprintf(plot_labels, "Jet p_{T} Delta for %s #left[p_{T}^{true} = p_{T}^{raw} #Sigma p_{T, const. i}^{PYTHIA} / #Sigma p_{T, const. i} #right]; p_{T, ch jet}^{reco} - p_{T, ch jet}^{true} [GeV]; Probability Density", pt_test_bin_title_arr[pt_test_bin_index]);
//
//        sprintf(th1d_feature_importance, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[15]);
//        th1d_feature_importance_list    = (TH1D*) input_file->Get(th1d_feature_importance);
//
//        sprintf(th1d_simple_correction, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[16]);
//        th1d_simple_correction_plot     = (TH1D*) input_file->Get(th1d_simple_correction);
//
//        sprintf(th1d_linear_regression, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[17]);
//        th1d_linear_regression_plot     = (TH1D*) input_file->Get(th1d_linear_regression);
//
//        sprintf(th1d_random_forest, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[18]);
//        th1d_random_forest_plot         = (TH1D*) input_file->Get(th1d_random_forest);
//
//        sprintf(th1d_neural_network, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[19]);
//        th1d_neural_network_plot        = (TH1D*) input_file->Get(th1d_neural_network);
//
//        th1d_simple_correction_plot     ->Scale( 1. / th1d_simple_correction_plot->Integral(),"WIDTH");
//        th1d_linear_regression_plot     ->Scale( 1. / th1d_linear_regression_plot->Integral(),"WIDTH");
//        th1d_random_forest_plot         ->Scale( 1. / th1d_random_forest_plot->Integral(),"WIDTH");
//        th1d_neural_network_plot        ->Scale( 1. / th1d_neural_network_plot->Integral(),"WIDTH");
//
//        Plot_ML_pT_Comparison(
//            plot_file_name,
//            subdir_plots,
//            plot_labels,
//            th1d_simple_correction_plot,
//            th1d_linear_regression_plot,
//            th1d_random_forest_plot,
//            th1d_neural_network_plot,
//            th1d_feature_importance_list,
//            true, // show feature list
//            true); // compare to paper plots
        
        
//        // Plot only ptRaw, pt_true method A, with comparison to paper
//        sprintf(plot_file_name, "th1d_%s_jet_pt_delta_ptRaw_ptTrueA_compare.pdf", pt_test_bin_arr[pt_test_bin_index]);
//        std::cout << "----- Generating: " << plot_file_name << " -----" << std::endl;
//
//        sprintf(plot_labels, "Jet p_{T} Delta using p_{T}^{raw} for %s #left[p_{T}^{true} = p_{T, ch jet}^{PYTHIA} #right]; p_{T, ch jet}^{reco} - p_{T, ch jet}^{true} [GeV]; Probability Density", pt_test_bin_title_arr[pt_test_bin_index]);
//
//        sprintf(th1d_feature_importance, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[20]);
//        th1d_feature_importance_list    = (TH1D*) input_file->Get(th1d_feature_importance);
//
//        sprintf(th1d_simple_correction, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[21]);
//        th1d_simple_correction_plot     = (TH1D*) input_file->Get(th1d_simple_correction);
//
//        sprintf(th1d_linear_regression, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[22]);
//        th1d_linear_regression_plot     = (TH1D*) input_file->Get(th1d_linear_regression);
//
//        sprintf(th1d_random_forest, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[23]);
//        th1d_random_forest_plot         = (TH1D*) input_file->Get(th1d_random_forest);
//
//        sprintf(th1d_neural_network, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[24]);
//        th1d_neural_network_plot        = (TH1D*) input_file->Get(th1d_neural_network);
//
//        th1d_simple_correction_plot     ->Scale( 1. / th1d_simple_correction_plot->Integral(),"WIDTH");
//        th1d_linear_regression_plot     ->Scale( 1. / th1d_linear_regression_plot->Integral(),"WIDTH");
//        th1d_random_forest_plot         ->Scale( 1. / th1d_random_forest_plot->Integral(),"WIDTH");
//        th1d_neural_network_plot        ->Scale( 1. / th1d_neural_network_plot->Integral(),"WIDTH");
//
//        Plot_ML_pT_Comparison(
//            plot_file_name,
//            subdir_plots,
//            plot_labels,
//            th1d_simple_correction_plot,
//            th1d_linear_regression_plot,
//            th1d_random_forest_plot,
//            th1d_neural_network_plot,
//            th1d_feature_importance_list,
//            false, // show feature list
//            true); // compare to paper plots
//
//        // Plot only ptRaw, pt_true method A, with comparison to paper
//        sprintf(plot_file_name, "th1d_%s_jet_pt_delta_ptCorr_ptTrueA_compare.pdf", pt_test_bin_arr[pt_test_bin_index]);
//        std::cout << "----- Generating: " << plot_file_name << " -----" << std::endl;
//
//        sprintf(plot_labels, "Jet p_{T} Delta using p_{T}^{corr} for %s #left[p_{T}^{true} = p_{T, ch jet}^{PYTHIA} #right]; p_{T, ch jet}^{reco} - p_{T, ch jet}^{true} [GeV]; Probability Density", pt_test_bin_title_arr[pt_test_bin_index]);
//
//        sprintf(th1d_feature_importance, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[25]);
//        th1d_feature_importance_list    = (TH1D*) input_file->Get(th1d_feature_importance);
//
//        sprintf(th1d_simple_correction, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[26]);
//        th1d_simple_correction_plot     = (TH1D*) input_file->Get(th1d_simple_correction);
//
//        sprintf(th1d_linear_regression, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[27]);
//        th1d_linear_regression_plot     = (TH1D*) input_file->Get(th1d_linear_regression);
//
//        sprintf(th1d_random_forest, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[28]);
//        th1d_random_forest_plot         = (TH1D*) input_file->Get(th1d_random_forest);
//
//        sprintf(th1d_neural_network, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[29]);
//        th1d_neural_network_plot        = (TH1D*) input_file->Get(th1d_neural_network);
//
//        th1d_simple_correction_plot     ->Scale( 1. / th1d_simple_correction_plot->Integral(),"WIDTH");
//        th1d_linear_regression_plot     ->Scale( 1. / th1d_linear_regression_plot->Integral(),"WIDTH");
//        th1d_random_forest_plot         ->Scale( 1. / th1d_random_forest_plot->Integral(),"WIDTH");
//        th1d_neural_network_plot        ->Scale( 1. / th1d_neural_network_plot->Integral(),"WIDTH");
//
//        Plot_ML_pT_Comparison(
//            plot_file_name,
//            subdir_plots,
//            plot_labels,
//            th1d_simple_correction_plot,
//            th1d_linear_regression_plot,
//            th1d_random_forest_plot,
//            th1d_neural_network_plot,
//            th1d_feature_importance_list,
//            false, // show feature list
//            true); // compare to paper plots
        
        
//        // Plot 10 input features, pt_true method A, with comparison to paper
//        sprintf(plot_file_name, "th1d_%s_jet_pt_delta_10feature_ptTrueA_compare.pdf", pt_test_bin_arr[pt_test_bin_index]);
//        std::cout << "----- Generating: " << plot_file_name << " -----" << std::endl;
//
//        sprintf(plot_labels, "Jet p_{T} Delta for %s #left[p_{T}^{true} = p_{T, ch jet}^{PYTHIA} #right]; p_{T, ch jet}^{reco} - p_{T, ch jet}^{true} [GeV]; Probability Density", pt_test_bin_title_arr[pt_test_bin_index]);
//
//        sprintf(th1d_feature_importance, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[30]);
//        th1d_feature_importance_list    = (TH1D*) input_file->Get(th1d_feature_importance);
//
//        sprintf(th1d_simple_correction, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[31]);
//        th1d_simple_correction_plot     = (TH1D*) input_file->Get(th1d_simple_correction);
//
//        sprintf(th1d_linear_regression, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[32]);
//        th1d_linear_regression_plot     = (TH1D*) input_file->Get(th1d_linear_regression);
//
//        sprintf(th1d_random_forest, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[33]);
//        th1d_random_forest_plot         = (TH1D*) input_file->Get(th1d_random_forest);
//
//        sprintf(th1d_neural_network, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[34]);
//        th1d_neural_network_plot        = (TH1D*) input_file->Get(th1d_neural_network);
//
//        th1d_simple_correction_plot     ->Scale( 1. / th1d_simple_correction_plot->Integral(),"WIDTH");
//        th1d_linear_regression_plot     ->Scale( 1. / th1d_linear_regression_plot->Integral(),"WIDTH");
//        th1d_random_forest_plot         ->Scale( 1. / th1d_random_forest_plot->Integral(),"WIDTH");
//        th1d_neural_network_plot        ->Scale( 1. / th1d_neural_network_plot->Integral(),"WIDTH");
//
//        Plot_ML_pT_Comparison(
//            plot_file_name,
//            subdir_plots,
//            plot_labels,
//            th1d_simple_correction_plot,
//            th1d_linear_regression_plot,
//            th1d_random_forest_plot,
//            th1d_neural_network_plot,
//            th1d_feature_importance_list,
//            true, // show feature list
//            true); // compare to paper plots
//
//
//        // Plot 10 input features, pt_true method B, with comparison to paper
//        sprintf(plot_file_name, "th1d_%s_jet_pt_delta_10feature_ptTrueB_compare.pdf", pt_test_bin_arr[pt_test_bin_index]);
//        std::cout << "----- Generating: " << plot_file_name << " -----" << std::endl;
//
//        sprintf(plot_labels, "Jet p_{T} Delta for %s #left[p_{T}^{true} = p_{T}^{raw} #Sigma p_{T, const. i}^{PYTHIA} / #Sigma p_{T, const. i} #right]; p_{T, ch jet}^{reco} - p_{T, ch jet}^{true} [GeV]; Probability Density", pt_test_bin_title_arr[pt_test_bin_index]);
//
//        sprintf(th1d_feature_importance, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[35]);
//        th1d_feature_importance_list    = (TH1D*) input_file->Get(th1d_feature_importance);
//
//        sprintf(th1d_simple_correction, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[36]);
//        th1d_simple_correction_plot     = (TH1D*) input_file->Get(th1d_simple_correction);
//
//        sprintf(th1d_linear_regression, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[37]);
//        th1d_linear_regression_plot     = (TH1D*) input_file->Get(th1d_linear_regression);
//
//        sprintf(th1d_random_forest, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[38]);
//        th1d_random_forest_plot         = (TH1D*) input_file->Get(th1d_random_forest);
//
//        sprintf(th1d_neural_network, "th1d_%s_%s", pt_test_bin_arr[pt_test_bin_index], th1d_name_arr[39]);
//        th1d_neural_network_plot        = (TH1D*) input_file->Get(th1d_neural_network);
//
//        th1d_simple_correction_plot     ->Scale( 1. / th1d_simple_correction_plot->Integral(),"WIDTH");
//        th1d_linear_regression_plot     ->Scale( 1. / th1d_linear_regression_plot->Integral(),"WIDTH");
//        th1d_random_forest_plot         ->Scale( 1. / th1d_random_forest_plot->Integral(),"WIDTH");
//        th1d_neural_network_plot        ->Scale( 1. / th1d_neural_network_plot->Integral(),"WIDTH");
//
//        Plot_ML_pT_Comparison(
//            plot_file_name,
//            subdir_plots,
//            plot_labels,
//            th1d_simple_correction_plot,
//            th1d_linear_regression_plot,
//            th1d_random_forest_plot,
//            th1d_neural_network_plot,
//            th1d_feature_importance_list,
//            true, // show feature list
//            true); // compare to paper plots
    }
    
//    // Machine Learning Comparison Plots
//    Plot_ML_pT_Comparison(
//        plot_filename,
//        plot_directory,
//        plot_title_xlabel_ylabel,
//        th1d_simple_correction,
//        th1d_linear_regression,
//        th1d_random_forest,
//        th1d_neural_network,
//        th1d_feature_list,
//        show_feature_list,
//        show_paper_plots)

    input_file->Close();
}





//    // Generates a directory for output files
//    char subdir_plots[200];
//    sprintf(subdir_plots, "%s/FastJet_Subtraction", dir_plots);
//    std::__fs::filesystem::create_directories(subdir_plots);
//
//    // Gets histograms from output file
//    // Note: names are the assigned names from histogram creation
//    TH1D* th1d_jet_pt_delta         = (TH1D*) input_file->Get("th1d_jet_pt_delta");
//    TH1D* th1d_const_median_pt      = (TH1D*) input_file->Get("th1d_median_pt_background");
//    TH1D* th1d_jet_const_n          = (TH1D*) input_file->Get("th1d_jet_part_n");
//    TH1D* th1d_jet_const_n_mean     = (TH1D*) input_file->Get("th1d_jet_part_n_mean");
//    TH1D* th1d_jet_const_n_median   = (TH1D*) input_file->Get("th1d_jet_part_n_median");
//    TH1D* th1d_event_rho            = (TH1D*) input_file->Get("th1d_event_rho");
//    TH1D* th1d_jet_truth_low_bin    = (TH1D*) input_file->Get("th1d_jet_truth_low_bin");
//
//    TH1D* th1d_jet_pt_compare       = (TH1D*) input_file->Get("th1d_jet_pt_compare");
//    TH1D* th1d_jet_pt_true          = (TH1D*) input_file->Get("th1d_jet_pt_true");
//    TH1D* th1d_jet_pt_corr          = (TH1D*) input_file->Get("th1d_jet_pt_corr");
//    TH1D* th1d_jet_pt_raw           = (TH1D*) input_file->Get("th1d_jet_pt_raw");
//
//    TH1D* th1d_jet_area_compare     = (TH1D*) input_file->Get("th1d_jet_area_compare");
//    TH1D* th1d_jet_area_compare_norm = (TH1D*) input_file->Get("th1d_jet_area_compare_norm");
//    TH1D* th1d_jet_area_total       = (TH1D*) input_file->Get("th1d_jet_area_total");
//    TH1D* th1d_jet_area_hard        = (TH1D*) input_file->Get("th1d_jet_area_hard");
//    TH1D* th1d_jet_area_soft        = (TH1D*) input_file->Get("th1d_jet_area_soft");
//
//    // Basic Plots
//
//    const int label_arr_size = 2;
//    //char line2[100];
//    //int int1 = sprintf(line2, "2.76 TeV, N_{event} = %i", nEvent);
//    char* label_arr[label_arr_size] = {
//        "FastJet, R = 0.4", "p_{T min, jet} > 5.0 GeV"};
//
//    th1d_plotter(
//        th1d_jet_pt_delta,
//        "FastJet_Subtraction/th1d_jet_pt_delta.pdf",
//        jet_gen_line, jet_gen_mark, mark_circ_fill[0], mark_circ_fill[1],
//        0, false,
//        label_arr, label_arr_size, 0.65, 0.75);
//
//    th1d_plotter(
//        th1d_const_median_pt,
//        "FastJet_Subtraction/th1d_const_median_pt.pdf",
//        jet_gen_line, jet_gen_mark, mark_circ_fill[0], mark_circ_fill[1],
//        0, false,
//        label_arr, label_arr_size, 0.65, 0.75);
//
//    th1d_plotter(
//        th1d_jet_const_n,
//        "FastJet_Subtraction/th1d_jet_const_n.pdf",
//        jet_gen_line, jet_gen_mark, mark_circ_fill[0], mark_circ_fill[1],
//        0, false,
//        label_arr, label_arr_size, 0.65, 0.75);
//
//    th1d_plotter(
//        th1d_jet_const_n_mean,
//        "FastJet_Subtraction/th1d_jet_const_n_mean.pdf",
//        jet_gen_line, jet_gen_mark, mark_circ_fill[0], mark_circ_fill[1],
//        0, false,
//        label_arr, label_arr_size, 0.65, 0.75);
//
//    th1d_plotter(
//        th1d_jet_const_n_median,
//        "FastJet_Subtraction/th1d_jet_const_n_median.pdf",
//        jet_gen_line, jet_gen_mark, mark_circ_fill[0], mark_circ_fill[1],
//        0, false,
//        label_arr, label_arr_size, 0.65, 0.75);
//
//    th1d_plotter(
//        th1d_event_rho,
//        "FastJet_Subtraction/th1d_event_rho.pdf",
//        jet_gen_line, jet_gen_mark, mark_circ_fill[0], mark_circ_fill[1],
//        0, false,
//        label_arr, label_arr_size, 0.65, 0.75);
//
//    th1d_plotter(
//        th1d_jet_truth_low_bin,
//        "FastJet_Subtraction/th1d_jet_truth_low_bin.pdf",
//        jet_gen_line, jet_gen_mark, mark_circ_fill[0], mark_circ_fill[1],
//        0, false,
//        label_arr, label_arr_size, 0.65, 0.75);
//
//
//    // General Plot Variables
//    TCanvas* canvas;
//    TLatex* latex;
//    TLegend* legend;
//    char plot_output[200];
//
//    // Creating Combined Jet pT Plot
//
//    // TCanvas("name", "title", width (px), height (px))
//    canvas = new TCanvas("canvas", "", 800, 600);
//
//    // Turns on ticks, turns off stats
//    gPad->SetTicks();
//    gStyle->SetOptStat(0);
//
//    // Turns on Latex formatting relative to canvas
//    latex = new TLatex();
//    latex->SetNDC(kTRUE);
//
//    // Build histogram plot
//    th1d_jet_pt_true->SetLineColor(jet_red_mark);
//    th1d_jet_pt_true->SetMarkerColor(jet_red_line);
//    th1d_jet_pt_true->SetMarkerStyle(mark_circ_open[0]);
//    th1d_jet_pt_true->SetMarkerSize(mark_circ_open[1]);
//
//    th1d_jet_pt_corr->SetLineColor(jet_vio_mark);
//    th1d_jet_pt_corr->SetMarkerColor(jet_vio_line);
//    th1d_jet_pt_corr->SetMarkerStyle(mark_squa_open[0]);
//    th1d_jet_pt_corr->SetMarkerSize(mark_squa_open[1]);
//
//    th1d_jet_pt_raw->SetLineColor(jet_blu_mark);
//    th1d_jet_pt_raw->SetMarkerColor(jet_blu_line);
//    th1d_jet_pt_raw->SetMarkerStyle(mark_diam_open[0]);
//    th1d_jet_pt_raw->SetMarkerSize(mark_diam_open[1]);
//
//    th1d_jet_pt_compare->SetAxisRange(0, 1.2 * th1d_jet_pt_true->GetMaximum(), "Y");
//
//    th1d_jet_pt_compare->Draw();
//    th1d_jet_pt_raw->Draw("same");
//    th1d_jet_pt_corr->Draw("same");
//    th1d_jet_pt_true->Draw("same");
//
//    // TLegend(xlow, ylow, xup, yup)
//    legend = new TLegend(0.75, 0.65, 0.85, 0.85);
//    legend->AddEntry(th1d_jet_pt_true, "Truth", "lp");
//    legend->AddEntry(th1d_jet_pt_corr, "Corrected", "lp");
//    legend->AddEntry(th1d_jet_pt_raw, "Raw", "lp");
//    legend->SetLineWidth(0);
//    legend->Draw();
//
//    //    latex->DrawLatex(0.60, (ypos + 0.02), "#scale[0.6]{#bf{PYTHIA8 + Thermal Model}}");
//    //    latex->DrawLatex(0.60, (ypos - 0.02), Form("#scale[0.6]{#bf{2.76 TeV, %i events}}", nEvent));
//
//    sprintf(plot_output, "%s/th1d_jet_pt_compare.pdf", subdir_plots);
//    canvas->Print(plot_output);
//    std::cout << "Plotted th1d_jet_pt_compare.pdf" << std::endl;
//    gPad->SetLogy(0);
//
//
//
//    // Creating Combined Jet Area Plot
//
//    // TCanvas("name", "title", width (px), height (px))
//    canvas = new TCanvas("canvas", "", 800, 600);
//
//    // Turns on ticks, turns off stats
//    gPad->SetTicks();
//    gStyle->SetOptStat(0);
//
//    // Turns on Latex formatting relative to canvas
//    latex = new TLatex();
//    latex->SetNDC(kTRUE);
//
//    // Build histogram plot
//    th1d_jet_area_hard->SetLineColor(jet_red_mark);
//    th1d_jet_area_hard->SetMarkerColor(jet_red_line);
//    th1d_jet_area_hard->SetMarkerStyle(mark_triu_open[0]);
//    th1d_jet_area_hard->SetMarkerSize(mark_triu_open[1]);
//
//    th1d_jet_area_total->SetLineColor(jet_vio_mark);
//    th1d_jet_area_total->SetMarkerColor(jet_vio_line);
//    th1d_jet_area_total->SetMarkerStyle(mark_squa_open[0]);
//    th1d_jet_area_total->SetMarkerSize(mark_squa_open[1]);
//
//    th1d_jet_area_soft->SetLineColor(jet_blu_mark);
//    th1d_jet_area_soft->SetMarkerColor(jet_vio_line);
//    th1d_jet_area_soft->SetMarkerStyle(mark_trid_open[0]);
//    th1d_jet_area_soft->SetMarkerSize(mark_trid_open[1]);
//
//    th1d_jet_area_compare->SetAxisRange(0, 1.2 * th1d_jet_area_total->GetMaximum(), "Y");
//
//    th1d_jet_area_compare->Draw();
//    th1d_jet_area_total->Draw("same");
//    th1d_jet_area_soft->Draw("same");
//    th1d_jet_area_hard->Draw("same");
//
//    // TLegend(xlow, ylow, xup, yup)
//    legend = new TLegend(0.75, 0.65, 0.85, 0.85);
//    legend->AddEntry(th1d_jet_area_hard, "Hard Jets", "lp");
//    legend->AddEntry(th1d_jet_area_total, "All Areas", "lp");
//    legend->AddEntry(th1d_jet_area_soft, "Soft Jets", "lp");
//    legend->SetLineWidth(0);
//    legend->Draw();
//
//    //    latex->DrawLatex(0.60, (ypos + 0.02), "#scale[0.6]{#bf{PYTHIA8 + Thermal Model}}");
//    //    latex->DrawLatex(0.60, (ypos - 0.02), Form("#scale[0.6]{#bf{2.76 TeV, %i events}}", nEvent));
//
//    sprintf(plot_output, "%s/th1d_jet_area_compare.pdf", subdir_plots);
//    canvas->Print(plot_output);
//    std::cout << "Plotted th1d_jet_area_compare.pdf" << std::endl;
//    gPad->SetLogy(0);
//
//     // Rescale Area Hists
//
//    th1d_jet_area_hard  ->Scale( 1. / th1d_jet_area_hard->Integral(),"WIDTH");
//    th1d_jet_area_soft  ->Scale( 1. / th1d_jet_area_soft->Integral(),"WIDTH");
//    th1d_jet_area_total ->Scale( 1. / th1d_jet_area_total->Integral(),"WIDTH");
//
//    // Print Scaled Plot
//
//    th1d_jet_area_compare_norm->SetAxisRange(0, 1.2 * th1d_jet_area_hard->GetMaximum(), "Y");
//
//    th1d_jet_area_compare_norm->Draw();
//    th1d_jet_area_total->Draw("same");
//    th1d_jet_area_soft->Draw("same");
//    th1d_jet_area_hard->Draw("same");
//
//    // TLegend(xlow, ylow, xup, yup)
//    legend->Draw();
//
//    //    latex->DrawLatex(0.60, (ypos + 0.02), "#scale[0.6]{#bf{PYTHIA8 + Thermal Model}}");
//    //    latex->DrawLatex(0.60, (ypos - 0.02), Form("#scale[0.6]{#bf{2.76 TeV, %i events}}", nEvent));
//
//    sprintf(plot_output, "%s/th1d_jet_area_compare_norm.pdf", subdir_plots);
//    canvas->Print(plot_output);
//    std::cout << "Plotted th1d_jet_area_compare_norm.pdf" << std::endl;
//    gPad->SetLogy(0);

