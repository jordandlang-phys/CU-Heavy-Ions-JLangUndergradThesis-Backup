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
    tf1_paper_simple_correction_fit = new TF1("tf1_paper_simple_correction_fit", "gaus", -40., 40.);
    tf1_paper_simple_correction_fit->SetParameters(.1, 0., 5.);
    tgraph_paper_simple_correction->Fit(tf1_paper_simple_correction_fit,"RNL");
    
    double tgraph_paper_nn_x[16] = {-15., -13., -11., -9., -7., -5., -3., -1., 1., 3., 5., 7., 9., 11., 12., 15.};
    double tgraph_paper_nn_y[16] = {0.0, 0.0, 0.0007272727272727430, 0.003636363636363630, 0.015272727272727300, 0.03563636363636370, 0.07272727272727270, 0.10581818181818200, 0.11054545454545500, 0.07854545454545460, 0.04290909090909090, 0.02036363636363640, 0.00872727272727275, 0.0029090909090909200, 0.00036363636363637200, 0.0};
    
    TGraph* tgraph_paper_neural_network = new TGraph(16, tgraph_paper_nn_x, tgraph_paper_nn_y);
    tf1_paper_neural_network_fit = new TF1("tf1_paper_neural_network_fit", "gaus", -40., 40.);
    tf1_paper_neural_network_fit->SetParameters(.1, 0., 5.);
    tgraph_paper_neural_network->Fit(tf1_paper_neural_network_fit,"RNL");
    
    double tgraph_paper_rf_x[16] = {-15., -13., -11., -9., -7., -5., -3., -1., 1., 3., 5., 7., 9., 11., 12., 15.};
    double tgraph_paper_rf_y[16] = {0.0, 0.0, 0.0014545454545454600, 0.0058181818181818300, 0.02036363636363640, 0.04618181818181820, 0.07927272727272730, 0.09854545454545460, 0.09600000000000000, 0.0701818181818182, 0.042181818181818200, 0.023272727272727300, 0.009454545454545470, 0.003636363636363630, 0.0007272727272727430, 0.0};
    
    TGraph* tgraph_paper_random_forest = new TGraph(16, tgraph_paper_rf_x, tgraph_paper_rf_y);
    tf1_paper_random_forest_fit = new TF1("tf1_paper_random_forest_fit", "gaus", -40., 40.);
    tf1_paper_random_forest_fit->SetParameters(.1, 0., 5.);
    tgraph_paper_random_forest->Fit(tf1_paper_random_forest_fit,"RNL");
    
    double tgraph_paper_lr_x[16] = {-15., -13., -11., -9., -7., -5., -3., -1., 1., 3., 5., 7., 9., 11., 12., 15.};
    double tgraph_paper_lr_y[16] = {0.0, 0.00036363636363637200, 0.002545454545454570, 0.009818181818181840, 0.02872727272727270, 0.052727272727272700, 0.08145454545454550, 0.09454545454545460, 0.08836363636363640, 0.06618181818181820, 0.03927272727272730, 0.021090909090909100, 0.009818181818181840, 0.003636363636363630, 0.0007272727272727430, 0.0};
    
    TGraph* tgraph_paper_linear_regression = new TGraph(16, tgraph_paper_lr_x, tgraph_paper_lr_y);
    tf1_paper_linear_regression_fit = new TF1("tf1_paper_linear_regression_fit", "gaus", -40., 40.);
    tf1_paper_linear_regression_fit->SetParameters(.1, 0., 5.);
    tgraph_paper_linear_regression->Fit(tf1_paper_linear_regression_fit,"RNL");
}

void Plot_ML_pT_Comparison(
    char* plot_filename,
    char* plot_directory,
    char* plot_title_xlabel_ylabel,
    TH1D* th1d_simple_correction,
    TH1D* th1d_linear_regression,
    TH1D* th1d_random_forest,
    TH1D* th1d_neural_network,
    TH1D* th1d_feature_list,
    double y_max,
    bool bool_showLegend,
    char* showFeatures,
    bool bool_compareToPaper
    ) {
    
    char feature_list_1[2][100];
    sprintf(feature_list_1[0], "#scale[0.6]{Feature Importance}");
    sprintf(feature_list_1[1], "#scale[0.6]{#bf{Jet p_{T, raw}:}}");
    
    char feature_list_3[4][100];
    sprintf(feature_list_3[0], "#scale[0.6]{Feature Importance}");
    sprintf(feature_list_3[1], "#scale[0.6]{#bf{Jet p_{T, raw}:}}");
    sprintf(feature_list_3[2], "#scale[0.6]{#bf{Jet area:}}");
    sprintf(feature_list_3[3], "#scale[0.6]{#bf{Jet #rho:}}");

    char feature_list_8[9][100];
    sprintf(feature_list_8[0], "#scale[0.6]{Feature Importance}");
    sprintf(feature_list_8[1], "#scale[0.6]{#bf{Jet p_{T, raw}:}}");
    sprintf(feature_list_8[2], "#scale[0.6]{#bf{Jet p_{T, corr}:}}");
    sprintf(feature_list_8[3], "#scale[0.6]{#bf{Mean Const. p_{T}:}}");
    sprintf(feature_list_8[4], "#scale[0.6]{#bf{Median Const. p_{T}:}}");
    sprintf(feature_list_8[5], "#scale[0.6]{#bf{p_{T, Const.}^{1}:}}");
    sprintf(feature_list_8[6], "#scale[0.6]{#bf{p_{T, Const.}^{2}:}}");
    sprintf(feature_list_8[7], "#scale[0.6]{#bf{p_{T, Const.}^{3}:}}");
    sprintf(feature_list_8[8], "#scale[0.6]{#bf{p_{T, Const.}^{4}:}}");

    char feature_list_12[13][100];
    sprintf(feature_list_12[0], "#scale[0.6]{Feature Importance}");
    sprintf(feature_list_12[1], "#scale[0.6]{#bf{Jet p_{T, raw}:}}");
    sprintf(feature_list_12[2], "#scale[0.6]{#bf{Jet p_{T, corr}:}}");
    sprintf(feature_list_12[3], "#scale[0.6]{#bf{Jet Mass:}}");
    sprintf(feature_list_12[4], "#scale[0.6]{#bf{Jet Area:}}");
    sprintf(feature_list_12[5], "#scale[0.6]{#bf{N_{const}:}}");
    sprintf(feature_list_12[6], "#scale[0.6]{#bf{Mean Const. p_{T}:}}");
    sprintf(feature_list_12[7], "#scale[0.6]{#bf{p_{T, Const.}^{1}:}}");
    sprintf(feature_list_12[8], "#scale[0.6]{#bf{p_{T, Const.}^{2}:}}");
    sprintf(feature_list_12[9], "#scale[0.6]{#bf{p_{T, Const.}^{3}:}}");
    sprintf(feature_list_12[10], "#scale[0.6]{#bf{p_{T, Const.}^{4}:}}");
    sprintf(feature_list_12[11], "#scale[0.6]{#bf{Jet y:}}");
    sprintf(feature_list_12[12], "#scale[0.6]{#bf{Jet #rho:}}");
    
    double x_min  = th1d_simple_correction->GetXaxis()->GetXmin();
    double x_max  = th1d_simple_correction->GetXaxis()->GetXmax();
    double x_bins = th1d_simple_correction->GetNbinsX();
    double max_arr[4] = {th1d_simple_correction->GetMaximum(), th1d_neural_network->GetMaximum(), th1d_random_forest->GetMaximum(), th1d_linear_regression->GetMaximum()};
    if (y_max < 0 ) {
        for ( int i = 0 ; i < 4 ; i++ ) {
            if ( y_max < max_arr[i] ) y_max = max_arr[i];
        }
        if (bool_compareToPaper) y_max = 0.15;
    }
    
    TH1D* th1d_canvas_plot = new TH1D("th1d_canvas_plot", plot_title_xlabel_ylabel, x_bins, x_min, x_max);
    th1d_canvas_plot->SetAxisRange(0, 1.2 * y_max, "Y");
    
    // Fit to Each Plot
    TF1* tf1_simple_correction_fit  = new TF1("tf1_simple_correction_fit",  "gaus", x_min, x_max);
    TF1* tf1_linear_regression_fit  = new TF1("tf1_linear_regression_fit",  "gaus", x_min, x_max);
    TF1* tf1_random_forest_fit      = new TF1("tf1_random_forest_fit",      "gaus", x_min, x_max);
    TF1* tf1_neural_network_fit     = new TF1("tf1_neural_network_fit",     "gaus", x_min, x_max);

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
    
    th1d_linear_regression->SetLineColor(jet_blu_line);
    th1d_linear_regression->SetMarkerColor(jet_blu_mark);
    th1d_linear_regression->SetMarkerStyle(mark_squa_open[0]);
    th1d_linear_regression->SetMarkerSize(mark_squa_open[1]);
    tf1_linear_regression_fit->SetLineColor(jet_blu_line);
    tf1_linear_regression_fit->SetLineWidth(1);
    
    th1d_random_forest->SetLineColor(jet_tea_line);
    th1d_random_forest->SetMarkerColor(jet_tea_mark);
    th1d_random_forest->SetMarkerStyle(mark_diam_open[0]);
    th1d_random_forest->SetMarkerSize(mark_diam_open[1]);
    tf1_random_forest_fit->SetLineColor(jet_tea_line);
    tf1_random_forest_fit->SetLineWidth(1);

    th1d_neural_network->SetLineColor(jet_blk_line);
    th1d_neural_network->SetMarkerColor(jet_blk_mark);
    th1d_neural_network->SetMarkerStyle(mark_star_open[0]);
    th1d_neural_network->SetMarkerSize(mark_star_open[1]);
    tf1_neural_network_fit->SetLineColor(jet_blk_line);
    tf1_neural_network_fit->SetLineWidth(1);
    
    // Draws the histograms and lines
    th1d_canvas_plot->Draw();
    
    // Draws results from ML paper
    if (bool_compareToPaper) {
        Fit_Paper_Data();
        
        tf1_paper_neural_network_fit->SetLineColor(kGray+3);
        tf1_paper_neural_network_fit->SetLineWidth(1);
        tf1_paper_neural_network_fit->SetLineStyle(3);
    
        tf1_paper_random_forest_fit->SetLineColor(jet_tea_line);
        tf1_paper_random_forest_fit->SetLineWidth(1);
        tf1_paper_random_forest_fit->SetLineStyle(3);
    
        tf1_paper_linear_regression_fit->SetLineColor(jet_blu_line);
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
    
    // Draws plots and fit lines
    tf1_simple_correction_fit   ->Draw("same");
    tf1_neural_network_fit      ->Draw("same");
    tf1_random_forest_fit       ->Draw("same");
    tf1_linear_regression_fit   ->Draw("same");
    th1d_simple_correction      ->Draw("same");
    th1d_neural_network         ->Draw("same");
    th1d_random_forest          ->Draw("same");
    th1d_linear_regression      ->Draw("same");

    // Draws Legend as TLegend(xmin, ymin, xmax, ymax)
    if ( bool_showLegend ) {
        TLegend* legend;
        if (bool_compareToPaper) legend = new TLegend(0.15, 0.45, 0.40, 0.85);
        else legend = new TLegend(0.15, 0.65, 0.40, 0.85);
        legend->AddEntry(th1d_neural_network,       "Neural Network", "lp");
        legend->AddEntry(th1d_random_forest,        "Random Forest", "lp");
        legend->AddEntry(th1d_linear_regression,    "Linear Regression", "lp");
        legend->AddEntry(th1d_simple_correction,    "Area Correction", "lp");
        if (bool_compareToPaper) {
            legend->AddEntry(tf1_paper_neural_network_fit,      "NN (Paper)", "l");
            legend->AddEntry(tf1_paper_random_forest_fit,       "RF (Paper)", "l");
            legend->AddEntry(tf1_paper_linear_regression_fit,   "LR (Paper)", "l");
            legend->AddEntry(tf1_paper_simple_correction_fit,   "AC (Paper)", "l");
        }
        legend->SetLineWidth(0);
        legend->SetFillStyle(0);
        legend->Draw();
    }
    
    // Draws Feature List
    if ( showFeatures == "sigmas" ) {
        char sigma_vals[5][100];
        char sigma_list[5][100];
        
        sprintf(sigma_list[0], "#scale[0.6]{Fit Widths}");
        sprintf(sigma_list[1], "#scale[0.6]{#bf{#sigma_{Neural Network}}}");
        sprintf(sigma_list[2], "#scale[0.6]{#bf{#sigma_{Random Forest}}}");
        sprintf(sigma_list[3], "#scale[0.6]{#bf{#sigma_{Linear Regression}}}");
        sprintf(sigma_list[4], "#scale[0.6]{#bf{#sigma_{Area Correction}}}");
        
        sprintf(sigma_vals[0], "");
        sprintf(sigma_vals[1], "#scale[0.6]{#bf{%1.3f}}", tf1_neural_network_fit    ->GetParameter(2));
        sprintf(sigma_vals[2], "#scale[0.6]{#bf{%1.3f}}", tf1_random_forest_fit     ->GetParameter(2));
        sprintf(sigma_vals[3], "#scale[0.6]{#bf{%1.3f}}", tf1_linear_regression_fit ->GetParameter(2));
        sprintf(sigma_vals[4], "#scale[0.6]{#bf{%1.3f}}", tf1_simple_correction_fit ->GetParameter(2));
        
        for ( int i = 0 ; i <= 5 ; i++ ) {
            latex->DrawLatex(0.63, (.813 - 0.0494 * i), sigma_list[i]);
        }
        for ( int i = 0 ; i <= 5 ; i++ ) {
            latex->DrawLatex(0.80, (.813 - 0.0494 * i), sigma_vals[i]);
        }
    }
    if ( showFeatures == "features" ) {
        int  feature_count = th1d_feature_list->GetNbinsX();
        char feature_vals_1[2][100];
        char feature_vals_3[4][100];
        char feature_vals_8[9][100];
        char feature_vals_12[13][100];
        
        if (feature_count == 1) {
            sprintf(feature_vals_1[0], "#scale[0.6]{}");
            sprintf(feature_vals_1[1], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(1));
            
            for ( int i = 0 ; i <= 1 ; i++ ) {
                latex->DrawLatex(0.63, (.813 - 0.0494 * i), feature_list_1[i]);
            }
            for ( int i = 0 ; i <= 1 ; i++ ) {
                latex->DrawLatex(0.80, (.813 - 0.0494 * i), feature_vals_3[i]);
            }
        }
        
        if (feature_count == 3) {
            sprintf(feature_vals_3[0], "#scale[0.6]{}");
            sprintf(feature_vals_3[1], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(1));
            sprintf(feature_vals_3[2], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(2));
            sprintf(feature_vals_3[3], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(3));
            
            for ( int i = 0 ; i <= 3 ; i++ ) {
                latex->DrawLatex(0.63, (.813 - 0.0494 * i), feature_list_3[i]);
            }
            for ( int i = 0 ; i <= 3 ; i++ ) {
                latex->DrawLatex(0.80, (.813 - 0.0494 * i), feature_vals_3[i]);
            }
        }
        if (feature_count == 8) {
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
        if (feature_count == 12) {
            sprintf(feature_vals_12[0], "#scale[0.6]{}");
            sprintf(feature_vals_12[1], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(1));
            sprintf(feature_vals_12[2], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(2));
            sprintf(feature_vals_12[3], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(3));
            sprintf(feature_vals_12[4], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(4));
            sprintf(feature_vals_12[5], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(5));
            sprintf(feature_vals_12[6], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(6));
            sprintf(feature_vals_12[7], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(7));
            sprintf(feature_vals_12[8], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(8));
            sprintf(feature_vals_12[9], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(7));
            sprintf(feature_vals_12[10], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(8));
            sprintf(feature_vals_12[11], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(7));
            sprintf(feature_vals_12[12], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(8));
            
            for ( int i = 0 ; i <= 12 ; i++ ) {
                latex->DrawLatex(0.63, (.813 - 0.0494 * i), feature_list_12[i]);
            }
            for ( int i = 0 ; i <= 12 ; i++ ) {
                latex->DrawLatex(0.80, (.813 - 0.0494 * i), feature_vals_12[i]);
            }
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





void Jet_ML_Plotter(
    char* input_file_name,
    char* th1d_name_arr[5],
    char* plot_file_name,
    char* plot_dir_name,
    float pt_min, float pt_max,
    bool bool_showLegend,
    bool bool_showFeatures,
    bool bool_normalize,
    bool bool_compareToPaper,
    char* truth_source
    ) {
    
    // Opens and reads the Root output file
    char input_file_path[200];
    sprintf(input_file_path, "%s/%s", dir_data, input_file_name);
    TFile* input_file = new TFile(input_file_path, "READ");
    std::cout << "Input file read." << std::endl;
    
    char subdir_plots [200];
    sprintf(subdir_plots, "%s/MachineLearning/%s", dir_plots, plot_dir_name);
    std::__fs::filesystem::create_directories(subdir_plots);
    
    // EXAMINING FEATURES
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
    
    // Plot 3 input features, pt_true = PYTHIA, with comparison to paper
    std::cout << "----- Generating: " << plot_file_name << " -----" << std::endl;
    
    if ( bool_normalize ) {
        if ( truth_source == "Pythia" ) {
            sprintf(plot_labels, "Jet p_{T} Delta for %i-%i GeV #left[p_{T}^{true} = p_{T, ch jet}^{PYTHIA} #right]; p_{T, ch jet}^{reco} - p_{T, ch jet}^{true} [GeV]; Probability Density", int(pt_min), int(pt_max));
        }
        else if ( truth_source == "Paper" ) {
            sprintf(plot_labels, "Jet p_{T} Delta for %i-%i GeV #left[p_{T}^{true} = p_{T}^{raw} #Sigma p_{T, const. i}^{PYTHIA} / #Sigma p_{T, const. i} #right]; p_{T, ch jet}^{reco} - p_{T, ch jet}^{true} [GeV]; Probability Density", int(pt_min), int(pt_max));
        }
    }
    else {
        if ( truth_source == "Pythia" ) {
            sprintf(plot_labels, "Jet p_{T} Delta for %i-%i GeV #left[p_{T}^{true} = p_{T, ch jet}^{PYTHIA} #right]; p_{T, ch jet}^{reco} - p_{T, ch jet}^{true} [GeV]; N_{events}", int(pt_min), int(pt_max));
        }
        else if ( truth_source == "Paper" ) {
            sprintf(plot_labels, "Jet p_{T} Delta for %i-%i GeV #left[p_{T}^{true} = p_{T}^{raw} #Sigma p_{T, const. i}^{PYTHIA} / #Sigma p_{T, const. i} #right]; p_{T, ch jet}^{reco} - p_{T, ch jet}^{true} [GeV]; N_{events}", int(pt_min), int(pt_max));
        }
    }
    
    sprintf(th1d_feature_importance, "%s", th1d_name_arr[0]);
    th1d_feature_importance_list    = (TH1D*) input_file->Get(th1d_feature_importance);
    std::cout << "Imported " << th1d_name_arr[0] << std::endl;
    sprintf(th1d_simple_correction, "%s", th1d_name_arr[1]);
    th1d_simple_correction_plot     = (TH1D*) input_file->Get(th1d_simple_correction);
    
    TCanvas* canvas_temp = new TCanvas("canvas", "", 1000, 600);
    th1d_simple_correction_plot->Draw();
    
    sprintf(th1d_linear_regression, "%s", th1d_name_arr[2]);
    th1d_linear_regression_plot     = (TH1D*) input_file->Get(th1d_linear_regression);

    sprintf(th1d_random_forest, "%s", th1d_name_arr[3]);
    th1d_random_forest_plot         = (TH1D*) input_file->Get(th1d_random_forest);

    sprintf(th1d_neural_network, "%s", th1d_name_arr[4]);
    th1d_neural_network_plot        = (TH1D*) input_file->Get(th1d_neural_network);
    
    std::cout << "Histograms imported." << std::endl;
    
    if ( bool_normalize || bool_compareToPaper ) {
        th1d_simple_correction_plot     ->Scale( 1. / th1d_simple_correction_plot->Integral(),"WIDTH");
        th1d_linear_regression_plot     ->Scale( 1. / th1d_linear_regression_plot->Integral(),"WIDTH");
        th1d_random_forest_plot         ->Scale( 1. / th1d_random_forest_plot->Integral(),"WIDTH");
        th1d_neural_network_plot        ->Scale( 1. / th1d_neural_network_plot->Integral(),"WIDTH");
    }
    
    std::cout << "Histograms normalized." << std::endl;
    
    char* char_showFeatures = "";
    if ( bool_showFeatures ) char_showFeatures = "features";
    
    Plot_ML_pT_Comparison(
        plot_file_name,
        subdir_plots,
        plot_labels,
        th1d_simple_correction_plot,
        th1d_linear_regression_plot,
        th1d_random_forest_plot,
        th1d_neural_network_plot,
        th1d_feature_importance_list,
        -1,
        bool_showLegend, // show legend on plot
        char_showFeatures, // show feature list
        bool_compareToPaper); // compare to paper plots
    
    input_file->Close();
        
    delete input_file;
}

void Jet_Test_Range_Plotter(
    char* input_file_name,
    char* plot_file_name,
    char* plot_dir_name,
    float train_pt_min, float train_pt_max,
    float test_pt_min, float test_pt_max,
    char* truth_source
    ) {
    
    float x_min = train_pt_min - 20.;
    float x_max = train_pt_max + 20.;
    int x_bins = int((x_max - x_min) / 1);
    
    // Opens and reads the Root output file
    char input_file_path[200];
    sprintf(input_file_path, "%s/%s", dir_data, input_file_name);
    TFile* input_file = new TFile(input_file_path, "READ");
    TTree* input_tree = (TTree*) input_file->Get("Tree_ML");
    std::cout << "Input file read." << std::endl;
    
    char subdir_plots [200];
    sprintf(subdir_plots, "%s/MachineLearning/%s", dir_plots, plot_dir_name);
    std::__fs::filesystem::create_directories(subdir_plots);
    
    float jet_area;
    float jet_pt_raw;
    float jet_pt_true;
    float jet_pt_corr;
    float jet_pt_ml_lr;
    float jet_pt_ml_rf;
    float jet_pt_ml_nn;
    
    input_tree->SetBranchAddress("jet_area", &jet_area);
    input_tree->SetBranchAddress("jet_pt_raw", &jet_pt_raw);
    input_tree->SetBranchAddress("jet_pt_true", &jet_pt_true);
    input_tree->SetBranchAddress("jet_pt_corr", &jet_pt_corr);
    input_tree->SetBranchAddress("jet_pt_ml_lr", &jet_pt_ml_lr);
    input_tree->SetBranchAddress("jet_pt_ml_rf", &jet_pt_ml_rf);
    input_tree->SetBranchAddress("jet_pt_ml_nn", &jet_pt_ml_nn);
    
    std::cout << "MADE IT HERE!" << std::endl;
    
    char plot_labels[200];
    if ( truth_source == "Pythia" ) {
        sprintf(plot_labels, "Jet p_{T} Delta for %i < p_{T, ch jet}^{true} < %i GeV (Trained on %i < p_{T, ch jet}^{true} < %i GeV); p_{T, ch jet}^{reco} [GeV]; Probability Density", int(test_pt_min), int(test_pt_max), int(train_pt_min), int(train_pt_max));
    }
    else if ( truth_source == "Paper" ) {
        sprintf(plot_labels, "Jet p_{T} Delta for %i < p_{T, ch jet}^{true} < %i GeV (Trained on %i < p_{T, ch jet}^{true} < %i GeV); p_{T, ch jet}^{reco} [GeV]; Probability Density", int(test_pt_min), int(test_pt_max), int(train_pt_min), int(train_pt_max));
    }
    
    TH1D* th1d_simple_correction_plot   = new TH1D("th1d_simple_correction_plot", "Machine Learning vs. Simple Correction; p_{T, ch jet}^{reco} [GeV]; Probability Density", x_bins, x_min, x_max);
    TH1D* th1d_linear_regression_plot   = new TH1D("th1d_linear_regression_plot", "Machine Learning vs. Simple Correction; p_{T, ch jet}^{reco} [GeV]; Probability Density", x_bins, x_min, x_max);
    TH1D* th1d_random_forest_plot       = new TH1D("th1d_random_forest_plot", "Machine Learning vs. Simple Correction; p_{T, ch jet}^{reco} [GeV]; Probability Density", x_bins, x_min, x_max);
    TH1D* th1d_neural_network_plot      = new TH1D("th1d_neural_network_plot", "Machine Learning vs. Simple Correction; p_{T, ch jet}^{reco} [GeV]; Probability Density", x_bins, x_min, x_max);
    TH1D* th1d_truth_plot               = new TH1D("th1d_truth_plot", "Machine Learning vs. Simple Correction; p_{T, ch jet}^{reco} [GeV]; Probability Density", x_bins, x_min, x_max);
    
    // Empty placeholder
    TH1D* th1d_feature_importance_list = new TH1D("", "", 1, 0, 1);
    
    // Fill histograms
    int event_n = input_tree->GetEntries();
    
    for ( int e = 0 ; e < event_n ; e++ ) {
        input_tree->GetEntry(e);
        if ( test_pt_min < jet_pt_true && test_pt_max > jet_pt_true ){
            th1d_simple_correction_plot ->Fill(jet_pt_corr);
            th1d_linear_regression_plot ->Fill(jet_pt_ml_lr);
            th1d_random_forest_plot     ->Fill(jet_pt_ml_rf);
            th1d_neural_network_plot    ->Fill(jet_pt_ml_nn);
            th1d_truth_plot             ->Fill(jet_pt_true);
        }
    }
    
    th1d_simple_correction_plot     ->Scale( 1. / th1d_simple_correction_plot->Integral(),"WIDTH");
    th1d_linear_regression_plot     ->Scale( 1. / th1d_linear_regression_plot->Integral(),"WIDTH");
    th1d_random_forest_plot         ->Scale( 1. / th1d_random_forest_plot->Integral(),"WIDTH");
    th1d_neural_network_plot        ->Scale( 1. / th1d_neural_network_plot->Integral(),"WIDTH");
    th1d_truth_plot                 ->Scale( 1. / th1d_truth_plot->Integral(),"WIDTH");
    
    Plot_ML_pT_Comparison(
        plot_file_name,
        subdir_plots,
        plot_labels,
        th1d_simple_correction_plot,
        th1d_linear_regression_plot,
        th1d_random_forest_plot,
        th1d_neural_network_plot,
        th1d_feature_importance_list,
        0.22, // ymax
        true, // show legend on plot
        "sigmas", // show feature list
        false); // compare to paper plots
    
}




void Jet_Truth_Plotter(
    char* input_file_name,
    char* plot_file_name,
    char* plot_dir_name,
    float train_pt_min, float train_pt_max,
    char* truth_source
    ) {
    
    float x_min = train_pt_min - 10.;
    float x_max = train_pt_max + 10.;
    int x_bins = int((x_max - x_min) / 1);
    
    // Opens and reads the Root output file
    char input_file_path[200];
    sprintf(input_file_path, "%s/%s", dir_data, input_file_name);
    TFile* input_file = new TFile(input_file_path, "READ");
    TTree* input_tree = (TTree*) input_file->Get("Tree_ML");
    std::cout << "Input file read." << std::endl;
    
    char subdir_plots [200];
    sprintf(subdir_plots, "%s/MachineLearning/%s", dir_plots, plot_dir_name);
    std::__fs::filesystem::create_directories(subdir_plots);
    
    float jet_pt_true;
    input_tree->SetBranchAddress("jet_pt_true", &jet_pt_true);
    
    std::cout << "MADE IT HERE!" << std::endl;
    
    char plot_labels[200];
    if ( truth_source == "Pythia" ) {
        sprintf(plot_labels, "Distribution of Jet p_{T}^{True} for %i < p_{T, ch jet}^{true} < %i GeV; p_{T, ch jet}^{reco} [GeV]; Probability Density", int(train_pt_min), int(train_pt_max));
    }
    else if ( truth_source == "Paper" ) {
        sprintf(plot_labels, "Distribution of Jet p_{T}^{True} for %i < p_{T, ch jet}^{true} < %i GeV; p_{T, ch jet}^{reco} [GeV]; Probability Density", int(train_pt_min), int(train_pt_max));
    }
    
    TH1D* th1d_truth_plot   = new TH1D("th1d_pt_true_plot", "Distribution of p_{T}^{True}; p_{T, ch jet}^{True} [GeV]; Probability Density", x_bins, x_min, x_max);
    
    // Fill histograms
    int event_n = input_tree->GetEntries();
    
    for ( int e = 0 ; e < event_n ; e++ ) {
        input_tree->GetEntry(e);
        th1d_truth_plot ->Fill(jet_pt_true);
    }
    
    th1d_truth_plot ->Scale( 1. / th1d_truth_plot->Integral(),"WIDTH");
    
    TH1D* th1d_canvas_plot = new TH1D("th1d_canvas_plot", plot_labels, x_bins, x_min, x_max);
    
    // TCanvas("name", "title", width (px), height (px))
    TCanvas* canvas = new TCanvas("canvas", "", 1000, 600);
    th1d_canvas_plot->SetAxisRange(0, .15, "Y");

    // Turns on ticks, turns off stats
    gPad->SetTicks();
    gStyle->SetOptStat(0);

    // Turns on Latex formatting relative to canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(kTRUE);
    
    // Build histogram plot
    th1d_truth_plot->SetLineColor(jet_blk_line);
    th1d_truth_plot->SetMarkerColor(jet_blk_mark);
    th1d_truth_plot->SetMarkerStyle(mark_circ_open[0]);
    th1d_truth_plot->SetMarkerSize(mark_circ_open[1]);
    
    // Draws the histograms and lines
    th1d_canvas_plot->Draw();
    
    th1d_truth_plot ->Draw("same");

    // Draws Legend as TLegend(xmin, ymin, xmax, ymax)
    TLegend* legend;
    legend = new TLegend(0.15, 0.80, 0.3, 0.85);
    legend->AddEntry(th1d_truth_plot, "p_{T}^{True}", "p");
    legend->SetLineWidth(0);
    legend->SetFillStyle(0);
    legend->Draw();
    
    // Prints out the plot
    char plot_output[300];
    sprintf(plot_output, "%s/MachineLearning/%s/%s", dir_plots, plot_dir_name, plot_file_name);
    canvas->Print(plot_output);
    std::cout << "Plotted " << plot_file_name << std::endl;
    gPad->SetLogy(0);
    
    return 0;
}



void Jet_ptTrue_Distribution_Plotter (
    char* input_file_name,
    char* plot_tree_name,
    char* plot_dir_name,
    float train_pt_min, float train_pt_max
    ) {

    float x_min = train_pt_min;
    float x_max = train_pt_max;
    int x_bins = int((x_max - x_min) / 1);
    char plot_file_name[100];
    sprintf(plot_file_name, "%s.pdf", plot_tree_name);

    // Opens and reads the Root output file
    char input_file_path[200];
    sprintf(input_file_path, "%s/%s", dir_data, input_file_name);
    TFile* input_file = new TFile(input_file_path, "READ");
    TTree* input_tree = (TTree*) input_file->Get("Tree_ML");
    std::cout << "Input file read." << std::endl;

    char subdir_plots [200];
    sprintf(subdir_plots, "%s/MachineLearning/%s", dir_plots, plot_dir_name);
    std::__fs::filesystem::create_directories(subdir_plots);

    TH1D* th1d_truth_distribution_plot     = (TH1D*) input_file->Get(plot_tree_name);

    std::cout << "MADE IT HERE!" << std::endl;

    // TCanvas("name", "title", width (px), height (px))
    TCanvas* canvas = new TCanvas("canvas", "", 1000, 600);
    th1d_truth_distribution_plot->SetAxisRange(100, 1000000, "Y");

    // Turns on ticks, turns off stats
    gPad->SetTicks();
    gStyle->SetOptStat(0);
    gPad->SetLogy(1);

    // Turns on Latex formatting relative to canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(kTRUE);

    // Build histogram plot
    th1d_truth_distribution_plot->SetLineColor(jet_blk_line);
    th1d_truth_distribution_plot->SetMarkerColor(jet_blk_mark);
    th1d_truth_distribution_plot->SetMarkerStyle(mark_circ_open[0]);
    th1d_truth_distribution_plot->SetMarkerSize(mark_circ_open[1]);

    // Draws the histograms and lines
    th1d_truth_distribution_plot->Draw();

    // Draws Legend as TLegend(xmin, ymin, xmax, ymax)
    TLegend* legend;
    legend = new TLegend(0.15, 0.80, 0.3, 0.85);
    legend->AddEntry(th1d_truth_distribution_plot, "p_{T}^{True}", "lp");
    legend->SetLineWidth(0);
    legend->SetFillStyle(0);
    legend->Draw();

    // Prints out the plot
    char plot_output[300];
    sprintf(plot_output, "%s/MachineLearning/%s/%s", dir_plots, plot_dir_name, plot_file_name);
    canvas->Print(plot_output);
    std::cout << "Plotted " << plot_file_name << std::endl;
    gPad->SetLogy(0);

    return 0;
}



void Jet_ML_Analyzer_ROOT() {
    
    // --- Training with NO Bias (B0)
    
    // 1 Feature Test
    
    char* th1d_name_arr_1[5] = {
        "th1d_Tree_40_60_B0_Test_1feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B0_Test_1feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B0_Test_1feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B0_Test_1feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B0_Test_1feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B0_Train.root", // input file name
        th1d_name_arr_1, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_1feat_B0_Train_40_60_B0_Test.pdf", // plot file name
        "1feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    char* th1d_name_arr_2[5] = {
        "th1d_Tree_40_60_B4_Test_1feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B4_Test_1feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B4_Test_1feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B4_Test_1feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B4_Test_1feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B0_Train.root", // input file name
        th1d_name_arr_2, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_1feat_B0_Train_40_60_B4_Test.pdf", // plot file name
        "1feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    char* th1d_name_arr_3[5] = {
        "th1d_Tree_40_60_B8_Test_1feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B8_Test_1feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B8_Test_1feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B8_Test_1feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B8_Test_1feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B0_Train.root", // input file name
        th1d_name_arr_3, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_1feat_B0_Train_40_60_B4_Test.pdf", // plot file name
        "1feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    // 3 Feature Test
    
    char* th1d_name_arr_4[5] = {
        "th1d_Tree_40_60_B0_Test_3feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B0_Test_3feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B0_Test_3feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B0_Test_3feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B0_Test_3feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B0_Train.root", // input file name
        th1d_name_arr_4, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_3feat_B0_Train_40_60_B0_Test.pdf", // plot file name
        "3feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    char* th1d_name_arr_5[5] = {
        "th1d_Tree_40_60_B4_Test_3feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B4_Test_3feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B4_Test_3feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B4_Test_3feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B4_Test_3feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B0_Train.root", // input file name
        th1d_name_arr_5, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_3feat_B0_Train_40_60_B4_Test.pdf", // plot file name
        "3feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    char* th1d_name_arr_6[5] = {
        "th1d_Tree_40_60_B8_Test_3feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B8_Test_3feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B8_Test_3feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B8_Test_3feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B8_Test_3feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B0_Train.root", // input file name
        th1d_name_arr_6, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_3feat_B0_Train_40_60_B8_Test.pdf", // plot file name
        "3feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    // 12 Feature Test
    
    char* th1d_name_arr_7[5] = {
        "th1d_Tree_40_60_B0_Test_12feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B0_Test_12feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B0_Test_12feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B0_Test_12feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B0_Test_12feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B0_Train.root", // input file name
        th1d_name_arr_7, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_12feat_B0_Train_40_60_B0_Test.pdf", // plot file name
        "12feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    char* th1d_name_arr_8[5] = {
        "th1d_Tree_40_60_B4_Test_12feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B4_Test_12feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B4_Test_12feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B4_Test_12feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B4_Test_12feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B0_Train.root", // input file name
        th1d_name_arr_8, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_12feat_B0_Train_40_60_B4_Test.pdf", // plot file name
        "12feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    char* th1d_name_arr_9[5] = {
        "th1d_Tree_40_60_B8_Test_12feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B8_Test_12feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B8_Test_12feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B8_Test_12feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B8_Test_12feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B0_Train.root", // input file name
        th1d_name_arr_9, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_12feat_B0_Train_40_60_B8_Test.pdf", // plot file name
        "12feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    // --- Training with pT^4 Bias (B4)
    
    // 1 Feature Test
    
    char* th1d_name_arr_11[5] = {
        "th1d_Tree_40_60_B0_Test_1feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B0_Test_1feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B0_Test_1feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B0_Test_1feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B0_Test_1feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B4_Train.root", // input file name
        th1d_name_arr_11, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_1feat_B4_Train_40_60_B0_Test.pdf", // plot file name
        "1feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    char* th1d_name_arr_12[5] = {
        "th1d_Tree_40_60_B4_Test_1feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B4_Test_1feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B4_Test_1feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B4_Test_1feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B4_Test_1feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B4_Train.root", // input file name
        th1d_name_arr_12, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_1feat_B4_Train_40_60_B4_Test.pdf", // plot file name
        "1feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    char* th1d_name_arr_13[5] = {
        "th1d_Tree_40_60_B8_Test_1feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B8_Test_1feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B8_Test_1feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B8_Test_1feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B8_Test_1feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B4_Train.root", // input file name
        th1d_name_arr_13, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_1feat_B4_Train_40_60_B8_Test.pdf", // plot file name
        "1feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    // 3 Feature Test
    
    char* th1d_name_arr_14[5] = {
        "th1d_Tree_40_60_B0_Test_3feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B0_Test_3feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B0_Test_3feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B0_Test_3feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B0_Test_3feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B4_Train.root", // input file name
        th1d_name_arr_14, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_3feat_B4_Train_40_60_B0_Test.pdf", // plot file name
        "3feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    char* th1d_name_arr_15[5] = {
        "th1d_Tree_40_60_B4_Test_3feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B4_Test_3feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B4_Test_3feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B4_Test_3feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B4_Test_3feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B4_Train.root", // input file name
        th1d_name_arr_15, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_3feat_B4_Train_40_60_B4_Test.pdf", // plot file name
        "3feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    char* th1d_name_arr_16[5] = {
        "th1d_Tree_40_60_B8_Test_3feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B8_Test_3feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B8_Test_3feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B8_Test_3feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B8_Test_3feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B4_Train.root", // input file name
        th1d_name_arr_16, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_3feat_B4_Train_40_60_B8_Test.pdf", // plot file name
        "3feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    // 12 Feature Test
    
    char* th1d_name_arr_17[5] = {
        "th1d_Tree_40_60_B0_Test_12feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B0_Test_12feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B0_Test_12feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B0_Test_12feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B0_Test_12feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B4_Train.root", // input file name
        th1d_name_arr_17, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_12feat_B4_Train_40_60_B0_Test.pdf", // plot file name
        "12feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    char* th1d_name_arr_18[5] = {
        "th1d_Tree_40_60_B4_Test_12feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B4_Test_12feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B4_Test_12feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B4_Test_12feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B4_Test_12feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B4_Train.root", // input file name
        th1d_name_arr_18, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_12feat_B4_Train_40_60_B4_Test.pdf", // plot file name
        "12feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    char* th1d_name_arr_19[5] = {
        "th1d_Tree_40_60_B8_Test_12feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B8_Test_12feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B8_Test_12feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B8_Test_12feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B8_Test_12feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B4_Train.root", // input file name
        th1d_name_arr_19, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_12feat_B4_Train_40_60_B8_Test.pdf", // plot file name
        "12feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    // --- Training with pT^8 Bias (B48)
    
    // 1 Feature Test
    
    char* th1d_name_arr_21[5] = {
        "th1d_Tree_40_60_B0_Test_1feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B0_Test_1feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B0_Test_1feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B0_Test_1feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B0_Test_1feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B8_Train.root", // input file name
        th1d_name_arr_21, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_1feat_B8_Train_40_60_B0_Test.pdf", // plot file name
        "1feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    char* th1d_name_arr_22[5] = {
        "th1d_Tree_40_60_B4_Test_1feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B4_Test_1feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B4_Test_1feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B4_Test_1feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B4_Test_1feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B8_Train.root", // input file name
        th1d_name_arr_22, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_1feat_B8_Train_40_60_B4_Test.pdf", // plot file name
        "1feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    char* th1d_name_arr_23[5] = {
        "th1d_Tree_40_60_B8_Test_1feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B8_Test_1feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B8_Test_1feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B8_Test_1feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B8_Test_1feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B8_Train.root", // input file name
        th1d_name_arr_23, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_1feat_B8_Train_40_60_B8_Test.pdf", // plot file name
        "1feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    // 3 Feature Test
    
    char* th1d_name_arr_24[5] = {
        "th1d_Tree_40_60_B0_Test_3feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B0_Test_3feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B0_Test_3feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B0_Test_3feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B0_Test_3feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B8_Train.root", // input file name
        th1d_name_arr_24, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_3feat_B8_Train_40_60_B0_Test.pdf", // plot file name
        "3feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    char* th1d_name_arr_25[5] = {
        "th1d_Tree_40_60_B4_Test_3feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B4_Test_3feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B4_Test_3feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B4_Test_3feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B4_Test_3feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B8_Train.root", // input file name
        th1d_name_arr_25, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_3feat_B8_Train_40_60_B4_Test.pdf", // plot file name
        "3feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    char* th1d_name_arr_26[5] = {
        "th1d_Tree_40_60_B8_Test_3feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B8_Test_3feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B8_Test_3feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B8_Test_3feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B8_Test_3feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B8_Train.root", // input file name
        th1d_name_arr_26, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_3feat_B8_Train_40_60_B8_Test.pdf", // plot file name
        "3feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    // 12 Feature Test
    
    char* th1d_name_arr_27[5] = {
        "th1d_Tree_40_60_B0_Test_12feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B0_Test_12feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B0_Test_12feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B0_Test_12feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B0_Test_12feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B8_Train.root", // input file name
        th1d_name_arr_27, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_12feat_B8_Train_40_60_B0_Test.pdf", // plot file name
        "12feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    char* th1d_name_arr_28[5] = {
        "th1d_Tree_40_60_B4_Test_12feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B4_Test_12feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B4_Test_12feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B4_Test_12feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B4_Test_12feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B8_Train.root", // input file name
        th1d_name_arr_28, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_12feat_B8_Train_40_60_B4_Test.pdf", // plot file name
        "12feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
    char* th1d_name_arr_29[5] = {
        "th1d_Tree_40_60_B8_Test_12feat_ptTruePythia_feature_importance",
        "th1d_Tree_40_60_B8_Test_12feat_ptTruePythia_simple_correction",
        "th1d_Tree_40_60_B8_Test_12feat_ptTruePythia_linear_regression",
        "th1d_Tree_40_60_B8_Test_12feat_ptTruePythia_random_forest",
        "th1d_Tree_40_60_B8_Test_12feat_ptTruePythia_neural_network"};
    Jet_ML_Plotter(
        "ResultsML_Results_10_90_B8_Train.root", // input file name
        th1d_name_arr_29, // th1d_name_arr array of tree names
        "th1d_jet_pt_delta_12feat_B8_Train_40_60_B8_Test.pdf", // plot file name
        "12feature_40-60", // pt_test_bin base name
        40., 60., // pt_min, pt_max
        true, // bool_showLegend
        true, // bool_showFeatures
        true, // bool_normalize,
        true, // bool_compareToPaper
        "Pythia" // truth source
        );
    
//    Jet_ptTrue_Distribution_Plotter(
//        "ML_Results_20_80.root", // input file name
//        "th1_20_80_Train_ptTrue_distribution", // plot tree name
//        "12feature_20-80", // pt_test_bin base name
//        20., 80. // pt_min, pt_max
//        );
//
//    Jet_Truth_Plotter(
//        "ML_Results_10_90_Tree_40_60_Test_12feat_ptTruePythia.root", //  input file name
//        "th1d_10_90_Tree_40_60_Test_12feat_ptTruePythia.pdf", // plot file name
//        "12feature_ptTruePythia", // pt_test_bin base name
//        40., 60., // pt_min, pt_max
//        "Pythia" // truth source
//        );

}


