#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <vector>
#include <sstream>

#include "jet_ml_constants.h"
using namespace Jet_ML_Constants;

namespace Jet_Plotter {

void Build_Results_TTree_FromCSV(
    char input_file_path[400],  // input file name
    char output_file_path[400], // output file name
    char output_tree_name[100]  // name for TTree in output file
    ) {
    
    fstream csv_file;
    csv_file.open( input_file_path, ios::in );
    
    std::cout << "Accessed input file." << std::endl;
    
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



void Build_Coeffs_TTree_FromCSV(
    char input_file_path[400],  // input file name
    char output_file_path[400], // output file name
    char output_tree_name[100]  // name for TTree in output file
    ) {
    
    fstream csv_file;
    csv_file.open( input_file_path, ios::in );
    
    std::cout << "Accessed input file." << std::endl;
    
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



void Plot_JetPt_ML_Corr_True_Difference(
    char  input_file_path[400],     // Full file name with path for input .root file
    char  results_tree_name[200],   // Name of TTree to get ML data from
    char  plot_file_name[400],      // Full output file name path WITH NO FILE EXTENSION
    char  plot_title[200],          // Title to be printed on plot
    float test_pt_min,
    float test_pt_max,
    bool  use_delta = false,
    bool  use_normalized = true,
    float y_max = 0.2
    ) {
    
    TFile* input_file = new TFile(input_file_path, "READ");
    std::cout << "Input file accessed: " << input_file_path << std::endl;
    TTree* tree_results = (TTree*) input_file->Get(results_tree_name);
    std::cout << "Input tree accessed: " << results_tree_name << std::endl;
    
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
    th1f_canvas_plot->SetAxisRange(0., y_max, "Y");
    
    char true_name[100];
    sprintf(true_name, "th1f_jet_pt_true_test_%.0f_%.0f%s", test_pt_min, test_pt_max, plot_type);
    TH1F* th1f_jet_pt_true = new TH1F(true_name, "True Jet pT; p_{T, ch jet} [GeV]; Probability Density", x_bins, x_min, x_max);
    
    char reco_name[100];
    sprintf(reco_name, "th1f_jet_pt_reco_test_%.0f_%.0f%s", test_pt_min, test_pt_max, plot_type);
    TH1F* th1f_jet_pt_reco = new TH1F(reco_name, "Jet pT Reconstructed by ML; p_{T, ch jet} [GeV]; Probability Density", x_bins, x_min, x_max);
    
    char corr_name[100];
    sprintf(corr_name, "th1f_jet_pt_corr_test_%.0f_%.0f%s", test_pt_min, test_pt_max, plot_type);
    TH1F* tg_jet_pt_corr = new TH1F(corr_name, "Jet pT from Simple Correction; p_{T, ch jet} [GeV]; Probability Density", x_bins, x_min, x_max);
    
    th1f_jet_pt_true->Sumw2();
    th1f_jet_pt_reco->Sumw2();
    tg_jet_pt_corr->Sumw2();
    
    for ( int j=0 ; j < tree_results->GetEntries() ; j++ ) {
        tree_results->GetEntry(j);
        th1f_jet_pt_true->Fill(jet_pt_true);
        
        if (use_delta) {
            th1f_jet_pt_reco->Fill(jet_pt_reco - jet_pt_true);
            tg_jet_pt_corr->Fill(jet_pt_corr - jet_pt_true);
        }
        else {
            th1f_jet_pt_reco->Fill(jet_pt_reco);
            tg_jet_pt_corr->Fill(jet_pt_corr);
        }
        
//        std::cout << jet_pt_true << ", " << jet_pt_reco << ", " << jet_pt_corr << std::endl;
    }
    
    // Uncomment to plot probability density
    if (!use_delta) th1f_jet_pt_true->Scale( 0.5 / th1f_jet_pt_true->Integral(),"WIDTH");
    th1f_jet_pt_reco->Scale( 1. / th1f_jet_pt_reco->Integral(),"WIDTH");
    tg_jet_pt_corr->Scale( 1. / tg_jet_pt_corr->Integral(),"WIDTH");
    
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
    th1f_jet_pt_true->SetLineColor(plot_black);
    th1f_jet_pt_true->SetMarkerColor(plot_black);
    th1f_jet_pt_true->SetMarkerStyle(mark_circ_open[0]);
    th1f_jet_pt_true->SetMarkerSize(mark_circ_open[1]);
    
    th1f_jet_pt_reco->SetLineColor(plot_blue);
    th1f_jet_pt_reco->SetMarkerColor(plot_blue);
    th1f_jet_pt_reco->SetMarkerStyle(mark_star_open[0]);
    th1f_jet_pt_reco->SetMarkerSize(mark_star_open[1]);
    
    tg_jet_pt_corr->SetLineColor(plot_red);
    tg_jet_pt_corr->SetMarkerColor(plot_red);
    tg_jet_pt_corr->SetMarkerStyle(mark_squa_open[0]);
    tg_jet_pt_corr->SetMarkerSize(mark_squa_open[1]);
    
    TLegend* legend = new TLegend(0.15, 0.65, 0.35, 0.85);
    char reco_stats[100];
    float reco_mean   = th1f_jet_pt_reco->GetMean();
    float reco_stddev = th1f_jet_pt_reco->GetStdDev();
    sprintf(reco_stats, "Mean: %.2f, StdDev: %.2f", reco_mean, reco_stddev);
    
    char corr_stats[100];
    float corr_mean   = tg_jet_pt_corr->GetMean();
    float corr_stddev = tg_jet_pt_corr->GetStdDev();
    sprintf(corr_stats, "Mean: %.2f, StdDev: %.2f", corr_mean, corr_stddev);
    
//    if (!use_delta) legend->AddEntry(th1f_jet_pt_true,  "p_{T}^{True} (0.5 scale)", "lp");
    legend->AddEntry(th1f_jet_pt_reco,  "Linear Regression", "lp");
    legend->AddEntry((TObject*)0, reco_stats, "");
    legend->AddEntry(tg_jet_pt_corr,  "Area Correction", "lp");
    legend->AddEntry((TObject*)0, corr_stats, "");
    legend->SetLineWidth(0);
    legend->SetFillStyle(0);
    
    th1f_canvas_plot->Draw();
    
//    if (!use_delta) th1f_jet_pt_true->Draw("same");
    tg_jet_pt_corr->Draw("same");
    th1f_jet_pt_reco->Draw("same");
    
    legend->Draw("same");
    
    canvas->Print(plot_file_name_ext);
    std::cout << "Plotted " << plot_file_name_ext << std::endl;
    
    delete legend;
    delete latex;
    delete canvas;
    delete th1f_jet_pt_true;
    delete tg_jet_pt_corr;
    delete th1f_jet_pt_reco;
    delete input_file;
    std::cout << "Made it to the end" << std::endl;
}



void Plot_LR_Coefficients(
    char  input_file_path[500],
    char  input_tree_name[10][100],
    char  plot_file_name[200],
    float test_min_max_array[20][2],
    int   bin_count,
    bool  show_intercept
    ) {
    
    TFile* input_file = new TFile(input_file_path, "READ");
    
    std::cout << "Creating graphs..." << std::endl;
    
    TMultiGraph* tmg_canvas_plot      = new TMultiGraph;
    char canvas_plot_title [200];
    snprintf(canvas_plot_title, 200, "Linear Regression Coefficients, Test/Train with Bins of %.0f GeV between 10-90 GeV; Average Train/Test p_{T} [GeV]; Coefficient Value", (test_min_max_array[0][1] - test_min_max_array[0][0]) );
    tmg_canvas_plot->SetTitle(canvas_plot_title);
    
    std::cout << "Histograms made. Filling histograms..." << std::endl;
    
    float arr_test_pt_avg[20];
    float arr_jet_pt_raw[20];
//    float arr_jet_pt_corr[20];
    float arr_jet_mass[20];
    float arr_jet_area[20];
    float arr_jet_const_n[20];
    float arr_const_pt_mean[20];
    float arr_const_pt_1[20];
    float arr_const_pt_2[20];
    float arr_const_pt_3[20];
    float arr_const_pt_4[20];
    float arr_jet_y[20];
    float arr_jet_rho[20];
    float arr_lr_intercept[20];
    
    for ( int i = 0 ; i < bin_count ; i++ ) {
        float test_pt_min = test_min_max_array[i][0];
        float test_pt_max = test_min_max_array[i][1];
        float test_pt_avg = 0.5 * (test_pt_min + test_pt_max);
        
        std::cout << "Adding data for range " << test_pt_min << " to " << test_pt_max << std::endl;
        
        TTree* coeff_tree = (TTree*) input_file->Get(input_tree_name[i]);
        std::cout << "TTree: " << input_tree_name[i] << std::endl;
        
        float jet_pt_raw;
//        float jet_pt_corr;
        float jet_mass;
        float jet_area;
        float jet_const_n;
        float const_pt_mean;
        float const_pt_1;
        float const_pt_2;
        float const_pt_3;
        float const_pt_4;
        float jet_y;
        float jet_rho;
        float lr_intercept;
        
        std::cout << "Made variables..." << std::endl;
        
        coeff_tree->SetBranchAddress("jet_pt_raw",    &jet_pt_raw);
        std::cout << "Made it here!" << std::endl;
//        coeff_tree->SetBranchAddress("jet_pt_corr",   &jet_pt_corr);
        coeff_tree->SetBranchAddress("jet_mass",      &jet_mass);
        coeff_tree->SetBranchAddress("jet_area",      &jet_area);
        coeff_tree->SetBranchAddress("jet_const_n",   &jet_const_n);
        coeff_tree->SetBranchAddress("const_pt_mean", &const_pt_mean);
        coeff_tree->SetBranchAddress("const_1_pt",    &const_pt_1);
        coeff_tree->SetBranchAddress("const_2_pt",    &const_pt_2);
        coeff_tree->SetBranchAddress("const_3_pt",    &const_pt_3);
        coeff_tree->SetBranchAddress("const_4_pt",    &const_pt_4);
        coeff_tree->SetBranchAddress("jet_y",         &jet_y);
        coeff_tree->SetBranchAddress("jet_rho",       &jet_rho);
        coeff_tree->SetBranchAddress("lr_intercept",  &lr_intercept);
        
        std::cout << "Set branches..." << std::endl;
        
        coeff_tree->GetEntry(0);
        
        arr_test_pt_avg[i]    = test_pt_avg;
        arr_jet_pt_raw[i]     = jet_pt_raw;
//        arr_jet_pt_corr[i]    = jet_pt_corr;
        arr_jet_mass[i]       = jet_mass;
        arr_jet_area[i]       = jet_area;
        arr_jet_const_n[i]    = jet_const_n;
        arr_const_pt_mean[i]  = const_pt_mean;
        arr_const_pt_1[i]     = const_pt_1;
        arr_const_pt_2[i]     = const_pt_2;
        arr_const_pt_3[i]     = const_pt_3;
        arr_const_pt_4[i]     = const_pt_4;
        arr_jet_y[i]          = jet_y;
        arr_jet_rho[i]        = jet_rho;
        arr_lr_intercept[i]   = lr_intercept;
        
        std::cout << "Filled branches..." << std::endl;
        
        delete coeff_tree;
    }
    
    TGraph* tg_jet_pt_raw       = new TGraph(bin_count, arr_test_pt_avg, arr_jet_pt_raw);
    tg_jet_pt_raw->SetTitle("Jet p_{T}^{Raw}; Average Train/Test p_{T} [GeV]; Coefficient Value");
//    TGraph* tg_jet_pt_corr      = new TGraph(bin_count, arr_test_pt_avg, arr_jet_pt_corr);
//    tg_jet_pt_corr->SetTitle("Jet p_{T}^{Corrected}; Average Train/Test p_{T} [GeV]; Coefficient Value");
    TGraph* tg_jet_mass         = new TGraph(bin_count, arr_test_pt_avg, arr_jet_mass);
    tg_jet_mass->SetTitle("Jet Mass; Average Train/Test p_{T} [GeV]; Coefficient Value");
    TGraph* tg_jet_area         = new TGraph(bin_count, arr_test_pt_avg, arr_jet_area);
    tg_jet_area->SetTitle("Jet Area; Average Train/Test p_{T} [GeV]; Coefficient Value");
    TGraph* tg_jet_const_n      = new TGraph(bin_count, arr_test_pt_avg, arr_jet_const_n);
    tg_jet_const_n->SetTitle("Number of Constituents in Jet; Average Train/Test p_{T} [GeV]; Coefficient Value");
    TGraph* tg_const_pt_mean    = new TGraph(bin_count, arr_test_pt_avg, arr_const_pt_mean);
    tg_const_pt_mean->SetTitle("Jet Constituent Mean p_{T}; Average Train/Test p_{T} [GeV]; Coefficient Value");
    TGraph* tg_const_pt_1       = new TGraph(bin_count, arr_test_pt_avg, arr_const_pt_1);
    tg_const_pt_1->SetTitle("p_{T}^{Const. 1}; Average Train/Test p_{T} [GeV]; Coefficient Value");
    TGraph* tg_const_pt_2       = new TGraph(bin_count, arr_test_pt_avg, arr_const_pt_2);
    tg_const_pt_2->SetTitle("p_{T}^{Const. 2}; Average Train/Test p_{T} [GeV]; Coefficient Value");
    TGraph* tg_const_pt_3       = new TGraph(bin_count, arr_test_pt_avg, arr_const_pt_3);
    tg_const_pt_3->SetTitle("p_{T}^{Const. 3}; Average Train/Test p_{T} [GeV]; Coefficient Value");
    TGraph* tg_const_pt_4       = new TGraph(bin_count, arr_test_pt_avg, arr_const_pt_4);
    tg_const_pt_4->SetTitle("p_{T}^{Const. 4}; Average Train/Test p_{T} [GeV]; Coefficient Value");
    TGraph* tg_jet_y            = new TGraph(bin_count, arr_test_pt_avg, arr_jet_y);
    tg_jet_y->SetTitle("Jet Rapidity, y; Average Train/Test p_{T} [GeV]; Coefficient Value");
    TGraph* tg_jet_rho          = new TGraph(bin_count, arr_test_pt_avg, arr_jet_rho);
    tg_jet_rho->SetTitle("Jet $rho; Average Train/Test p_{T} [GeV]; Coefficient Value");
    TGraph* tg_lr_intercept     = new TGraph(bin_count, arr_test_pt_avg, arr_lr_intercept);
    tg_lr_intercept->SetTitle("Linear Regression y-intercept; Average Train/Test p_{T} [GeV]; Coefficient Value");
    
    std::cout << "Histograms filled. Making plot..." << std::endl;
    
    // TCanvas("name", "title", width (px), height (px))
    TCanvas* canvas = new TCanvas("canvas", "", 1000, 600);

    // Turns on ticks, turns off stats
    gPad->SetTicks();
    gStyle->SetOptStat(0);

    // Turns on Latex formatting relative to canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(kTRUE);
    
    // Style histograms
    
    tg_jet_pt_raw->SetLineColor(kViolet+5);
    tg_jet_pt_raw->SetMarkerColor(kViolet+5);
    tg_jet_pt_raw->SetMarkerStyle(mark_circ_open[0]);
    tg_jet_pt_raw->SetMarkerSize(mark_circ_open[1]);
    
//    tg_jet_pt_corr->SetLineColor(kViolet-5);
//    tg_jet_pt_corr->SetMarkerColor(kViolet-5);
//    tg_jet_pt_corr->SetMarkerStyle(mark_squa_open[0]);
//    tg_jet_pt_corr->SetMarkerSize(mark_squa_open[1]);
    
    tg_jet_mass->SetLineColor(kPink+5);
    tg_jet_mass->SetMarkerColor(kPink+5);
    tg_jet_mass->SetMarkerStyle(mark_circ_open[0]);
    tg_jet_mass->SetMarkerSize(mark_circ_open[1]);
    
    tg_jet_area->SetLineColor(kPink-5);
    tg_jet_area->SetMarkerColor(kPink-5);
    tg_jet_area->SetMarkerStyle(mark_squa_open[0]);
    tg_jet_area->SetMarkerSize(mark_squa_open[1]);
    
    tg_jet_const_n->SetLineColor(kOrange+5);
    tg_jet_const_n->SetMarkerColor(kOrange+5);
    tg_jet_const_n->SetMarkerStyle(mark_circ_open[0]);
    tg_jet_const_n->SetMarkerSize(mark_circ_open[1]);
    
    tg_const_pt_mean->SetLineColor(kOrange-5);
    tg_const_pt_mean->SetMarkerColor(kOrange-5);
    tg_const_pt_mean->SetMarkerStyle(mark_squa_open[0]);
    tg_const_pt_mean->SetMarkerSize(mark_squa_open[1]);
    
    tg_const_pt_1->SetLineColor(kSpring+5);
    tg_const_pt_1->SetMarkerColor(kSpring+5);
    tg_const_pt_1->SetMarkerStyle(mark_circ_open[0]);
    tg_const_pt_1->SetMarkerSize(mark_circ_open[1]);
    
    tg_const_pt_2->SetLineColor(kSpring-5);
    tg_const_pt_2->SetMarkerColor(kSpring-5);
    tg_const_pt_2->SetMarkerStyle(mark_squa_open[0]);
    tg_const_pt_2->SetMarkerSize(mark_squa_open[1]);
    
    tg_const_pt_3->SetLineColor(kTeal+5);
    tg_const_pt_3->SetMarkerColor(kTeal+5);
    tg_const_pt_3->SetMarkerStyle(mark_circ_open[0]);
    tg_const_pt_3->SetMarkerSize(mark_circ_open[1]);
    
    tg_const_pt_4->SetLineColor(kTeal-5);
    tg_const_pt_4->SetMarkerColor(kTeal-5);
    tg_const_pt_4->SetMarkerStyle(mark_squa_open[0]);
    tg_const_pt_4->SetMarkerSize(mark_squa_open[1]);
    
    tg_jet_y->SetLineColor(kAzure+5);
    tg_jet_y->SetMarkerColor(kAzure+5);
    tg_jet_y->SetMarkerStyle(mark_circ_open[0]);
    tg_jet_y->SetMarkerSize(mark_circ_open[1]);
    
    tg_jet_rho->SetLineColor(kAzure-5);
    tg_jet_rho->SetMarkerColor(kAzure-5);
    tg_jet_rho->SetMarkerStyle(mark_squa_open[0]);
    tg_jet_rho->SetMarkerSize(mark_squa_open[1]);
    
    tg_lr_intercept->SetLineColor(kGray+3);
    tg_lr_intercept->SetMarkerColor(kGray+3);
    tg_lr_intercept->SetMarkerStyle(mark_circ_open[0]);
    tg_lr_intercept->SetMarkerSize(mark_circ_open[1]);
    
    // Plot coefficient histograms onto canvas
    if ( show_intercept ) tmg_canvas_plot->Add(tg_lr_intercept, "PL"); // Draws intercept values first to be in background
    tmg_canvas_plot->Add(tg_jet_pt_raw, "PL");
//    tmg_canvas_plot->Add(tg_jet_pt_corr, "PL");
    tmg_canvas_plot->Add(tg_jet_mass, "PL");
    tmg_canvas_plot->Add(tg_jet_area, "PL");
    tmg_canvas_plot->Add(tg_jet_const_n, "PL");
    tmg_canvas_plot->Add(tg_const_pt_mean, "PL");
    tmg_canvas_plot->Add(tg_const_pt_1, "PL");
    tmg_canvas_plot->Add(tg_const_pt_2, "PL");
    tmg_canvas_plot->Add(tg_const_pt_3, "PL");
    tmg_canvas_plot->Add(tg_const_pt_4, "PL");
    tmg_canvas_plot->Add(tg_jet_y, "PL");
    tmg_canvas_plot->Add(tg_jet_rho, "PL");
    tmg_canvas_plot->Draw("A");
    
    // Make and plot legend
    std::cout << "Making legend..." << std::endl;
    TLegend* legend = new TLegend(0.65, 0.50, 0.85, 0.85);
    
    legend->AddEntry(tg_jet_pt_raw,       "Jet p_{T}^{Raw}", "lp");
//    legend->AddEntry(tg_jet_pt_corr,      "Jet p_{T}^{Area Correction}", "lp");
    legend->AddEntry(tg_jet_mass,         "Jet Mass", "lp");
    legend->AddEntry(tg_jet_area,         "Jet Area", "lp");
    legend->AddEntry(tg_jet_const_n,      "Number of Constituents", "lp");
    legend->AddEntry(tg_const_pt_mean,    "Mean Constituent p_{T}", "lp");
    legend->AddEntry(tg_const_pt_1,       "Constituent 1 p_{T}", "lp");
    legend->AddEntry(tg_const_pt_2,       "Constituent 2 p_{T}", "lp");
    legend->AddEntry(tg_const_pt_3,       "Constituent 3 p_{T}", "lp");
    legend->AddEntry(tg_const_pt_4,       "Constituent 4 p_{T}", "lp");
    legend->AddEntry(tg_jet_y,            "Jet y", "lp");
    legend->AddEntry(tg_jet_rho,          "Jet #rho", "lp");
    if ( show_intercept ) legend->AddEntry(tg_lr_intercept,     "Lin. Reg. Intercept", "lp");
    
    legend->SetLineWidth(0);
    legend->SetFillStyle(0);
    
    legend->Draw("same");
    
    canvas->Print(plot_file_name);
    std::cout << "Plotted " << plot_file_name << std::endl;
}
    
    
double Gaussian_Func(double *x,double *par) {
    double arg = 0;
    if (par[2]!=0) arg = (x[0] - par[1])/par[2];
    double fitval = par[0]*TMath::Exp(-0.5*arg*arg);
    return fitval;
}

void Fit_Paper_Data(
    TF1* tf1_paper_simple_correction_fit,
    TF1* tf1_paper_neural_network_fit,
    TF1* tf1_paper_random_forest_fit,
    TF1* tf1_paper_linear_regression_fit
    ) {
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

void Plot_JetPt_ML_Comparison(
    char  input_file_path[400],
    char  plot_file_dir[400],
    char  plot_title[100],
    char  lr_tree_name[100],
    char  rf_tree_name[100],
    char  mlp_tree_name[100],
    char  label_ptbias[100],
    char  label_feature[100],
    char  coeff_tree_name[100],
    char  coeff_label[100],
    float train_pt_min,
    float train_pt_max,
    float test_pt_min,
    float test_pt_max,
    bool  normalize = true,
    bool  show_legend = true,
    bool  show_widths = true,
    bool  show_coeffs = false,
    bool  show_paper_plots = false,
    bool  use_delta = true,
    float y_max = 0.2
    ) {
    
    TFile* input_file = new TFile(input_file_path, "READ");
    TTree* lr_tree      = (TTree*) input_file->Get(lr_tree_name);
    TTree* rf_tree      = (TTree*) input_file->Get(rf_tree_name);
    TTree* mlp_tree     = (TTree*) input_file->Get(mlp_tree_name);
    TTree* coeff_tree   = (TTree*) input_file->Get(coeff_tree_name);
    
    TObjArray* coeff_branch_list = coeff_tree->GetListOfBranches();
    vector<string>  branch_label_array;
    vector<float>   branch_coeff_array;
    
//    for ( int i=0 ; i < coeff_branch_list->IndexOf(coeff_branch_list->Last()) ; i++ ) {
//        branch_label_array[i] = (string) coeff_branch_list[i];
//        float jet_feature;
//        coeff_tree->SetBranchAddress(branch_label_array[i].c_str(), &jet_feature);
//        coeff_tree->GetEntry(0);
//        branch_coeff_array[i] = jet_feature;
//    }
    float x_min = test_pt_min - 20.;
    float x_max = test_pt_max + 20.;
    if ( use_delta ) {
        x_min = -30.; // (test_pt_min-test_pt_max)*0.5 - 20.;
        x_max =  30.; // (test_pt_max-test_pt_min)*0.5 + 20.;
    }
    int x_bins = int((x_max - x_min) / 1);
    
    char canvas_title_units[200];
    string x_label = "p_{T} [GeV]";
    string y_label = "N_{Jets}";
    if ( use_delta ) x_label = "p_{T}^{Reco} - p_{T}^{True} [GeV]";
    if ( normalize ) y_label = "Probability Density";
    snprintf(canvas_title_units, 200, "%s;%s;%s", plot_title, x_label.c_str(), y_label.c_str());
    TH1F* th1_canvas_plot       = new TH1F("canvas_plot", canvas_title_units, x_bins, x_min, x_max);
    th1_canvas_plot->SetAxisRange(0, y_max, "Y");
    TH1F* th1_simple_correction = new TH1F("sc_plot", "", x_bins, x_min, x_max);
    TH1F* th1_linear_regression = new TH1F("lr_plot", "", x_bins, x_min, x_max);
    TH1F* th1_random_forest     = new TH1F("rf_plot", "", x_bins, x_min, x_max);
    TH1F* th1_neural_network    = new TH1F("mlp_plot", "", x_bins, x_min, x_max);
    
    float jet_pt_true;
    float jet_pt_corr;
    float jet_pt_reco_lr;
    float jet_pt_reco_rf;
    float jet_pt_reco_mlp;
    
    lr_tree->SetBranchAddress("jet_pt_true", &jet_pt_true);
    lr_tree->SetBranchAddress("jet_pt_corr", &jet_pt_corr);
    lr_tree->SetBranchAddress("jet_pt_reco", &jet_pt_reco_lr);
    rf_tree->SetBranchAddress("jet_pt_reco", &jet_pt_reco_rf);
    mlp_tree->SetBranchAddress("jet_pt_reco", &jet_pt_reco_mlp);
    
    if ( use_delta ) {
        for ( int e=0 ; e < lr_tree->GetEntries() ; e++ ) {
            lr_tree->GetEvent(e);
            rf_tree->GetEvent(e);
            mlp_tree->GetEvent(e);
            th1_simple_correction->Fill(jet_pt_corr - jet_pt_true);
            th1_linear_regression->Fill(jet_pt_reco_lr - jet_pt_true);
            th1_random_forest->Fill(jet_pt_reco_rf - jet_pt_true);
            th1_neural_network->Fill(jet_pt_reco_mlp - jet_pt_true);
        }
    }
    else {
        for ( int e=0 ; e < lr_tree->GetEntries() ; e++ ) {
            lr_tree->GetEvent(e);
            rf_tree->GetEvent(e);
            mlp_tree->GetEvent(e);
            th1_simple_correction->Fill(jet_pt_corr);
            th1_linear_regression->Fill(jet_pt_reco_lr);
            th1_random_forest->Fill(jet_pt_reco_rf);
            th1_neural_network->Fill(jet_pt_reco_mlp);
        }
    }
    std::cout << "Histograms filled." << std::endl;
    
    if ( normalize || show_paper_plots ) {
        th1_simple_correction  ->Scale( 1. / th1_simple_correction->Integral(),"WIDTH");
        th1_linear_regression  ->Scale( 1. / th1_linear_regression->Integral(),"WIDTH");
        th1_random_forest      ->Scale( 1. / th1_random_forest->Integral(),"WIDTH");
        th1_neural_network     ->Scale( 1. / th1_neural_network->Integral(),"WIDTH");
    }
    std::cout << "Histograms normalized." << std::endl;
    
    // Apply fit to each plot
    TF1* tf1_simple_correction_fit  = new TF1("tf1_simple_correction_fit",  "gaus", x_min, x_max);
    TF1* tf1_linear_regression_fit  = new TF1("tf1_linear_regression_fit",  "gaus", x_min, x_max);
    TF1* tf1_random_forest_fit      = new TF1("tf1_random_forest_fit",      "gaus", x_min, x_max);
    TF1* tf1_neural_network_fit     = new TF1("tf1_neural_network_fit",     "gaus", x_min, x_max);

    th1_simple_correction   ->Fit(tf1_simple_correction_fit,"RNL");
    th1_linear_regression   ->Fit(tf1_linear_regression_fit,"RNL");
    th1_random_forest       ->Fit(tf1_random_forest_fit,"RNL");
    th1_neural_network      ->Fit(tf1_neural_network_fit,"RNL");
    
    // TCanvas("name", "title", width (px), height (px))
    TCanvas* canvas = new TCanvas("canvas", "", 1000, 600);

    // Turns on ticks, turns off stats
    gPad->SetTicks();
    gStyle->SetOptStat(0);

    // Turns on Latex formatting relative to canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(kTRUE);
    
    // Build histogram plot
    th1_simple_correction->SetLineColor(plot_red);
    th1_simple_correction->SetMarkerColor(plot_red);
    th1_simple_correction->SetMarkerStyle(mark_circ_open[0]);
    th1_simple_correction->SetMarkerSize(mark_circ_open[1]);
    tf1_simple_correction_fit->SetLineColor(plot_red);
    tf1_simple_correction_fit->SetLineWidth(1);
    
    th1_linear_regression->SetLineColor(plot_blue);
    th1_linear_regression->SetMarkerColor(plot_blue);
    th1_linear_regression->SetMarkerStyle(mark_squa_open[0]);
    th1_linear_regression->SetMarkerSize(mark_squa_open[1]);
    tf1_linear_regression_fit->SetLineColor(plot_blue);
    tf1_linear_regression_fit->SetLineWidth(1);
    
    th1_random_forest->SetLineColor(plot_green);
    th1_random_forest->SetMarkerColor(plot_green);
    th1_random_forest->SetMarkerStyle(mark_diam_open[0]);
    th1_random_forest->SetMarkerSize(mark_diam_open[1]);
    tf1_random_forest_fit->SetLineColor(plot_green);
    tf1_random_forest_fit->SetLineWidth(1);

    th1_neural_network->SetLineColor(plot_black);
    th1_neural_network->SetMarkerColor(plot_black);
    th1_neural_network->SetMarkerStyle(mark_star_open[0]);
    th1_neural_network->SetMarkerSize(mark_star_open[1]);
    tf1_neural_network_fit->SetLineColor(plot_black);
    tf1_neural_network_fit->SetLineWidth(1);
    
    // Draws the histograms and lines
    th1_canvas_plot->Draw();
    
    // Draws results from ML paper
    TF1* tf1_paper_simple_correction_fit;
    TF1* tf1_paper_neural_network_fit;
    TF1* tf1_paper_random_forest_fit;
    TF1* tf1_paper_linear_regression_fit;
    if ( show_paper_plots ) {
        Fit_Paper_Data(
            tf1_paper_simple_correction_fit,
            tf1_paper_neural_network_fit,
            tf1_paper_random_forest_fit,
            tf1_paper_linear_regression_fit
        );
        
        tf1_paper_neural_network_fit->SetLineColor(plot_black);
        tf1_paper_neural_network_fit->SetLineWidth(1);
        tf1_paper_neural_network_fit->SetLineStyle(3);
    
        tf1_paper_random_forest_fit->SetLineColor(plot_green);
        tf1_paper_random_forest_fit->SetLineWidth(1);
        tf1_paper_random_forest_fit->SetLineStyle(3);
    
        tf1_paper_linear_regression_fit->SetLineColor(plot_blue);
        tf1_paper_linear_regression_fit->SetLineWidth(1);
        tf1_paper_linear_regression_fit->SetLineStyle(3);
    
        tf1_paper_simple_correction_fit->SetLineColor(plot_red);
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
    th1_simple_correction       ->Draw("same");
    th1_neural_network          ->Draw("same");
    th1_random_forest           ->Draw("same");
    th1_linear_regression       ->Draw("same");
    
    string legend_title_line1_str;
    string legend_title_line2_str = "Train " + std::to_string(int(train_pt_min)) + "-" + std::to_string(int(train_pt_max)) + " GeV, Test " + std::to_string(int(test_pt_min)) + "-" + std::to_string(int(test_pt_max)) + " GeV";
    string title_bias;
    if      ( !string(label_ptbias).compare("B0") ) title_bias = "No Bias";
    else if ( !string(label_ptbias).compare("B4") ) title_bias = "p_{T}^{4} Bias";
    else if ( !string(label_ptbias).compare("B8") ) title_bias = "p_{T}^{8} Bias";
    else if ( !string(label_ptbias).compare("B8_Flat") ) title_bias = "Flat";
    else if ( !string(label_ptbias).compare("B8_Flat_T2") ) title_bias = "Flat";
    else title_bias = "Unknown";
    string title_feature;
    if      ( !string(label_feature).compare("F1") ) title_feature = "1 Feature";
    else if ( !string(label_feature).compare("F3") ) title_feature = "3 Features";
    else if ( !string(label_feature).compare("F11") ) title_feature = "11 Features";
    else if ( !string(label_feature).compare("F12") ) title_feature = "12 Features";
    else title_feature = "Unknown";
    legend_title_line1_str = title_bias + " Distribution, " + title_feature;
    string legend_title_all = legend_title_line1_str + "\n" + legend_title_line2_str;
    
    
    char  legend_title_1[100];
    snprintf(legend_title_1, 100, "%s", legend_title_line1_str.c_str());
    char  legend_title_2[100];
    snprintf(legend_title_2, 100, "%s", legend_title_line2_str.c_str());
    char  nn_stats[100];
    float nn_mean   = th1_neural_network->GetMean();
    float nn_stddev = th1_neural_network->GetStdDev();
    snprintf(nn_stats, 100, "#Deltap_{T}: %.2f, #sigma: %.2f", nn_mean, nn_stddev);
    char  rf_stats[100];
    float rf_mean   = th1_random_forest->GetMean();
    float rf_stddev = th1_random_forest->GetStdDev();
    snprintf(rf_stats, 100, "#Deltap_{T}: %.2f, #sigma: %.2f", rf_mean, rf_stddev);
    char  lr_stats[100];
    float lr_mean   = th1_linear_regression->GetMean();
    float lr_stddev = th1_linear_regression->GetStdDev();
    snprintf(lr_stats, 100, "#Deltap_{T}: %.2f, #sigma: %.2f", lr_mean, lr_stddev);
    char  sc_stats[100];
    float sc_mean   = th1_simple_correction->GetMean();
    float sc_stddev = th1_simple_correction->GetStdDev();
    snprintf(sc_stats, 100, "#Deltap_{T}: %.2f, #sigma: %.2f", sc_mean, sc_stddev);
    
    // Draws Legend as TLegend(xmin, ymin, xmax, ymax)
    if ( show_legend ) {
        TLegend* legend;
        if      ( show_widths and show_paper_plots ) legend = new TLegend(0.60, 0.25, 0.85, 0.85);
        else if ( show_widths xor show_paper_plots ) legend = new TLegend(0.60, 0.45, 0.85, 0.85);
        else legend = new TLegend(0.60, 0.65, 0.85, 0.85);
//        legend->AddEntry((TObject*)0, legend_title_1, ""); // Adds a title to the legend
//        legend->AddEntry((TObject*)0, legend_title_2, ""); // Adds a title to the legend
        legend->AddEntry(th1_neural_network,       "Multilayer Perceptron", "pl");
        if ( show_widths ) legend->AddEntry((TObject*)0, nn_stats, "");
        legend->AddEntry(th1_random_forest,        "Random Forest", "pl");
        if ( show_widths ) legend->AddEntry((TObject*)0, rf_stats, "");
        legend->AddEntry(th1_linear_regression,    "Linear Regression", "pl");
        if ( show_widths ) legend->AddEntry((TObject*)0, lr_stats, "");
        legend->AddEntry(th1_simple_correction,    "Area Correction", "pl");
        if ( show_widths ) legend->AddEntry((TObject*)0, sc_stats, "");
        if ( show_paper_plots) {
            legend->AddEntry(tf1_paper_neural_network_fit,      "NN (Paper)", "l");
            legend->AddEntry(tf1_paper_random_forest_fit,       "RF (Paper)", "l");
            legend->AddEntry(tf1_paper_linear_regression_fit,   "LR (Paper)", "l");
            legend->AddEntry(tf1_paper_simple_correction_fit,   "AC (Paper)", "l");
        }
        legend->SetLineWidth(0);
        legend->SetFillStyle(0);
        legend->Draw();
        
        latex->DrawLatex(0.15, 0.81, ("#scale[0.75]{" + legend_title_line1_str + "}").c_str());
        latex->DrawLatex(0.15, 0.75, ("#scale[0.75]{" + legend_title_line2_str + "}").c_str());
    }
    
    // Draws Feature List
//    if ( show_widths ) {
//        char sigma_list[5][100];
//        sprintf(sigma_list[0], "#scale[0.6]{Fit Widths}");
//        sprintf(sigma_list[1], "#scale[0.6]{#bf{#sigma_{Neural Network}}}");
//        sprintf(sigma_list[2], "#scale[0.6]{#bf{#sigma_{Random Forest}}}");
//        sprintf(sigma_list[3], "#scale[0.6]{#bf{#sigma_{Linear Regression}}}");
//        sprintf(sigma_list[4], "#scale[0.6]{#bf{#sigma_{Area Correction}}}");
//        char sigma_vals[5][100];
//        sprintf(sigma_vals[0], "");
//        sprintf(sigma_vals[1], "#scale[0.6]{#bf{%1.3f}}", tf1_neural_network_fit    ->GetParameter(2));
//        sprintf(sigma_vals[2], "#scale[0.6]{#bf{%1.3f}}", tf1_random_forest_fit     ->GetParameter(2));
//        sprintf(sigma_vals[3], "#scale[0.6]{#bf{%1.3f}}", tf1_linear_regression_fit ->GetParameter(2));
//        sprintf(sigma_vals[4], "#scale[0.6]{#bf{%1.3f}}", tf1_simple_correction_fit ->GetParameter(2));
//        for ( int i = 0 ; i < 5 ; i++ ) latex->DrawLatex(0.63, (.813 - 0.0494 * i), sigma_list[i]);
//        for ( int i = 0 ; i < 5 ; i++ ) latex->DrawLatex(0.80, (.813 - 0.0494 * i), sigma_vals[i]);
//    }
    
    // Prints out the plot
    char plot_output[500];
    snprintf(plot_output, 500, "%s%s.pdf", plot_file_dir, plot_title);
    canvas->Print(plot_output);
    std::cout << "Plotted " << plot_output << std::endl;
    gPad->SetLogy(0);
}

//
//void Plot_ML_pT_Comparison(
//    char* plot_filename,
//    char* plot_directory,
//    char* plot_title_xlabel_ylabel,
//    TH1D* th1d_simple_correction,
//    TH1D* th1d_linear_regression,
//    TH1D* th1d_random_forest,
//    TH1D* th1d_neural_network,
//    TH1D* th1d_feature_list,
//    double y_max,
//    bool bool_showLegend,
//    char* showFeatures,
//    bool bool_compareToPaper
//    ) {
//
//    char feature_list_1[2][100];
//    sprintf(feature_list_1[0], "#scale[0.6]{Feature Importance}");
//    sprintf(feature_list_1[1], "#scale[0.6]{#bf{Jet p_{T, raw}:}}");
//
//    char feature_list_3[4][100];
//    sprintf(feature_list_3[0], "#scale[0.6]{Feature Importance}");
//    sprintf(feature_list_3[1], "#scale[0.6]{#bf{Jet p_{T, raw}:}}");
//    sprintf(feature_list_3[2], "#scale[0.6]{#bf{Jet area:}}");
//    sprintf(feature_list_3[3], "#scale[0.6]{#bf{Jet #rho:}}");
//
//    char feature_list_11[12][100];
//    sprintf(feature_list_11[0], "#scale[0.6]{Feature Importance}");
//    sprintf(feature_list_11[1], "#scale[0.6]{#bf{Jet p_{T, raw}:}}");
//    sprintf(feature_list_11[2], "#scale[0.6]{#bf{Jet Mass:}}");
//    sprintf(feature_list_11[3], "#scale[0.6]{#bf{Jet Area:}}");
//    sprintf(feature_list_11[4], "#scale[0.6]{#bf{N_{const}:}}");
//    sprintf(feature_list_11[5], "#scale[0.6]{#bf{Mean Const. p_{T}:}}");
//    sprintf(feature_list_11[6], "#scale[0.6]{#bf{p_{T, Const.}^{1}:}}");
//    sprintf(feature_list_11[7], "#scale[0.6]{#bf{p_{T, Const.}^{2}:}}");
//    sprintf(feature_list_11[8], "#scale[0.6]{#bf{p_{T, Const.}^{3}:}}");
//    sprintf(feature_list_11[9], "#scale[0.6]{#bf{p_{T, Const.}^{4}:}}");
//    sprintf(feature_list_11[10], "#scale[0.6]{#bf{Jet y:}}");
//    sprintf(feature_list_11[11], "#scale[0.6]{#bf{Jet #rho:}}");
//
//    char feature_list_12[13][100];
//    sprintf(feature_list_12[0], "#scale[0.6]{Feature Importance}");
//    sprintf(feature_list_12[1], "#scale[0.6]{#bf{Jet p_{T, raw}:}}");
//    sprintf(feature_list_12[2], "#scale[0.6]{#bf{Jet p_{T, corr}:}}");
//    sprintf(feature_list_12[3], "#scale[0.6]{#bf{Jet Mass:}}");
//    sprintf(feature_list_12[4], "#scale[0.6]{#bf{Jet Area:}}");
//    sprintf(feature_list_12[5], "#scale[0.6]{#bf{N_{const}:}}");
//    sprintf(feature_list_12[6], "#scale[0.6]{#bf{Mean Const. p_{T}:}}");
//    sprintf(feature_list_12[7], "#scale[0.6]{#bf{p_{T, Const.}^{1}:}}");
//    sprintf(feature_list_12[8], "#scale[0.6]{#bf{p_{T, Const.}^{2}:}}");
//    sprintf(feature_list_12[9], "#scale[0.6]{#bf{p_{T, Const.}^{3}:}}");
//    sprintf(feature_list_12[10], "#scale[0.6]{#bf{p_{T, Const.}^{4}:}}");
//    sprintf(feature_list_12[11], "#scale[0.6]{#bf{Jet y:}}");
//    sprintf(feature_list_12[12], "#scale[0.6]{#bf{Jet #rho:}}");
//
//    double x_min  = th1d_simple_correction->GetXaxis()->GetXmin();
//    double x_max  = th1d_simple_correction->GetXaxis()->GetXmax();
//    double x_bins = th1d_simple_correction->GetNbinsX();
//    double max_arr[4] = {th1d_simple_correction->GetMaximum(), th1d_neural_network->GetMaximum(), th1d_random_forest->GetMaximum(), th1d_linear_regression->GetMaximum()};
//    if (y_max < 0 ) {
//        for ( int i = 0 ; i < 4 ; i++ ) {
//            if ( y_max < max_arr[i] ) y_max = max_arr[i];
//        }
//        if (bool_compareToPaper) y_max = 0.15;
//    }
//
//    TH1D* th1d_canvas_plot = new TH1D("th1d_canvas_plot", plot_title_xlabel_ylabel, x_bins, x_min, x_max);
//    th1d_canvas_plot->SetAxisRange(0, 1.2 * y_max, "Y");
//
//    // Fit to Each Plot
//    TF1* tf1_simple_correction_fit  = new TF1("tf1_simple_correction_fit",  "gaus", x_min, x_max);
//    TF1* tf1_linear_regression_fit  = new TF1("tf1_linear_regression_fit",  "gaus", x_min, x_max);
//    TF1* tf1_random_forest_fit      = new TF1("tf1_random_forest_fit",      "gaus", x_min, x_max);
//    TF1* tf1_neural_network_fit     = new TF1("tf1_neural_network_fit",     "gaus", x_min, x_max);
//
//    th1d_simple_correction  ->Fit(tf1_simple_correction_fit,"RNL");
//    th1d_linear_regression  ->Fit(tf1_linear_regression_fit,"RNL");
//    th1d_random_forest      ->Fit(tf1_random_forest_fit,"RNL");
//    th1d_neural_network     ->Fit(tf1_neural_network_fit,"RNL");
//
//    // TCanvas("name", "title", width (px), height (px))
//    TCanvas* canvas = new TCanvas("canvas", "", 1000, 600);
//
//    // Turns on ticks, turns off stats
//    gPad->SetTicks();
//    gStyle->SetOptStat(0);
//
//    // Turns on Latex formatting relative to canvas
//    TLatex* latex = new TLatex();
//    latex->SetNDC(kTRUE);
//
//    // Build histogram plot
//    th1d_simple_correction->SetLineColor(plot_red);
//    th1d_simple_correction->SetMarkerColor(plot_red);
//    th1d_simple_correction->SetMarkerStyle(mark_circ_open[0]);
//    th1d_simple_correction->SetMarkerSize(mark_circ_open[1]);
//    tf1_simple_correction_fit->SetLineColor(jet_red_line);
//    tf1_simple_correction_fit->SetLineWidth(1);
//
//    th1d_linear_regression->SetLineColor(plot_blue);
//    th1d_linear_regression->SetMarkerColor(plot_blue);
//    th1d_linear_regression->SetMarkerStyle(mark_squa_open[0]);
//    th1d_linear_regression->SetMarkerSize(mark_squa_open[1]);
//    tf1_linear_regression_fit->SetLineColor(jet_blu_line);
//    tf1_linear_regression_fit->SetLineWidth(1);
//
//    th1d_random_forest->SetLineColor(plot_green);
//    th1d_random_forest->SetMarkerColor(plot_green);
//    th1d_random_forest->SetMarkerStyle(mark_diam_open[0]);
//    th1d_random_forest->SetMarkerSize(mark_diam_open[1]);
//    tf1_random_forest_fit->SetLineColor(jet_tea_line);
//    tf1_random_forest_fit->SetLineWidth(1);
//
//    th1d_neural_network->SetLineColor(plot_black);
//    th1d_neural_network->SetMarkerColor(plot_black);
//    th1d_neural_network->SetMarkerStyle(mark_star_open[0]);
//    th1d_neural_network->SetMarkerSize(mark_star_open[1]);
//    tf1_neural_network_fit->SetLineColor(jet_blk_line);
//    tf1_neural_network_fit->SetLineWidth(1);
//
//    // Draws the histograms and lines
//    th1d_canvas_plot->Draw();
//
//    // Draws results from ML paper
//    TF1* tf1_paper_simple_correction_fit;
//    TF1* tf1_paper_neural_network_fit;
//    TF1* tf1_paper_random_forest_fit;
//    TF1* tf1_paper_linear_regression_fit;
//    if (bool_compareToPaper) {
//        Fit_Paper_Data(
//            tf1_paper_simple_correction_fit,
//            tf1_paper_neural_network_fit,
//            tf1_paper_random_forest_fit,
//            tf1_paper_linear_regression_fit
//        );
//
//        tf1_paper_neural_network_fit->SetLineColor(plot_black);
//        tf1_paper_neural_network_fit->SetLineWidth(1);
//        tf1_paper_neural_network_fit->SetLineStyle(3);
//
//        tf1_paper_random_forest_fit->SetLineColor(plot_green);
//        tf1_paper_random_forest_fit->SetLineWidth(1);
//        tf1_paper_random_forest_fit->SetLineStyle(3);
//
//        tf1_paper_linear_regression_fit->SetLineColor(plot_blue);
//        tf1_paper_linear_regression_fit->SetLineWidth(1);
//        tf1_paper_linear_regression_fit->SetLineStyle(3);
//
//        tf1_paper_simple_correction_fit->SetLineColor(plot_red);
//        tf1_paper_simple_correction_fit->SetLineWidth(1);
//        tf1_paper_simple_correction_fit->SetLineStyle(3);
//
//        tf1_paper_simple_correction_fit ->Draw("same");
//        tf1_paper_neural_network_fit    ->Draw("same");
//        tf1_paper_random_forest_fit     ->Draw("same");
//        tf1_paper_linear_regression_fit ->Draw("same");
//    }
//
//    // Draws plots and fit lines
//    tf1_simple_correction_fit   ->Draw("same");
//    tf1_neural_network_fit      ->Draw("same");
//    tf1_random_forest_fit       ->Draw("same");
//    tf1_linear_regression_fit   ->Draw("same");
//    th1d_simple_correction      ->Draw("same");
//    th1d_neural_network         ->Draw("same");
//    th1d_random_forest          ->Draw("same");
//    th1d_linear_regression      ->Draw("same");
//
//    // Draws Legend as TLegend(xmin, ymin, xmax, ymax)
//    if ( bool_showLegend ) {
//        TLegend* legend;
//        if (bool_compareToPaper) legend = new TLegend(0.15, 0.45, 0.40, 0.85);
//        else legend = new TLegend(0.15, 0.65, 0.40, 0.85);
//        legend->AddEntry(th1d_neural_network,       "Neural Network", "lp");
//        legend->AddEntry(th1d_random_forest,        "Random Forest", "lp");
//        legend->AddEntry(th1d_linear_regression,    "Linear Regression", "lp");
//        legend->AddEntry(th1d_simple_correction,    "Area Correction", "lp");
//        if (bool_compareToPaper) {
//            legend->AddEntry(tf1_paper_neural_network_fit,      "NN (Paper)", "l");
//            legend->AddEntry(tf1_paper_random_forest_fit,       "RF (Paper)", "l");
//            legend->AddEntry(tf1_paper_linear_regression_fit,   "LR (Paper)", "l");
//            legend->AddEntry(tf1_paper_simple_correction_fit,   "AC (Paper)", "l");
//        }
//        legend->SetLineWidth(0);
//        legend->SetFillStyle(0);
//        legend->Draw();
//    }
//
//    // Draws Feature List
//    if ( showFeatures == "sigmas" ) {
//        char sigma_vals[5][100];
//        char sigma_list[5][100];
//
//        sprintf(sigma_list[0], "#scale[0.6]{Fit Widths}");
//        sprintf(sigma_list[1], "#scale[0.6]{#bf{#sigma_{Neural Network}}}");
//        sprintf(sigma_list[2], "#scale[0.6]{#bf{#sigma_{Random Forest}}}");
//        sprintf(sigma_list[3], "#scale[0.6]{#bf{#sigma_{Linear Regression}}}");
//        sprintf(sigma_list[4], "#scale[0.6]{#bf{#sigma_{Area Correction}}}");
//
//        sprintf(sigma_vals[0], "");
//        sprintf(sigma_vals[1], "#scale[0.6]{#bf{%1.3f}}", tf1_neural_network_fit    ->GetParameter(2));
//        sprintf(sigma_vals[2], "#scale[0.6]{#bf{%1.3f}}", tf1_random_forest_fit     ->GetParameter(2));
//        sprintf(sigma_vals[3], "#scale[0.6]{#bf{%1.3f}}", tf1_linear_regression_fit ->GetParameter(2));
//        sprintf(sigma_vals[4], "#scale[0.6]{#bf{%1.3f}}", tf1_simple_correction_fit ->GetParameter(2));
//
//        for ( int i = 0 ; i <= 5 ; i++ ) {
//            latex->DrawLatex(0.63, (.813 - 0.0494 * i), sigma_list[i]);
//        }
//        for ( int i = 0 ; i <= 5 ; i++ ) {
//            latex->DrawLatex(0.80, (.813 - 0.0494 * i), sigma_vals[i]);
//        }
//    }
//    if ( showFeatures == "features" ) {
//
//        TTree* coeff_tree = new TTree("LR_Coeffs");
//        vector<string> branch_array = coeff_tree->GetListOfBranches();
//        vector<auto>   coeff_array;
//        for ( int i=0 ; i < branch_array.size() ; i++ ) {
//            coeff_tree->;
//        }
//
//
//        int  feature_count = th1d_feature_list->GetNbinsX();
//        char feature_vals_1[2][100];
//        char feature_vals_3[4][100];
//        char feature_vals_8[9][100];
//        char feature_vals_12[13][100];
//
//        if (feature_count == 1) {
//            sprintf(feature_vals_1[0], "#scale[0.6]{}");
//            sprintf(feature_vals_1[1], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(1));
//
//            for ( int i = 0 ; i <= 1 ; i++ ) {
//                latex->DrawLatex(0.63, (.813 - 0.0494 * i), feature_list_1[i]);
//            }
//            for ( int i = 0 ; i <= 1 ; i++ ) {
//                latex->DrawLatex(0.80, (.813 - 0.0494 * i), feature_vals_3[i]);
//            }
//        }
//
//        if (feature_count == 3) {
//            sprintf(feature_vals_3[0], "#scale[0.6]{}");
//            sprintf(feature_vals_3[1], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(1));
//            sprintf(feature_vals_3[2], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(2));
//            sprintf(feature_vals_3[3], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(3));
//
//            for ( int i = 0 ; i <= 3 ; i++ ) {
//                latex->DrawLatex(0.63, (.813 - 0.0494 * i), feature_list_3[i]);
//            }
//            for ( int i = 0 ; i <= 3 ; i++ ) {
//                latex->DrawLatex(0.80, (.813 - 0.0494 * i), feature_vals_3[i]);
//            }
//        }
//        if (feature_count == 8) {
//            sprintf(feature_vals_8[0], "#scale[0.6]{}");
//            sprintf(feature_vals_8[1], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(1));
//            sprintf(feature_vals_8[2], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(2));
//            sprintf(feature_vals_8[3], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(3));
//            sprintf(feature_vals_8[4], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(4));
//            sprintf(feature_vals_8[5], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(5));
//            sprintf(feature_vals_8[6], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(6));
//            sprintf(feature_vals_8[7], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(7));
//            sprintf(feature_vals_8[8], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(8));
//
//            for ( int i = 0 ; i <= 8 ; i++ ) {
//                latex->DrawLatex(0.63, (.813 - 0.0494 * i), feature_list_8[i]);
//            }
//            for ( int i = 0 ; i <= 8 ; i++ ) {
//                latex->DrawLatex(0.80, (.813 - 0.0494 * i), feature_vals_8[i]);
//            }
//        }
//        if (feature_count == 12) {
//            sprintf(feature_vals_12[0], "#scale[0.6]{}");
//            sprintf(feature_vals_12[1], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(1));
//            sprintf(feature_vals_12[2], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(2));
//            sprintf(feature_vals_12[3], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(3));
//            sprintf(feature_vals_12[4], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(4));
//            sprintf(feature_vals_12[5], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(5));
//            sprintf(feature_vals_12[6], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(6));
//            sprintf(feature_vals_12[7], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(7));
//            sprintf(feature_vals_12[8], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(8));
//            sprintf(feature_vals_12[9], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(7));
//            sprintf(feature_vals_12[10], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(8));
//            sprintf(feature_vals_12[11], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(7));
//            sprintf(feature_vals_12[12], "#scale[0.6]{#bf{%1.3f}}", th1d_feature_list->GetBinContent(8));
//
//            for ( int i = 0 ; i <= 12 ; i++ ) {
//                latex->DrawLatex(0.63, (.813 - 0.0494 * i), feature_list_12[i]);
//            }
//            for ( int i = 0 ; i <= 12 ; i++ ) {
//                latex->DrawLatex(0.80, (.813 - 0.0494 * i), feature_vals_12[i]);
//            }
//        }
//    }
//
//    // Prints out the plot
//    char plot_output[300];
//    sprintf(plot_output, "%s/%s", plot_directory, plot_filename);
//    canvas->Print(plot_output);
//    std::cout << "Plotted " << plot_filename << std::endl;
//    gPad->SetLogy(0);
//
//    return 0;
//}
//
//
//
//
//
//void Jet_ML_Plotter(
//    char* input_file_name,
//    char* th1d_name_arr[5],
//    char* plot_file_name,
//    char* plot_dir_name,
//    float pt_min, float pt_max,
//    bool bool_showLegend,
//    bool bool_showFeatures,
//    bool bool_normalize,
//    bool bool_compareToPaper,
//    char* truth_source
//    ) {
//
//    // Opens and reads the Root output file
//    char input_file_path[200];
//    sprintf(input_file_path, "%s/%s", dir_data, input_file_name);
//    TFile* input_file = new TFile(input_file_path, "READ");
//    std::cout << "Input file read." << std::endl;
//
//    char subdir_plots [200];
//    sprintf(subdir_plots, "%s/MachineLearning/%s", dir_plots, plot_dir_name);
//    std::__fs::filesystem::create_directories(subdir_plots);
//
//    // EXAMINING FEATURES
//    char  plot_labels[400];
//    char  th1d_feature_importance[200];
//    TH1D* th1d_feature_importance_list;
//    char  th1d_simple_correction[200];
//    TH1D* th1d_simple_correction_plot;
//    char  th1d_linear_regression[200];
//    TH1D* th1d_linear_regression_plot;
//    char  th1d_random_forest[200];
//    TH1D* th1d_random_forest_plot;
//    char  th1d_neural_network[200];
//    TH1D* th1d_neural_network_plot;
//
//    // Plot 3 input features, pt_true = PYTHIA, with comparison to paper
//    std::cout << "----- Generating: " << plot_file_name << " -----" << std::endl;
//
//    if ( bool_normalize ) {
//        if ( truth_source == "Pythia" ) {
//            sprintf(plot_labels, "Jet p_{T} Delta for %i-%i GeV #left[p_{T}^{true} = p_{T, ch jet}^{PYTHIA} #right]; p_{T, ch jet}^{reco} - p_{T, ch jet}^{true} [GeV]; Probability Density", int(pt_min), int(pt_max));
//        }
//        else if ( truth_source == "Paper" ) {
//            sprintf(plot_labels, "Jet p_{T} Delta for %i-%i GeV #left[p_{T}^{true} = p_{T}^{raw} #Sigma p_{T, const. i}^{PYTHIA} / #Sigma p_{T, const. i} #right]; p_{T, ch jet}^{reco} - p_{T, ch jet}^{true} [GeV]; Probability Density", int(pt_min), int(pt_max));
//        }
//    }
//    else {
//        if ( truth_source == "Pythia" ) {
//            sprintf(plot_labels, "Jet p_{T} Delta for %i-%i GeV #left[p_{T}^{true} = p_{T, ch jet}^{PYTHIA} #right]; p_{T, ch jet}^{reco} - p_{T, ch jet}^{true} [GeV]; N_{events}", int(pt_min), int(pt_max));
//        }
//        else if ( truth_source == "Paper" ) {
//            sprintf(plot_labels, "Jet p_{T} Delta for %i-%i GeV #left[p_{T}^{true} = p_{T}^{raw} #Sigma p_{T, const. i}^{PYTHIA} / #Sigma p_{T, const. i} #right]; p_{T, ch jet}^{reco} - p_{T, ch jet}^{true} [GeV]; N_{events}", int(pt_min), int(pt_max));
//        }
//    }
//
//    sprintf(th1d_feature_importance, "%s", th1d_name_arr[0]);
//    th1d_feature_importance_list    = (TH1D*) input_file->Get(th1d_feature_importance);
//    std::cout << "Imported " << th1d_name_arr[0] << std::endl;
//    sprintf(th1d_simple_correction, "%s", th1d_name_arr[1]);
//    th1d_simple_correction_plot     = (TH1D*) input_file->Get(th1d_simple_correction);
//
//    TCanvas* canvas_temp = new TCanvas("canvas", "", 1000, 600);
//    th1d_simple_correction_plot->Draw();
//
//    sprintf(th1d_linear_regression, "%s", th1d_name_arr[2]);
//    th1d_linear_regression_plot     = (TH1D*) input_file->Get(th1d_linear_regression);
//
//    sprintf(th1d_random_forest, "%s", th1d_name_arr[3]);
//    th1d_random_forest_plot         = (TH1D*) input_file->Get(th1d_random_forest);
//
//    sprintf(th1d_neural_network, "%s", th1d_name_arr[4]);
//    th1d_neural_network_plot        = (TH1D*) input_file->Get(th1d_neural_network);
//
//    std::cout << "Histograms imported." << std::endl;
//
//    if ( bool_normalize || bool_compareToPaper ) {
//        th1d_simple_correction_plot     ->Scale( 1. / th1d_simple_correction_plot->Integral(),"WIDTH");
//        th1d_linear_regression_plot     ->Scale( 1. / th1d_linear_regression_plot->Integral(),"WIDTH");
//        th1d_random_forest_plot         ->Scale( 1. / th1d_random_forest_plot->Integral(),"WIDTH");
//        th1d_neural_network_plot        ->Scale( 1. / th1d_neural_network_plot->Integral(),"WIDTH");
//    }
//
//    std::cout << "Histograms normalized." << std::endl;
//
//    char* char_showFeatures = "";
//    if ( bool_showFeatures ) char_showFeatures = "features";
//
//    Plot_ML_pT_Comparison(
//        plot_file_name,
//        subdir_plots,
//        plot_labels,
//        th1d_simple_correction_plot,
//        th1d_linear_regression_plot,
//        th1d_random_forest_plot,
//        th1d_neural_network_plot,
//        th1d_feature_importance_list,
//        -1,
//        bool_showLegend, // show legend on plot
//        char_showFeatures, // show feature list
//        bool_compareToPaper); // compare to paper plots
//
//    input_file->Close();
//
//    delete input_file;
//}

//void Jet_Test_Range_Plotter(
//    char  input_file_path[400],
//    char  plot_file_path[400],
//    char  plot_title[100],
//    char  lr_tree_name[100],
//    char  rf_tree_name[100],
//    char  mlp_tree_name[100],
//    char  coeff_tree_name[100],
//    float train_pt_min,
//    float train_pt_max,
//    float test_pt_min,
//    float test_pt_max
//    ) {
//
//    float x_min = train_pt_min - 20.;
//    float x_max = train_pt_max + 20.;
//    int x_bins = int((x_max - x_min) / 1);
//
//    // Opens and reads the Root output file
//    char input_file_path[200];
//    sprintf(input_file_path, "%s/%s", dir_data, input_file_name);
//    TFile* input_file = new TFile(input_file_path, "READ");
//    TTree* input_tree = (TTree*) input_file->Get("Tree_ML");
//    std::cout << "Input file read." << std::endl;
//
//    char subdir_plots [200];
//    sprintf(subdir_plots, "%s/MachineLearning/%s", dir_plots, plot_dir_name);
//    std::__fs::filesystem::create_directories(subdir_plots);
//
//    float jet_area;
//    float jet_pt_raw;
//    float jet_pt_true;
//    float jet_pt_corr;
//    float jet_pt_ml_lr;
//    float jet_pt_ml_rf;
//    float jet_pt_ml_nn;
//
//    input_tree->SetBranchAddress("jet_area", &jet_area);
//    input_tree->SetBranchAddress("jet_pt_raw", &jet_pt_raw);
//    input_tree->SetBranchAddress("jet_pt_true", &jet_pt_true);
//    input_tree->SetBranchAddress("jet_pt_corr", &jet_pt_corr);
//    input_tree->SetBranchAddress("jet_pt_ml_lr", &jet_pt_ml_lr);
//    input_tree->SetBranchAddress("jet_pt_ml_rf", &jet_pt_ml_rf);
//    input_tree->SetBranchAddress("jet_pt_ml_nn", &jet_pt_ml_nn);
//
//    std::cout << "MADE IT HERE!" << std::endl;
//
//    char plot_labels[200];
//    if ( truth_source == "Pythia" ) {
//        sprintf(plot_labels, "Jet p_{T} Delta for %i < p_{T, ch jet}^{true} < %i GeV (Trained on %i < p_{T, ch jet}^{true} < %i GeV); p_{T, ch jet}^{reco} [GeV]; Probability Density", int(test_pt_min), int(test_pt_max), int(train_pt_min), int(train_pt_max));
//    }
//    else if ( truth_source == "Paper" ) {
//        sprintf(plot_labels, "Jet p_{T} Delta for %i < p_{T, ch jet}^{true} < %i GeV (Trained on %i < p_{T, ch jet}^{true} < %i GeV); p_{T, ch jet}^{reco} [GeV]; Probability Density", int(test_pt_min), int(test_pt_max), int(train_pt_min), int(train_pt_max));
//    }
//
//    TH1D* th1d_simple_correction_plot   = new TH1D("th1d_simple_correction_plot", "Machine Learning vs. Simple Correction; p_{T, ch jet}^{reco} [GeV]; Probability Density", x_bins, x_min, x_max);
//    TH1D* th1d_linear_regression_plot   = new TH1D("th1d_linear_regression_plot", "Machine Learning vs. Simple Correction; p_{T, ch jet}^{reco} [GeV]; Probability Density", x_bins, x_min, x_max);
//    TH1D* th1d_random_forest_plot       = new TH1D("th1d_random_forest_plot", "Machine Learning vs. Simple Correction; p_{T, ch jet}^{reco} [GeV]; Probability Density", x_bins, x_min, x_max);
//    TH1D* th1d_neural_network_plot      = new TH1D("th1d_neural_network_plot", "Machine Learning vs. Simple Correction; p_{T, ch jet}^{reco} [GeV]; Probability Density", x_bins, x_min, x_max);
//    TH1D* th1d_truth_plot               = new TH1D("th1d_truth_plot", "Machine Learning vs. Simple Correction; p_{T, ch jet}^{reco} [GeV]; Probability Density", x_bins, x_min, x_max);
//
//    // Empty placeholder
//    TH1D* th1d_feature_importance_list = new TH1D("", "", 1, 0, 1);
//
//    // Fill histograms
//    int event_n = input_tree->GetEntries();
//
//    for ( int e = 0 ; e < event_n ; e++ ) {
//        input_tree->GetEntry(e);
//        if ( test_pt_min < jet_pt_true && test_pt_max > jet_pt_true ){
//            th1d_simple_correction_plot ->Fill(jet_pt_corr);
//            th1d_linear_regression_plot ->Fill(jet_pt_ml_lr);
//            th1d_random_forest_plot     ->Fill(jet_pt_ml_rf);
//            th1d_neural_network_plot    ->Fill(jet_pt_ml_nn);
//            th1d_truth_plot             ->Fill(jet_pt_true);
//        }
//    }
//
//    th1d_simple_correction_plot     ->Scale( 1. / th1d_simple_correction_plot->Integral(),"WIDTH");
//    th1d_linear_regression_plot     ->Scale( 1. / th1d_linear_regression_plot->Integral(),"WIDTH");
//    th1d_random_forest_plot         ->Scale( 1. / th1d_random_forest_plot->Integral(),"WIDTH");
//    th1d_neural_network_plot        ->Scale( 1. / th1d_neural_network_plot->Integral(),"WIDTH");
//    th1d_truth_plot                 ->Scale( 1. / th1d_truth_plot->Integral(),"WIDTH");
//
//    Plot_ML_pT_Comparison(
//        plot_file_name,
//        subdir_plots,
//        plot_labels,
//        th1d_simple_correction_plot,
//        th1d_linear_regression_plot,
//        th1d_random_forest_plot,
//        th1d_neural_network_plot,
//        th1d_feature_importance_list,
//        0.22, // ymax
//        true, // show legend on plot
//        "sigmas", // show feature list
//        false); // compare to paper plots
//
//}
//
//
//
//void Jet_Truth_Plotter(
//    char* input_file_name,
//    char* plot_file_name,
//    char* plot_dir_name,
//    float train_pt_min, float train_pt_max,
//    char* truth_source
//    ) {
//
//    float x_min = train_pt_min - 10.;
//    float x_max = train_pt_max + 10.;
//    int x_bins = int((x_max - x_min) / 1);
//
//    // Opens and reads the Root output file
//    char input_file_path[200];
//    sprintf(input_file_path, "%s/%s", dir_data, input_file_name);
//    TFile* input_file = new TFile(input_file_path, "READ");
//    TTree* input_tree = (TTree*) input_file->Get("Tree_ML");
//    std::cout << "Input file read." << std::endl;
//
//    char subdir_plots [200];
//    sprintf(subdir_plots, "%s/MachineLearning/%s", dir_plots, plot_dir_name);
//    std::__fs::filesystem::create_directories(subdir_plots);
//
//    float jet_pt_true;
//    input_tree->SetBranchAddress("jet_pt_true", &jet_pt_true);
//
//    std::cout << "MADE IT HERE!" << std::endl;
//
//    char plot_labels[200];
//    if ( truth_source == "Pythia" ) {
//        sprintf(plot_labels, "Distribution of Jet p_{T}^{True} for %i < p_{T, ch jet}^{true} < %i GeV; p_{T, ch jet}^{reco} [GeV]; Probability Density", int(train_pt_min), int(train_pt_max));
//    }
//    else if ( truth_source == "Paper" ) {
//        sprintf(plot_labels, "Distribution of Jet p_{T}^{True} for %i < p_{T, ch jet}^{true} < %i GeV; p_{T, ch jet}^{reco} [GeV]; Probability Density", int(train_pt_min), int(train_pt_max));
//    }
//
//    TH1D* th1d_truth_plot   = new TH1D("th1d_pt_true_plot", "Distribution of p_{T}^{True}; p_{T, ch jet}^{True} [GeV]; Probability Density", x_bins, x_min, x_max);
//
//    // Fill histograms
//    int event_n = input_tree->GetEntries();
//
//    for ( int e = 0 ; e < event_n ; e++ ) {
//        input_tree->GetEntry(e);
//        th1d_truth_plot ->Fill(jet_pt_true);
//    }
//
//    th1d_truth_plot ->Scale( 1. / th1d_truth_plot->Integral(),"WIDTH");
//
//    TH1D* th1d_canvas_plot = new TH1D("th1d_canvas_plot", plot_labels, x_bins, x_min, x_max);
//
//    // TCanvas("name", "title", width (px), height (px))
//    TCanvas* canvas = new TCanvas("canvas", "", 1000, 600);
//    th1d_canvas_plot->SetAxisRange(0, .15, "Y");
//
//    // Turns on ticks, turns off stats
//    gPad->SetTicks();
//    gStyle->SetOptStat(0);
//
//    // Turns on Latex formatting relative to canvas
//    TLatex* latex = new TLatex();
//    latex->SetNDC(kTRUE);
//
//    // Build histogram plot
//    th1d_truth_plot->SetLineColor(jet_blk_line);
//    th1d_truth_plot->SetMarkerColor(jet_blk_mark);
//    th1d_truth_plot->SetMarkerStyle(mark_circ_open[0]);
//    th1d_truth_plot->SetMarkerSize(mark_circ_open[1]);
//
//    // Draws the histograms and lines
//    th1d_canvas_plot->Draw();
//
//    th1d_truth_plot ->Draw("same");
//
//    // Draws Legend as TLegend(xmin, ymin, xmax, ymax)
//    TLegend* legend;
//    legend = new TLegend(0.15, 0.80, 0.3, 0.85);
//    legend->AddEntry(th1d_truth_plot, "p_{T}^{True}", "p");
//    legend->SetLineWidth(0);
//    legend->SetFillStyle(0);
//    legend->Draw();
//
//    // Prints out the plot
//    char plot_output[300];
//    sprintf(plot_output, "%s/MachineLearning/%s/%s", dir_plots, plot_dir_name, plot_file_name);
//    canvas->Print(plot_output);
//    std::cout << "Plotted " << plot_file_name << std::endl;
//    gPad->SetLogy(0);
//
//    return 0;
//}
//
//
//
//void Jet_ptTrue_Distribution_Plotter (
//    char* input_file_name,
//    char* plot_tree_name,
//    char* plot_dir_name,
//    float train_pt_min, float train_pt_max
//    ) {
//
//    float x_min = train_pt_min;
//    float x_max = train_pt_max;
//    int x_bins = int((x_max - x_min) / 1);
//    char plot_file_name[100];
//    sprintf(plot_file_name, "%s.pdf", plot_tree_name);
//
//    // Opens and reads the Root output file
//    char input_file_path[200];
//    sprintf(input_file_path, "%s/%s", dir_data, input_file_name);
//    TFile* input_file = new TFile(input_file_path, "READ");
//    TTree* input_tree = (TTree*) input_file->Get("Tree_ML");
//    std::cout << "Input file read." << std::endl;
//
//    char subdir_plots [200];
//    sprintf(subdir_plots, "%s/MachineLearning/%s", dir_plots, plot_dir_name);
//    std::__fs::filesystem::create_directories(subdir_plots);
//
//    TH1D* th1d_truth_distribution_plot     = (TH1D*) input_file->Get(plot_tree_name);
//
//    std::cout << "MADE IT HERE!" << std::endl;
//
//    // TCanvas("name", "title", width (px), height (px))
//    TCanvas* canvas = new TCanvas("canvas", "", 1000, 600);
//    th1d_truth_distribution_plot->SetAxisRange(100, 1000000, "Y");
//
//    // Turns on ticks, turns off stats
//    gPad->SetTicks();
//    gStyle->SetOptStat(0);
//    gPad->SetLogy(1);
//
//    // Turns on Latex formatting relative to canvas
//    TLatex* latex = new TLatex();
//    latex->SetNDC(kTRUE);
//
//    // Build histogram plot
//    th1d_truth_distribution_plot->SetLineColor(jet_blk_line);
//    th1d_truth_distribution_plot->SetMarkerColor(jet_blk_mark);
//    th1d_truth_distribution_plot->SetMarkerStyle(mark_circ_open[0]);
//    th1d_truth_distribution_plot->SetMarkerSize(mark_circ_open[1]);
//
//    // Draws the histograms and lines
//    th1d_truth_distribution_plot->Draw();
//
//    // Draws Legend as TLegend(xmin, ymin, xmax, ymax)
//    TLegend* legend;
//    legend = new TLegend(0.15, 0.80, 0.3, 0.85);
//    legend->AddEntry(th1d_truth_distribution_plot, "p_{T}^{True}", "lp");
//    legend->SetLineWidth(0);
//    legend->SetFillStyle(0);
//    legend->Draw();
//
//    // Prints out the plot
//    char plot_output[300];
//    sprintf(plot_output, "%s/MachineLearning/%s/%s", dir_plots, plot_dir_name, plot_file_name);
//    canvas->Print(plot_output);
//    std::cout << "Plotted " << plot_file_name << std::endl;
//    gPad->SetLogy(0);
//
//    return 0;
//}

} // End of namespace for Jet_Plotter
