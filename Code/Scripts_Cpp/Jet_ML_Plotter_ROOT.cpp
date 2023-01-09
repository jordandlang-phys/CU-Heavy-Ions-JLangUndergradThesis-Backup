#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <vector>
#include <sstream>

#include "jet_ml_constants.h"
using namespace Jet_ML_Constants;

namespace Jet_Plotter {

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
    
    if (!use_delta) legend->AddEntry(th1f_jet_pt_true,  "p_{T}^{True} (0.5 scale)", "lp");
    legend->AddEntry(th1f_jet_pt_reco,  "Linear Regression", "lp");
    legend->AddEntry((TObject*)0, reco_stats, "");
    legend->AddEntry(tg_jet_pt_corr,  "Area Correction", "lp");
    legend->AddEntry((TObject*)0, corr_stats, "");
    legend->SetLineWidth(0);
    legend->SetFillStyle(0);
    
    th1f_canvas_plot->Draw();
    
    if (!use_delta) th1f_jet_pt_true->Draw("same");
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



void Plot_ML_LR_Coefficients_F12(
    char  input_file_path[500],
    char  input_tree_name[10][100],
    char  plot_file_name[100],
    float test_min_max_array[20][2],
    int   bin_count
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
    float arr_jet_pt_corr[20];
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
        
        float jet_pt_raw;
        float jet_pt_corr;
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
        
        coeff_tree->SetBranchAddress("jet_pt_raw",    &jet_pt_raw);
        coeff_tree->SetBranchAddress("jet_pt_corr",   &jet_pt_corr);
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
        coeff_tree->SetBranchAddress("lr_intercept",       &lr_intercept);
        
        coeff_tree->GetEntry(0);
        
        arr_test_pt_avg[i]    = test_pt_avg;
        arr_jet_pt_raw[i]     = jet_pt_raw;
        arr_jet_pt_corr[i]    = jet_pt_corr;
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
        
        delete coeff_tree;
    }
    
    TGraph* tg_jet_pt_raw       = new TGraph(bin_count, arr_test_pt_avg, arr_jet_pt_raw);
    tg_jet_pt_raw->SetTitle("Jet p_{T}^{Raw}; Average Train/Test p_{T} [GeV]; Coefficient Value");
    TGraph* tg_jet_pt_corr      = new TGraph(bin_count, arr_test_pt_avg, arr_jet_pt_corr);
    tg_jet_pt_corr->SetTitle("Jet p_{T}^{Corrected}; Average Train/Test p_{T} [GeV]; Coefficient Value");
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
    
    tg_jet_pt_corr->SetLineColor(kViolet-5);
    tg_jet_pt_corr->SetMarkerColor(kViolet-5);
    tg_jet_pt_corr->SetMarkerStyle(mark_squa_open[0]);
    tg_jet_pt_corr->SetMarkerSize(mark_squa_open[1]);
    
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
    tmg_canvas_plot->Add(tg_lr_intercept, "PL"); // Draws intercept values first to be in background
    tmg_canvas_plot->Add(tg_jet_pt_raw, "PL");
    tmg_canvas_plot->Add(tg_jet_pt_corr, "PL");
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
    TLegend* legend = new TLegend(0.62, 0.12, 0.77, 0.50);
    
    legend->AddEntry(tg_jet_pt_raw,       "Jet p_{T}^{Raw}", "lp");
    legend->AddEntry(tg_jet_pt_corr,      "Jet p_{T}^{Area Correction}", "lp");
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
    legend->AddEntry(tg_lr_intercept,     "Lin. Reg. Intercept", "lp");
    
    legend->SetLineWidth(0);
    legend->SetFillStyle(0);
    
    legend->Draw("same");
    
    canvas->Print(plot_file_name);
    std::cout << "Plotted " << plot_file_name << std::endl;
}

} // End of namespace for Jet_Plotter
