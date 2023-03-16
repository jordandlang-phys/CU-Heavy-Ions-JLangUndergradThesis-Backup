#include "TFile.h"
#include "TTree.h"
#include <string>

#include "../Scripts_Cpp/jet_ml_constants.h"
using namespace Jet_ML_Constants;

TF1* tf1_paper_simple_correction_fit;

double Gaussian_Func(double *x,double *par) {
    double arg = 0;
    if (par[2]!=0) arg = (x[0] - par[1])/par[2];
    double fitval = par[0]*TMath::Exp(-0.5*arg*arg);
    return fitval;
}

void Fit_Paper_Data(
    float x_min,
    float x_max
    ) {
    // Create Graphs to Fit Paper Data
    double tgraph_paper_sc_x[25] = {-21., -19., -17., -15., -13., -11., -9., -7., -5., -3., -1., 1., 3., 5., 7., 9., 11., 12., 15., 17., 19., 21., 23., 25., 27.};
    double tgraph_paper_sc_y[25] = {-0.0003636363636363440, 0.0007272727272727430, 0.0018181818181818300, 0.0040000000000000000, 0.009090909090909090, 0.015272727272727300, 0.024363636363636400, 0.034181818181818200, 0.04000000000000000, 0.048727272727272700, 0.049090909090909100, 0.05127272727272730, 0.048000000000000000, 0.04436363636363640, 0.03672727272727280, 0.029454545454545500, 0.02145454545454550, 0.015636363636363600, 0.010545454545454600, 0.00545454545454549, 0.0032727272727272900, 0.0029090909090909200, 0.00036363636363637200, 0.00036363636363637200, 0.0};
    
    TGraph* tgraph_paper_simple_correction = new TGraph(25, tgraph_paper_sc_x, tgraph_paper_sc_y);
    tf1_paper_simple_correction_fit = new TF1("tf1_paper_simple_correction_fit", "gaus", -40., 40.);
    tf1_paper_simple_correction_fit->SetParameters(.1, 0., 5.);
    tgraph_paper_simple_correction->Fit(tf1_paper_simple_correction_fit,"RNL");
}

void Plot_AreaCorr_Compare(
    vector<vector<string>> input_list, // file_name, tree_name, variable, plot_legend_label, plot_color, plot_symbol, use_integer_variable
    vector<string> plot_title,
    string plot_axes,
    string output_path,
    int    x_bins,
    float  x_min,
    float  x_max,
    bool   y_log,
    float  y_min,
    float  y_max,
    bool   use_normalized = false,
    bool   show_widths = false
    ) {
    
    TH1D* th1d_canvas_plot = new TH1D("canvas_plot", (plot_title[0] + ";" + plot_axes).c_str(), x_bins, x_min, x_max);
        th1d_canvas_plot->SetAxisRange(y_min, y_max, "Y");
    
    // TCanvas("name", "title", width (px), height (px))
    TCanvas* canvas = new TCanvas("canvas", "", 1000, 600);
    gPad->SetTicks();
    gStyle->SetOptStat(0);
    if ( y_log ) gPad->SetLogy(1);

    // Turns on Latex formatting relative to canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(kTRUE);
    TLegend* legend;
    if ( show_widths ) legend = new TLegend(0.60, 0.85 - (0.10 + 0.10 * input_list.size()), 0.85, 0.85);
    else legend = new TLegend(0.60, 0.85 - (0.05 + 0.05 * input_list.size()), 0.85, 0.85);
    th1d_canvas_plot->Draw();
    
    Fit_Paper_Data(x_min, x_max);
    std::cout << "Paper area-correction is fit!" << std::endl;
    tf1_paper_simple_correction_fit->SetLineColor(plot_black);
    tf1_paper_simple_correction_fit->SetLineWidth(1);
    tf1_paper_simple_correction_fit->SetLineStyle(3);
    tf1_paper_simple_correction_fit->Draw("same");
    std::cout << "Paper area-correction is styled!" << std::endl;
    
    // Add plots
    for ( int i=0 ; i < input_list.size() ; i++ ) {
        string input_file_path      = input_list[i][0];
        string input_tree_name      = input_list[i][1];
        string input_variable       = input_list[i][2];
        string plot_legend_label    = input_list[i][3];
        string plot_color_str       = input_list[i][4];
        string plot_symbol_str      = input_list[i][5];
        
        TFile* input_file = new TFile(input_file_path.c_str(), "READ");
        std::cout << "Input file accessed: " << input_file_path << std::endl;
        TTree* input_tree = (TTree*) input_file->Get(input_tree_name.c_str());
        std::cout << "Input tree accessed: " << input_tree_name << std::endl;
        
        float plot_var;
        float jet_pt_true;
        input_tree->SetBranchAddress(input_variable.c_str(), &plot_var);
        input_tree->SetBranchAddress("jet_pt_true", &jet_pt_true);
        std::cout << "Input variable accessed: " << input_variable << std::endl;
        
        TH1D* th1d_plot_var = new TH1D("plot_var", "", x_bins, x_min, x_max);
        th1d_plot_var->Sumw2();
        
        std::cout << "Filling histogram..." << std::endl;
        for ( int i=0 ; i < input_tree->GetEntries() ; i++ ) {
            input_tree->GetEntry(i);
            th1d_plot_var->Fill(plot_var-jet_pt_true);
        }
        if ( use_normalized ) th1d_plot_var->Scale( 1. / th1d_plot_var->Integral(),"WIDTH");
        std::cout << "Histogram filled!" << std::endl;
        
        int plot_color;
        if ( plot_color_str.compare("plot_red") == 0 ) plot_color = plot_red;
        else if ( plot_color_str.compare("plot_green") == 0 ) plot_color = plot_green;
        else if ( plot_color_str.compare("plot_blue") == 0 ) plot_color = plot_blue;
        else if ( plot_color_str.compare("plot_violet") == 0 ) plot_color = plot_violet;
        else plot_color = plot_black;
        
        float plot_symbol[2];
        if ( plot_symbol_str.compare("circ_fill") == 0 ) {
            plot_symbol[0] = mark_circ_fill[0];
            plot_symbol[1] = mark_circ_fill[1];
        }
        else {
            plot_symbol[0] = mark_circ_open[0];
            plot_symbol[1] = mark_circ_open[1];
        }
        
        th1d_plot_var->SetLineColor(plot_color);
        th1d_plot_var->SetMarkerColor(plot_color);
        th1d_plot_var->SetMarkerStyle(plot_symbol[0]);
        th1d_plot_var->SetMarkerSize(plot_symbol[1]);
        
        legend->AddEntry(th1d_plot_var, plot_legend_label.c_str(), "lp");
        if ( show_widths ) {
            char var_stats[100];
            float var_mean   = th1d_plot_var->GetMean();
            float var_stddev = th1d_plot_var->GetStdDev();
            snprintf(var_stats, 100, "#Delta: %.2f, #sigma: %.2f", var_mean, var_stddev);
            legend->AddEntry((TObject*)0, var_stats, "");
        }
        th1d_plot_var->Draw("same");
        std::cout << "Histogram added to plot.\n" << std::endl;
    }
    std::cout << "Adding Legend" << std::endl;
    
    legend->AddEntry(tf1_paper_simple_correction_fit, "p_{T, Jet}^{Area Corr.} (Haake, 2019)", "l");
    if ( show_widths ) {
        char var_stats[100];
        float var_mean   = tf1_paper_simple_correction_fit->GetParameter(1);
        float var_stddev = tf1_paper_simple_correction_fit->GetParameter(2);
        snprintf(var_stats, 100, "#Delta: %.2f, #sigma: %.2f", var_mean, var_stddev);
        legend->AddEntry((TObject*)0, var_stats, "");
    }
    
    legend->SetLineWidth(0);
    legend->SetFillStyle(0);
    legend->Draw("same");
    
    for ( int i=0 ; i < plot_title.size() ; i++ ) {
        latex->DrawLatex(0.15, 0.81 - (0.06*i), ("#scale[0.75]{" + plot_title[i] + "}").c_str());
    }
    
    canvas->Print(output_path.c_str());
    std::cout << "Plotted: " << output_path << std::endl;
    
    gPad->SetLogy(0);
    
    delete legend;
    delete latex;
    delete canvas;
    delete th1d_canvas_plot;
    std::cout << "Plotting complete!" << std::endl;
}

void Plot_Simple(
    vector<vector<string>> input_list, // file_name, tree_name, variable, plot_legend_label, plot_color, plot_symbol
    vector<string> plot_title,
    string plot_axes,
    string output_path,
    int    x_bins,
    float  x_min,
    float  x_max,
    bool   y_log,
    float  y_min,
    float  y_max,
    bool   use_normalized = false,
    bool   show_widths = false,
    int    plot_width = 1000,
    int    plot_height = 600
    ) {
    
    TH1D* th1d_canvas_plot = new TH1D("canvas_plot", (plot_title[0] + ";" + plot_axes).c_str(), x_bins, x_min, x_max);
        th1d_canvas_plot->SetAxisRange(y_min, y_max, "Y");
    
    // TCanvas("name", "title", width (px), height (px))
    TCanvas* canvas = new TCanvas("canvas", "", plot_width, plot_height);
    gPad->SetTicks();
    gStyle->SetOptStat(0);
    if ( y_log ) gPad->SetLogy(1);

    // Turns on Latex formatting relative to canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(kTRUE);
    TLegend* legend = new TLegend(0.60, 0.85 - (0.05 * input_list.size()), 0.85, 0.85);
    th1d_canvas_plot->Draw();
    
    // Add plots
    for ( int i=0 ; i < input_list.size() ; i++ ) {
        string input_file_path      = input_list[i][0];
        string input_tree_name      = input_list[i][1];
        string input_variable       = input_list[i][2];
        string plot_legend_label    = input_list[i][3];
        string plot_color_str       = input_list[i][4];
        string plot_symbol_str      = input_list[i][5];
        bool   use_var_int = false;
        if ( input_list[i][6].compare("true") == 0 ) use_var_int = true;
        
        TFile* input_file = new TFile(input_file_path.c_str(), "READ");
        std::cout << "Input file accessed: " << input_file_path << std::endl;
        TTree* input_tree = (TTree*) input_file->Get(input_tree_name.c_str());
        std::cout << "Input tree accessed: " << input_tree_name << std::endl;
        
        float plot_var;
        int   plot_var_int;
        if ( use_var_int ) input_tree->SetBranchAddress(input_variable.c_str(), &plot_var_int);
        else input_tree->SetBranchAddress(input_variable.c_str(), &plot_var);
        std::cout << "Input variable accessed: " << input_variable << std::endl;
        
        TH1D* th1d_plot_var = new TH1D("plot_var", "", x_bins, x_min, x_max);
        th1d_plot_var->Sumw2();
        
        std::cout << "Filling histogram..." << std::endl;
        for ( int i=0 ; i < input_tree->GetEntries() ; i++ ) {
            input_tree->GetEntry(i);
            if ( use_var_int ) th1d_plot_var->Fill(plot_var_int);
            else th1d_plot_var->Fill(plot_var);
        }
        if ( use_normalized ) th1d_plot_var->Scale( 1. / th1d_plot_var->Integral(),"WIDTH");
        std::cout << "Histogram filled!" << std::endl;
        
        int plot_color;
        if ( plot_color_str.compare("plot_red") == 0 ) plot_color = plot_red;
        else if ( plot_color_str.compare("plot_green") == 0 ) plot_color = plot_green;
        else if ( plot_color_str.compare("plot_blue") == 0 ) plot_color = plot_blue;
        else if ( plot_color_str.compare("plot_violet") == 0 ) plot_color = plot_violet;
        else plot_color = plot_black;
        
        float plot_symbol[2];
        if ( plot_symbol_str.compare("circ_fill") == 0 ) {
            plot_symbol[0] = mark_circ_fill[0];
            plot_symbol[1] = mark_circ_fill[1];
        }
        else {
            plot_symbol[0] = mark_circ_open[0];
            plot_symbol[1] = mark_circ_open[1];
        }
        
        th1d_plot_var->SetLineColor(plot_color);
        th1d_plot_var->SetMarkerColor(plot_color);
        th1d_plot_var->SetMarkerStyle(plot_symbol[0]);
        th1d_plot_var->SetMarkerSize(plot_symbol[1]);
        
        legend->AddEntry(th1d_plot_var, plot_legend_label.c_str(), "lp");
        if ( show_widths ) {
            char var_stats[100];
            float var_mean   = th1d_plot_var->GetMean();
            float var_stddev = th1d_plot_var->GetStdDev();
            snprintf(var_stats, 100, "#Delta: %.2f, #sigma: %.2f", var_mean, var_stddev);
            legend->AddEntry((TObject*)0, var_stats, "");
        }
        
        th1d_plot_var->Draw("same");
        std::cout << "Histogram added to plot.\n" << std::endl;
    }
    
    legend->SetLineWidth(0);
    legend->SetFillStyle(0);
    
    legend->Draw("same");
    
    for ( int i=0 ; i < plot_title.size() ; i++ ) {
        latex->DrawLatex(0.15, 0.81 - (0.06*i), ("#scale[0.75]{" + plot_title[i] + "}").c_str());
    }
    
    canvas->Print(output_path.c_str());
    std::cout << "Plotted: " << output_path << std::endl;
    
    gPad->SetLogy(0);
    
    delete legend;
    delete latex;
    delete canvas;
    delete th1d_canvas_plot;
    std::cout << "Plotting complete!" << std::endl;
}

const int plot_colors[35] = {kOrange+10, kPink+10, kViolet,
    kAzure+1, kAzure+2, kAzure-5, kAzure-4, kAzure+6, kAzure+5, kAzure-8, kAzure-9,
    kAzure+1, kAzure+2, kAzure-5, kAzure-4, kAzure+6, kAzure+5, kAzure-8, kAzure-9,
    kAzure+1, kAzure+2, kAzure-5, kAzure-4, kAzure+6, kAzure+5, kAzure-8, kAzure-9,
    kAzure+1, kAzure+2, kAzure-5, kAzure-4, kAzure+6, kAzure+5, kAzure-8, kAzure-9};


void Event_2D_Plotter(
    string comb_file_path,
    string comb_tree_name,
    string pyth_file_path,
    string pyth_tree_name,
    string output_dir,
    int first_event = 0,
    int last_event = 10
    ) {
    
    // Generates a directory for output files
    std::__fs::filesystem::create_directories(output_dir.c_str());
    char plot_dir[200];
    sprintf(plot_dir, "%s/Jet_2D_Plots/", output_dir.c_str());
    std::__fs::filesystem::create_directories(plot_dir);
    
    // Variables
    
    const int max_jets = 100;
    const int max_parts = 400;
    
    // Combined file of jets from thermal and PYTHIA data
    TFile* combJet_file = new TFile(comb_file_path.c_str(), "READ");
    TTree* combJet_tree = (TTree*) combJet_file->Get(comb_tree_name.c_str());
    std::cout << "Reading combined data file." << std::endl;
    int     c_jet_n;
    float   c_jet_pt[max_jets];
    float   c_jet_y[max_jets];
    float   c_jet_phi[max_jets];
    int     c_jet_part_n[max_jets];
    float   c_jet_part_pt[max_jets][max_parts];
    float   c_jet_part_eta[max_jets][max_parts];
    float   c_jet_part_phi[max_jets][max_parts];
    combJet_tree->SetBranchAddress("jet_n",         &c_jet_n);
    combJet_tree->SetBranchAddress("jet_pt",        c_jet_pt);
    combJet_tree->SetBranchAddress("jet_y",         c_jet_y);
    combJet_tree->SetBranchAddress("jet_phi",       c_jet_phi);
    combJet_tree->SetBranchAddress("jet_const_n",   c_jet_part_n);
    combJet_tree->SetBranchAddress("jet_const_pt",  c_jet_part_pt);
    combJet_tree->SetBranchAddress("jet_const_eta", c_jet_part_eta);
    combJet_tree->SetBranchAddress("jet_const_phi", c_jet_part_phi);
    
    // File of jets from PYTHIA
    TFile* pythJet_file = new TFile(pyth_file_path.c_str(), "READ");
    TTree* pythJet_tree = (TTree*) pythJet_file->Get(pyth_tree_name.c_str());
    std::cout << "Reading pythia data file." << std::endl;
    int     p_jet_n;
    float   p_jet_pt[max_jets];
    float   p_jet_y[max_jets];
    float   p_jet_phi[max_jets];
    int     p_jet_part_n[max_jets];
    float   p_jet_part_pt[max_jets][max_parts];
    float   p_jet_part_eta[max_jets][max_parts];
    float   p_jet_part_phi[max_jets][max_parts];
    pythJet_tree->SetBranchAddress("jet_n",         &p_jet_n);
    pythJet_tree->SetBranchAddress("jet_pt",        p_jet_pt);
    pythJet_tree->SetBranchAddress("jet_y",         p_jet_y);
    pythJet_tree->SetBranchAddress("jet_phi",       p_jet_phi);
    pythJet_tree->SetBranchAddress("jet_const_n",   p_jet_part_n);
    pythJet_tree->SetBranchAddress("jet_const_pt",  p_jet_part_pt);
    pythJet_tree->SetBranchAddress("jet_const_eta", p_jet_part_eta);
    pythJet_tree->SetBranchAddress("jet_const_phi", p_jet_part_phi);
    
    std::cout << "Input variables initialized." << std::endl;
    
    // TGraph of Jets
    for ( int e = first_event ; e < last_event ; e++ ) {
        combJet_tree->GetEntry(e);
        pythJet_tree->GetEntry(e);

        TGraph* tgraph_embed_jets[max_jets];
        TGraph* tgraph_pythia_jets[max_jets];
        TText*  jet_pt_text[max_jets];
        TText*  jet_pt_outline[max_jets];

        TCanvas* canvas = new TCanvas("canvas", "", 2000, 600);
        gPad->SetTicks();
        gStyle->SetOptStat(0);
        
        TH2D* th2d_jet_pt_distribution = new TH2D(Form("th2d_event%i_jet_pt_distribution", e), "Jet p_{T} and Constituent Distribution; #phi; #eta", 20, 0, 2*3.141593, 18, -0.9, 0.9);
        th2d_jet_pt_distribution->Draw("p, same");
        
        for ( int j = (c_jet_n - 1) ; j >= 0 ; j-- ) {
            float  particle_eta[max_parts], particle_phi[max_parts];

            for ( int p = 0 ; p < c_jet_part_n[j] ; p++ ) {
                particle_eta[p] = c_jet_part_eta[j][p];
                particle_phi[p] = c_jet_part_phi[j][p];
            }
            
            tgraph_embed_jets[j] = new TGraph(c_jet_part_n[j], particle_phi, particle_eta);
            tgraph_embed_jets[j]->SetMarkerStyle(24);
            if ( j == 0 ) tgraph_embed_jets[j]->SetMarkerColorAlpha(plot_colors[j], 0.90);
            else if ( j == 1 ) tgraph_embed_jets[j]->SetMarkerColorAlpha(plot_colors[j], 0.80);
            else if ( j == 2 ) tgraph_embed_jets[j]->SetMarkerColorAlpha(plot_colors[j], 0.65);
            else tgraph_embed_jets[j]->SetMarkerColorAlpha(plot_colors[j], 0.60);
            tgraph_embed_jets[j]->Draw("p, same");
        }
        
        for ( int j = (p_jet_n - 1) ; j >= 0 ; j-- ) {
            float  particle_eta[max_parts], particle_phi[max_parts];

            for ( int p = 0 ; p < p_jet_part_n[j] ; p++ ) {
                particle_eta[p] = p_jet_part_eta[j][p];
                particle_phi[p] = p_jet_part_phi[j][p];
            }
            
            tgraph_pythia_jets[j] = new TGraph(p_jet_part_n[j], particle_phi, particle_eta);
            tgraph_pythia_jets[j]->SetMarkerStyle(20);
            tgraph_pythia_jets[j]->SetMarkerSize(0.7);
            tgraph_pythia_jets[j]->SetMarkerColorAlpha(kBlack, 1);
//            if ( j == 0 ) tgraph_pythia_jets[j]->SetMarkerColorAlpha(plot_colors[j], 0.90);
//            else if ( j == 1 ) tgraph_pythia_jets[j]->SetMarkerColorAlpha(plot_colors[j], 0.80);
//            else if ( j == 2 ) tgraph_pythia_jets[j]->SetMarkerColorAlpha(plot_colors[j], 0.65);
//            else tgraph_pythia_jets[j]->SetMarkerColorAlpha(plot_colors[j], 0.50);
            tgraph_pythia_jets[j]->Draw("p, same");
        }
        
        for ( int j = (c_jet_n - 1) ; j >= 0 ; j-- ) {
            char jet_pt[20];
            float  x_offset = 0.005; // Outline offset in plot axes units
            float  y_offset = 0.005; // Outline offset in plot axes units
            snprintf(jet_pt, 20, "%3.2f GeV", c_jet_pt[j]);
            jet_pt_outline[j] = new TText(c_jet_phi[j] + x_offset, c_jet_y[j] + y_offset, jet_pt);
            jet_pt_outline[j]->SetTextColor(kWhite);
            jet_pt_outline[j]->SetTextSize(0.04);
            jet_pt_outline[j]->SetTextAlign(22);
            jet_pt_text[j] = new TText(c_jet_phi[j], c_jet_y[j], jet_pt);
            jet_pt_text[j]->SetTextColor(kBlack);
            jet_pt_text[j]->SetTextSize(0.04);
            jet_pt_text[j]->SetTextAlign(22);
            jet_pt_outline[j]->Draw("same");
            jet_pt_outline[j]->DrawText(c_jet_phi[j] + x_offset, c_jet_y[j] + y_offset, "same");
            jet_pt_outline[j]->DrawText(c_jet_phi[j] + x_offset, c_jet_y[j] - y_offset, "same");
            jet_pt_outline[j]->DrawText(c_jet_phi[j] - x_offset, c_jet_y[j] + y_offset, "same");
            jet_pt_outline[j]->DrawText(c_jet_phi[j] - x_offset, c_jet_y[j] - y_offset, "same");
            jet_pt_text[j]->Draw("same");
        }
        
        char dir_tgraph[200];
        sprintf(dir_tgraph, "%stgraph_event%i_jets_all.pdf", plot_dir, e);
        canvas->Print(dir_tgraph);
    }
    
    // TGraph of Jets
    for ( int e = first_event ; e < last_event ; e++ ) {
        combJet_tree->GetEntry(e);

        TGraph* tgraph_embed_jets[max_jets];
        TText*  jet_pt_text[max_jets];
        TText*  jet_pt_outline[max_jets];

        TCanvas* canvas = new TCanvas("canvas", "", 2000, 600);
        gPad->SetTicks();
        gStyle->SetOptStat(0);
        
        TH2D* th2d_jet_pt_distribution = new TH2D(Form("th2d_event%i_jet_pt_distribution", e), "Jet p_{T} and Constituent Distribution; #phi; #eta", 20, 0, 2*3.141593, 18, -0.9, 0.9);
        th2d_jet_pt_distribution->Draw("p, same");
        
        for ( int j = (c_jet_n - 1) ; j >= 0 ; j-- ) {
            float  particle_eta[max_parts], particle_phi[max_parts];

            for ( int p = 0 ; p < c_jet_part_n[j] ; p++ ) {
                particle_eta[p] = c_jet_part_eta[j][p];
                particle_phi[p] = c_jet_part_phi[j][p];
            }
            
            tgraph_embed_jets[j] = new TGraph(c_jet_part_n[j], particle_phi, particle_eta);
            tgraph_embed_jets[j]->SetMarkerStyle(24);
            if ( j == 0 ) tgraph_embed_jets[j]->SetMarkerColorAlpha(plot_colors[j], 0.90);
            else if ( j == 1 ) tgraph_embed_jets[j]->SetMarkerColorAlpha(plot_colors[j], 0.80);
            else if ( j == 2 ) tgraph_embed_jets[j]->SetMarkerColorAlpha(plot_colors[j], 0.65);
            else tgraph_embed_jets[j]->SetMarkerColorAlpha(plot_colors[j], 1.0);
            tgraph_embed_jets[j]->Draw("p, same");
        }
        
        for ( int j = (c_jet_n - 1) ; j >= 0 ; j-- ) {
            char jet_pt[20];
            float  x_offset = 0.005; // Outline offset in plot axes units
            float  y_offset = 0.005; // Outline offset in plot axes units
            snprintf(jet_pt, 20, "%3.2f GeV", c_jet_pt[j]);
            jet_pt_outline[j] = new TText(c_jet_phi[j] + x_offset, c_jet_y[j] + y_offset, jet_pt);
            jet_pt_outline[j]->SetTextColor(kWhite);
            jet_pt_outline[j]->SetTextSize(0.04);
            jet_pt_outline[j]->SetTextAlign(22);
            jet_pt_text[j] = new TText(c_jet_phi[j], c_jet_y[j], jet_pt);
            jet_pt_text[j]->SetTextColor(kBlack);
            jet_pt_text[j]->SetTextSize(0.04);
            jet_pt_text[j]->SetTextAlign(22);
            jet_pt_outline[j]->Draw("same");
            jet_pt_outline[j]->DrawText(c_jet_phi[j] + x_offset, c_jet_y[j] + y_offset, "same");
            jet_pt_outline[j]->DrawText(c_jet_phi[j] + x_offset, c_jet_y[j] - y_offset, "same");
            jet_pt_outline[j]->DrawText(c_jet_phi[j] - x_offset, c_jet_y[j] + y_offset, "same");
            jet_pt_outline[j]->DrawText(c_jet_phi[j] - x_offset, c_jet_y[j] - y_offset, "same");
            jet_pt_text[j]->Draw("same");
        }
        
        char dir_tgraph[200];
        sprintf(dir_tgraph, "%stgraph_event%i_jets_comb.pdf", plot_dir, e);
        canvas->Print(dir_tgraph);
    }
    
    // TGraph of Jets
    for ( int e = first_event ; e < last_event ; e++ ) {
        pythJet_tree->GetEntry(e);

        TGraph* tgraph_pythia_jets[max_jets];
        TText*  jet_pt_text[max_jets];
        TText*  jet_pt_outline[max_jets];

        TCanvas* canvas = new TCanvas("canvas", "", 2000, 600);
        gPad->SetTicks();
        gStyle->SetOptStat(0);
        
        TH2D* th2d_jet_pt_distribution = new TH2D(Form("th2d_event%i_jet_pt_distribution", e), "Jet p_{T} and Constituent Distribution; #phi; #eta", 20, 0, 2*3.141593, 18, -0.9, 0.9);
        th2d_jet_pt_distribution->Draw("p, same");
        
        for ( int j = (p_jet_n - 1) ; j >= 0 ; j-- ) {
            float  particle_eta[max_parts], particle_phi[max_parts];

            for ( int p = 0 ; p < p_jet_part_n[j] ; p++ ) {
                particle_eta[p] = p_jet_part_eta[j][p];
                particle_phi[p] = p_jet_part_phi[j][p];
            }
            
            tgraph_pythia_jets[j] = new TGraph(p_jet_part_n[j], particle_phi, particle_eta);
            tgraph_pythia_jets[j]->SetMarkerStyle(20);
            if ( j == 0 ) tgraph_pythia_jets[j]->SetMarkerColorAlpha(plot_colors[j], 0.90);
            else if ( j == 1 ) tgraph_pythia_jets[j]->SetMarkerColorAlpha(plot_colors[j], 0.80);
            else if ( j == 2 ) tgraph_pythia_jets[j]->SetMarkerColorAlpha(plot_colors[j], 0.65);
            else tgraph_pythia_jets[j]->SetMarkerColorAlpha(plot_colors[j], 1.0);
            tgraph_pythia_jets[j]->Draw("p, same");
        }
        
        for ( int j = (p_jet_n - 1) ; j >= 0 ; j-- ) {
            char   jet_pt[20];
            float  x_offset = 0.005; // Outline offset in plot axes units
            float  y_offset = 0.005; // Outline offset in plot axes units
            snprintf(jet_pt, 20, "%3.2f GeV", p_jet_pt[j]);
            jet_pt_outline[j] = new TText(p_jet_phi[j] + x_offset, p_jet_y[j] + y_offset, jet_pt);
            jet_pt_outline[j]->SetTextColor(kWhite);
            jet_pt_outline[j]->SetTextSize(0.04);
            jet_pt_outline[j]->SetTextAlign(22);
            jet_pt_text[j] = new TText(p_jet_phi[j], p_jet_y[j], jet_pt);
            jet_pt_text[j]->SetTextColor(kBlack);
            jet_pt_text[j]->SetTextSize(0.04);
            jet_pt_text[j]->SetTextAlign(22);
            jet_pt_outline[j]->Draw("same");
            jet_pt_outline[j]->DrawText(p_jet_phi[j] + x_offset, p_jet_y[j] + y_offset, "same");
            jet_pt_outline[j]->DrawText(p_jet_phi[j] + x_offset, p_jet_y[j] - y_offset, "same");
            jet_pt_outline[j]->DrawText(p_jet_phi[j] - x_offset, p_jet_y[j] + y_offset, "same");
            jet_pt_outline[j]->DrawText(p_jet_phi[j] - x_offset, p_jet_y[j] - y_offset, "same");
            jet_pt_text[j]->Draw("same");
        }
        
        char dir_tgraph[200];
        sprintf(dir_tgraph, "%stgraph_event%i_jets_pyth.pdf", plot_dir, e);
        canvas->Print(dir_tgraph);
    }
}

double Func_ModifiedHagedorn(double* x, double* par) {
    double arg1 = par[0] * pow(x[0], 2) / pow( (pow(x[0], 2) + pow(par[1], 2)), 0.5);
    double arg2 = pow( 1 + (x[0] /par[2]), par[3] );
    double arg3 = 1 / x[0];
    double fitval =  arg1 * arg2 * arg3;
    return fitval;
}

void Plot_HEP_Data() {
    string hep_file_path = "../../External/HEP_Data/HEPData-ins879583-v1-root.root";
    string hep_hist_name = "Table 1/Hist1D_y1";
    string hep_hist_err1 = "Table 1/Hist1D_y1_e1";
    string hep_hist_err2 = "Table 1/Hist1D_y1_e2";
    string output_path = "../../Files/Thesis_Data/Plots_Thesis/Modified_Hagedorn_HEP_Data.pdf";
    vector<string> plot_title = {"Particle p_{T} for Pb+Pb at #sqrt{s_{NN}} = 2.76 TeV", "HEP Data Fit with Modified Hagedorn"};
    
//    TH1D* th1d_canvas_plot = new TH1D("canvas_plot", (plot_title[0] + ";" + plot_axes).c_str(), x_bins, x_min, x_max);
//        th1d_canvas_plot->SetAxisRange(y_min, y_max, "Y");
    TF1* thermal_pt_func  = new TF1("thermal_pt_func",  Func_ModifiedHagedorn, 0.5, 20., 4);
//    thermal_pt_func->SetParameters(64547, 3.076, 1.126, -8.491);
    thermal_pt_func->SetParameters(d_MH_par_1, d_MH_par_2, d_MH_par_3, d_MH_par_4);
    
    TFile* hep_file = new TFile(hep_file_path.c_str(), "READ");
    TH1* th1_hepdata_plot = (TH1*) hep_file->Get(hep_hist_name.c_str());
    TH1* th1_hepdata_err1 = (TH1*) hep_file->Get(hep_hist_err1.c_str());
    TH1* th1_hepdata_err2 = (TH1*) hep_file->Get(hep_hist_err2.c_str());
    
    for ( int b = 0 ; b < th1_hepdata_plot->GetNbinsX() ; b++ ) {
        double data_val = th1_hepdata_plot->GetBinContent(b);
        double stat_err = th1_hepdata_err1->GetBinContent(b);
        double syst_err = th1_hepdata_err2->GetBinContent(b);
        double total_err = sqrt( stat_err*stat_err + syst_err*syst_err);
        if (total_err / data_val < 0.4) total_err = 0.4 * data_val;
        th1_hepdata_plot->SetBinError( b + 1 , total_err );
    }
    
    // TCanvas("name", "title", width (px), height (px))
    TCanvas* canvas = new TCanvas("canvas", "", 1000, 600);
    gPad->SetTicks();
    gStyle->SetOptStat(0);
    gPad->SetLogy(1);
    
    th1_hepdata_plot->SetLineColor(plot_black);
    th1_hepdata_plot->SetMarkerColor(plot_black);
    th1_hepdata_plot->SetMarkerStyle(mark_circ_open[0]);
    th1_hepdata_plot->SetMarkerSize(mark_circ_open[1]);
    
    thermal_pt_func->SetLineColor(plot_red);
    thermal_pt_func->SetLineWidth(2);

    // Turns on Latex formatting relative to canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(kTRUE);
    
    th1_hepdata_plot->Draw();
    thermal_pt_func->Draw("same");
    
    TLegend* legend = new TLegend(0.60, 0.75, 0.85, 0.85);
    legend->AddEntry(th1_hepdata_plot, "Particle p_{T} for Pb+Pb", "lp");
    legend->AddEntry(thermal_pt_func, "Modified Hagedorn Fit", "lp");
    legend->SetLineWidth(0);
    legend->SetFillStyle(0);
    legend->Draw("same");
    
    for ( int i=0 ; i < plot_title.size() ; i++ ) {
        latex->DrawLatex(0.15, 0.81 - (0.06*i), ("#scale[0.75]{" + plot_title[i] + "}").c_str());
    }
    
    canvas->Print(output_path.c_str());
    
    gPad->SetLogy(0);
    
    delete legend;
    delete latex;
    delete canvas;
    std::cout << "Plotting complete!" << std::endl;
}




void MACRO_Thesis_Plots() {
    
//    // Toy Model - HEP Data Comparison
//    Plot_HEP_Data();

    // Jet Truth Distributions
    Plot_Simple(
        {   {"../../Files/Thesis_Data/Full_Train_B8_10_90_N500000_ML_Prep.root", "ML_Train_B8_10_90_N500000", "jet_pt_true", "p_{T}^{8} Bias Dist.", "plot_blue", "circ_open", "false"},
            {"../../Files/Thesis_Data/Full_Train_B8_10_90_N500000_ML_Prep.root", "ML_Train_B8_10_90_N500000_Flat", "jet_pt_true", "Flattened Dist.", "plot_red", "circ_open", "false"}
        }, // file_name(.root), tree_name, variable, plot_legend_label, plot_color, plot_symbol, use_var_int
        {"Flat Distribution (from p_{T}^{8} Bias)", "p_{T, Jet}^{True} in 10-90 GeV"},
        "p_{T, Jet}^{True} [GeV]; dN_{Jets}/dp_{T}",
        "../../Files/Thesis_Data/Plots_Thesis/ptTrue_Distribution_Flat.pdf",
        50, // x_bins (int)
        0., // x_min
        100., // x_max
        true, // y_log (bool)
        1., // y_min
        10000000., // y_max
        false, // use_normalized (bool)
        false, // show_widths
        700, // plot_width
        600 // plot_height
    );
    Plot_Simple(
        {
            {"../../Files/Thesis_Data/Full_Train_B4_10_90_N500000_ML_Prep.root", "ML_Train_B4_10_90_N500000", "jet_pt_true", "p_{T}^{4} Bias Dist.", "plot_blue", "circ_open", "false"}
        }, // file_name(.root), tree_name, variable, plot_legend_label, plot_color, plot_symbol, use_var_int
        {"p_{T}^{4} Bias Distribution", "p_{T, Jet}^{True} in 10-90 GeV"},
        "p_{T, Jet}^{True} [GeV]; dN_{Jets}/dp_{T}",
        "../../Files/Thesis_Data/Plots_Thesis/ptTrue_Distribution_B4.pdf",
        50, // x_bins (int)
        0., // x_min
        100., // x_max
        true, // y_log (bool)
        1., // y_min
        10000000., // y_max
        false, // use_normalized (bool)
        false, // show_widths
        700, // plot_width
        600 // plot_height
    );
    Plot_Simple(
        {
            {"../../Files/Thesis_Data/Full_Train_B0_10_90_N500000_ML_Prep.root", "ML_Train_B0_10_90_N500000", "jet_pt_true", "No Bias Dist.", "plot_blue", "circ_open", "false"}
        }, // file_name(.root), tree_name, variable, plot_legend_label, plot_color, plot_symbol, use_var_int
        {"No p_{T}-Bias Distribution", "p_{T, Jet}^{True} in 10-90 GeV"},
        "p_{T, Jet}^{True} [GeV]; dN_{Jets}/dp_{T}",
        "../../Files/Thesis_Data/Plots_Thesis/ptTrue_Distribution_B0.pdf",
        50, // x_bins (int)
        0., // x_min
        100., // x_max
        true, // y_log (bool)
        1., // y_min
        10000000., // y_max
        false, // use_normalized (bool)
        false, // show_widths
        700, // plot_width
        600 // plot_height
    );


//    Event_2D_Plotter(
//        "../../Files/Thesis_Data/Generator_Data.nosync/Full_Train_T2_B8_10_90_N500000.root", // comb_file_path
//        "Combined_Jet_Tree", // comb_tree_name
//        "../../Files/Thesis_Data/Generator_Data.nosync/Full_Train_T2_B8_10_90_N500000.root", // pyth_file_path
//        "PYTHIA_Jet_Tree", // pyth_tree_name
//        "../../Files/Thesis_Data/Plots_Thesis/", // output_dir
//        20, // first_event
//        40 // last_event
//    );
    
    
    // Area-Correction Comparison
    Plot_AreaCorr_Compare(
        {{"../../Files/Thesis_Data/Train_B8_Flat/Test_B8_Flat/Root/Train_B8_Flat_F11_10_90_Test_Centered_Wide_Bins.root", "Test_40_60_LR", "jet_pt_corr", "p_{T, Jet}^{Area Corr.} (Thesis)", "plot_red", "circ_open", "false"}
        }, // file_name(.root), tree_name, variable, plot_legend_label, plot_color, plot_symbol, variable_is_integer("true"/"false")
        {"Flat Distribution for p_{T, Jet}^{True} in 40-60 GeV", "Area-Correction Comparison"},
        "p_{T, Jet}^{Corr.} - p_{T, Jet}^{True} [GeV]; 1/N dN_{Jets}/dp_{T} [GeV^{-1}]",
        "../../Files/Thesis_Data/Plots_Thesis/B8_Flat_40_60_areacorr_comparison.pdf",
        50, // x_bins (int)
        -40., // x_min
        40., // x_max
        false, // y_log (bool)
        0., // y_min
        .08, // y_max
        true, // use_normalized (bool)
        true // show_widths
    );


    // Constituent pT
    Plot_Simple(
        {   {"../../Files/Thesis_Data/Full_Train_B8_10_90_N500000_ML_Prep.root", "ML_Train_B8_10_90_N500000_Flat", "const_4_pt", "p_{T, Const. 4}", "plot_green", "circ_open", "false"},
            {"../../Files/Thesis_Data/Full_Train_B8_10_90_N500000_ML_Prep.root", "ML_Train_B8_10_90_N500000_Flat", "const_3_pt", "p_{T, Const. 3}", "plot_blue", "circ_open", "false"},
            {"../../Files/Thesis_Data/Full_Train_B8_10_90_N500000_ML_Prep.root", "ML_Train_B8_10_90_N500000_Flat", "const_2_pt", "p_{T, Const. 2}", "plot_violet", "circ_open", "false"},
            {"../../Files/Thesis_Data/Full_Train_B8_10_90_N500000_ML_Prep.root", "ML_Train_B8_10_90_N500000_Flat", "const_1_pt", "p_{T, Const. 1}", "plot_red", "circ_open", "false"}
        }, // file_name(.root), tree_name, variable, plot_legend_label, plot_color, plot_symbol, use_var_int
        {"Flat Distribution for p_{T, Jet}^{True} in 10-90 GeV", "p_{T} for Constituents 1, 2, 3, 4"},
        "p_{T, Const.} [GeV]; dN_{Jets}/dp_{T}",
        "../../Files/Thesis_Data/Plots_Thesis/B8_Flat_10_90_const_pt_1234.pdf",
        60, // x_bins (int)
        0., // x_min
        40., // x_max
        false, // y_log (bool)
        0., // y_min
        70000., // y_max
        false, // use_normalized (bool)
        false // show_widths
    );

    // Jet Const. Mean, Median
    Plot_Simple(
        {   {"../../Files/Thesis_Data/Full_Train_B8_10_90_N500000_ML_Prep.root", "ML_Train_B8_10_90_N500000_Flat", "const_pt_mean", "Mean p_{T, Const.}", "plot_red", "circ_open", "false"},
            {"../../Files/Thesis_Data/Full_Train_B8_10_90_N500000_ML_Prep.root", "ML_Train_B8_10_90_N500000_Flat", "const_pt_median", "Median p_{T, Const.}", "plot_blue", "circ_open", "false"}
        }, // file_name(.root), tree_name, variable, plot_legend_label, plot_color, plot_symbol
        {"Flat Distribution for p_{T, Jet}^{True} in 10-90 GeV", "Constituent p_{T} Mean and Median"},
        "p_{T, Const.} [GeV]; dN_{Jets}/dp_{T}",
        "../../Files/Thesis_Data/Plots_Thesis/B8_Flat_10_90_const_pt_mean_median.pdf",
        60, // x_bins (int)
        0., // x_min
        3., // x_max
        false, // y_log (bool)
        0., // y_min
        100000., // y_max
        false, // use_normalized (bool)
        false // show_widths
    );

    // Jet Area
    Plot_Simple(
        {{"../../Files/Thesis_Data/Full_Train_B8_10_90_N500000_ML_Prep.root", "ML_Train_B8_10_90_N500000_Flat", "jet_area", "Jet Area", "plot_blue", "circ_open", "false"}
        }, // file_name(.root), tree_name, variable, plot_legend_label, plot_color, plot_symbol
        {"Flat Distribution for p_{T, Jet}^{True} in 10-90 GeV", "Jet Area"},
        "A_{Jet}; dN_{Jets}/dA",
        "../../Files/Thesis_Data/Plots_Thesis/B8_Flat_10_90_jet_area.pdf",
        40, // x_bins (int)
        0., // x_min
        0.8, // x_max
        false, // y_log (bool)
        0., // y_min
        90000, // y_max
        false, // use_normalized (bool)
        false // show_widths
    );

    // Jet Background Density
    Plot_Simple(
        {{"../../Files/Thesis_Data/Full_Train_B8_10_90_N500000_ML_Prep.root", "ML_Train_B8_10_90_N500000_Flat", "jet_rho", "Background p_{T} Density", "plot_blue", "circ_open", "false"}
        }, // file_name(.root), tree_name, variable, plot_legend_label, plot_color, plot_symbol
        {"Flat Distribution for p_{T, Jet}^{True} in 10-90 GeV", "Background p_{T} Density"},
        "#rho_{Jet Event} [GeV/A^{2}]; dN_{Jets}/d#rho",
        "../../Files/Thesis_Data/Plots_Thesis/B8_Flat_10_90_jet_bgrho.pdf",
        60, // x_bins (int)
        45., // x_min
        175, // x_max
        false, // y_log (bool)
        0., // y_min
        30000, // y_max
        false, // use_normalized (bool)
        false // show_widths
    );

    // Jet Const. N
    Plot_Simple(
        {{"../../Files/Thesis_Data/Full_Train_B8_10_90_N500000_ML_Prep.root", "ML_Train_B8_10_90_N500000_Flat", "jet_const_n", "Jet Constituents", "plot_blue", "circ_open", "true"}
        }, // file_name(.root), tree_name, variable, plot_legend_label, plot_color, plot_symbol
        {"Flat Distribution for p_{T, Jet}^{True} in 10-90 GeV", "Jet Constituents"},
        "N_{Const.}; dN_{Jets}/dN_{Const.}",
        "../../Files/Thesis_Data/Plots_Thesis/B8_Flat_10_90_jet_const_n.pdf",
        60, // x_bins (int)
        0., // x_min
        180., // x_max
        false, // y_log (bool)
        0., // y_min
        30000., // y_max
        false, // use_normalized (bool)
        false // show_widths
    );

    // Jet Mass
    Plot_Simple(
        {{"../../Files/Thesis_Data/Full_Train_B8_10_90_N500000_ML_Prep.root", "ML_Train_B8_10_90_N500000_Flat", "jet_mass", "Jet Mass", "plot_blue", "circ_open", "false"}
        }, // file_name(.root), tree_name, variable, plot_legend_label, plot_color, plot_symbol
        {"Flat Distribution for p_{T, Jet}^{True} in 10-90 GeV", "Jet Mass"},
        "m_{Jet} [MeV/c^{2}]; dN_{Jets}/dm",
        "../../Files/Thesis_Data/Plots_Thesis/B8_Flat_10_90_jet_mass.pdf",
        60, // x_bins (int)
        0., // x_min
        60., // x_max
        false, // y_log (bool)
        0., // y_min
        30000, // y_max
        false, // use_normalized (bool)
        false // show_widths
    );

    // Jet pT True
    Plot_Simple(
        {{"../../Files/Thesis_Data/Full_Train_B8_10_90_N500000_ML_Prep.root", "ML_Train_B8_10_90_N500000_Flat", "jet_pt_true", "p_{T, Jet}^{True}", "plot_blue", "circ_open", "false"}
        }, // file_name(.root), tree_name, variable, plot_legend_label, plot_color, plot_symbol
        {"Flat Distribution for p_{T, Jet}^{True} in 10-90 GeV", "Jet p_{T}^{True}"},
        "p_{T, Jet} [GeV]; dN_{Jets}/dp_{T}",
        "../../Files/Thesis_Data/Plots_Thesis/B8_Flat_10_90_jet_pt_true.pdf",
        100, // x_bins (int)
        0., // x_min
        100., // x_max
        true, // y_log (bool)
        10., // y_min
        100000., // y_max
        false, // use_normalized (bool)
        false // show_widths
    );

    // Jet pT Raw
    Plot_Simple(
        {{"../../Files/Thesis_Data/Full_Train_B8_10_90_N500000_ML_Prep.root", "ML_Train_B8_10_90_N500000_Flat", "jet_pt_raw", "p_{T, Jet}^{Raw}", "plot_blue", "circ_open", "false"}
        }, // file_name(.root), tree_name, variable, plot_legend_label, plot_color, plot_symbol
        {"Flat Distribution for p_{T, Jet}^{True} in 10-90 GeV", "Jet p_{T}^{Raw}"},
        "p_{T, Jet} [GeV]; dN_{Jets}/dp_{T}",
        "../../Files/Thesis_Data/Plots_Thesis/B8_Flat_10_90_jet_pt_raw.pdf",
        100, // x_bins (int)
        0., // x_min
        200., // x_max
        true, // y_log (bool)
        10., // y_min
        100000., // y_max
        false, // use_normalized (bool)
        false // show_widths
    );

    // Jet Rapidity
    Plot_Simple(
        {{"../../Files/Thesis_Data/Full_Train_B8_10_90_N500000_ML_Prep.root", "ML_Train_B8_10_90_N500000_Flat", "jet_y", "Jet Rapidity", "plot_blue", "circ_open", "false"}
        }, // file_name(.root), tree_name, variable, plot_legend_label, plot_color, plot_symbol
        {"Flat Distribution for p_{T, Jet}^{True} in 10-90 GeV", "Jet Rapidity"},
        "y_{Jet}; dN_{Jets}/dy",
        "../../Files/Thesis_Data/Plots_Thesis/B8_Flat_10_90_jet_y.pdf",
        90, // x_bins (int)
        -0.9, // x_min
        0.9, // x_max
        false, // y_log (bool)
        0., // y_min
        10000, // y_max
        false, // use_normalized (bool)
        false // show_widths
    );

}
