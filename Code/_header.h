#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TPDF.h"
#include <stdio.h>
#include <cmath>

//#if !defined(MYLIB_ML_PROJECT_HEADER_H)
//#define MYLIB_ML_PROJECT_HEADER_H

//  File Name : _header.h
namespace Project_Constants {

// General Settings
const int    trialNum       = 2;
const double math_pi        = 3.14159265359;
const double math_e         = 2.71828182846;
const double m_pion         = 0.1396;           // Pion+/- mass in GeV

// PYTHIA Settings
const float  beamPower      = 2760.;            // Beam eCM, [GeV]
const float  ptJetMin       = 20.;              // Minimum jet pT, [GeV]
const float  ptHatMin       = 20.;              // Total pT minimum, [GeV]
const float  ptBiasPow      = 1.;                // pT is biased by this power
const bool   slimJets       = true;             // If true, slims jets to satsify conditions
const double slim_pt_min    = 40.;              // Min E for jet slimming, [GeV]
const double slim_pt_max    = 60.;              // Max E for jet slimming, [GeV]
const double slim_rap       = 0.5;              // Min/Max rapidity for jet slimming

// Thermal Settings
const int    gaus_mean      = 1800;
const int    gaus_sigma     = 200;
const double gaus_norm      = 1 / sqrt(2 * math_pi * pow(gaus_sigma,2));
const double power_exp      = -3.1;
const double power_zero     = -1.;
const double eta_max        = 0.9;

// FastJet Settings
const double fj_jetR        = 0.4;
const double fj_jetptmin    = 5.0;
const double fj_rSquared    = pow(0.3, 2);
// Set input_source to "particle", "thermal", or "combined"
const char* input_source = "combined";

// General Directories
char  dir_master[200];
int n1 = sprintf(dir_master, "../Files/ALICE_kinematics_trial%i", trialNum);
//int n1 = sprintf(dir_master, "../_Files/N-%i_beam-%i_ptMin-%i_ptHatMin-%i_ptBias-%i_gaus-%i_sigma-%i_modHag_trial%i",
//    nEvent, int(beamPower), int(ptJetMin), int(ptHatMin), int(ptBiasPow), gaus_mean, gaus_sigma, 1);

char  dir_data[200];
int n2 =  sprintf(dir_data, "%s/Data", dir_master);

char  dir_plots[200];
int n3 = sprintf(dir_plots, "%s/Plots", dir_master);

// PYTHIA Directories
const char* file_p_trees    = Form("%s/PYTHIA_Trees.root", dir_data);
const char* file_p_plots    = Form("%s/PYTHIA_Plots.root", dir_data);

// Thermal Directories
const char* file_t_trees    = Form("%s/Thermal_Trees.root", dir_data);
const char* file_t_plots    = Form("%s/Thermal_Plots.root", dir_data);

// Combined Directories
const char* file_c_trees    = Form("%s/Combined_Trees.root", dir_data);
const char* file_c_plots    = Form("%s/Combined_Plots.root", dir_data);

char file_fj_p_trees[200];
int n4 = sprintf(file_fj_p_trees, "%s/FastJet_PYTHIA_Trees.root", dir_data);
char file_fj_t_trees[200];
int n5 = sprintf(file_fj_t_trees, "%s/FastJet_Thermal_Trees.root", dir_data);
char file_fj_c_trees[200];
int n6 = sprintf(file_fj_c_trees, "%s/FastJet_Combined_Trees.root", dir_data);

char file_fj_p_plots[200];
int n7 = sprintf(file_fj_p_plots, "%s/FastJet_PYTHIA_Plots.root", dir_data);
char file_fj_t_plots[200];
int n8 = sprintf(file_fj_t_plots, "%s/FastJet_Thermal_Plots.root", dir_data);
char file_fj_c_plots[200];
int n9 = sprintf(file_fj_c_plots, "%s/FastJet_Combined_Plots.root", dir_data);

char file_fj_pc_plots[200];
int n10 = sprintf(file_fj_pc_plots, "%s/FastJet_Combined_PYTHIA_Plots.root", dir_data);

const char* file_fj_jet_root = Form("%s/Jet_Subtract_ML_TreesPlots.root", dir_data);

// Styles
const int particle_line_color   = kOrange+7;
const int particle_marker_color = kOrange+7;
const int particle_marker_style = 20;

const int jet_line_color        = kAzure-3;
const int jet_marker_color      = kAzure-3;
const int jet_marker_style      = 33;

const int thermal_line_color    = kPink+7;
const int thermal_marker_color  = kPink+7;
const int thermal_marker_style  = 20;

const int combined_line_color   = kTeal+2;
const int combined_marker_color = kTeal+2;
const int combined_marker_style = 20;

// Marker Styles
const double mark_circ_open[2] = {24, 1.0};
const double mark_circ_fill[2] = {20, 1.0};
const double mark_squa_open[2] = {25, 0.9};
const double mark_squa_fill[2] = {21, 0.9};
const double mark_diam_open[2] = {27, 1.5};
const double mark_diam_fill[2] = {33, 1.5};
const double mark_star_open[2] = {42, 1.4};
const double mark_star_fill[2] = {43, 1.4};
const double mark_plus_open[2] = {28, 0.9};
const double mark_plus_fill[2] = {34, 0.9};
const double mark_triu_open[2] = {26, 1.1};
const double mark_triu_fill[2] = {22, 1.1};
const double mark_trid_open[2] = {32, 1.1};
const double mark_trid_fill[2] = {23, 1.1};
// Jet Plotter Colors
int jet_gen_mark = kViolet+7;
int jet_gen_line = kViolet+7;
int jet_blk_mark = kGray+3;
int jet_blk_line = kGray+3;
int jet_red_mark = kPink;
int jet_red_line = kPink;
int jet_vio_mark = kViolet-1;
int jet_vio_line = kViolet-1;
int jet_blu_mark = kAzure;
int jet_blu_line = kAzure;
int jet_tea_mark = kTeal-6;
int jet_tea_line = kTeal-6;

void th1d_plotter(
    TH1D* histogram,
    char* plot_filename,
    int line_color,
    int marker_color,
    int marker_style,
    float marker_size,
    float ymin,
    bool logy,
    char* label_arr[10],
    int label_arr_size,
    float label_x,
    float label_y) {
    
    // TCanvas("name", "title", width (px), height (px))
    TCanvas* canvas = new TCanvas("canvas", "", 800, 600);
    
    // Turns on ticks, turns off stats
    gPad->SetTicks();
    gStyle->SetOptStat(0);
    gPad->SetLogy(0);
    
    // Turns on Latex formatting relative to canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(kTRUE);
    
    // Build histogram plot
    histogram->SetLineColor(line_color);
    histogram->SetMarkerColor(marker_color);
    histogram->SetMarkerStyle(int(marker_style));
    histogram->SetMarkerSize(marker_size);
    
    if (logy) gPad->SetLogy(1);
    histogram->SetMinimum(ymin);
//    else histogram->GetYaxis()->SetRangeUser(0, ymax);
    
    histogram->Draw();
    
    for ( int i = 0 ; i < label_arr_size ; i ++) {
        float label_ymax = ( 0.04 * label_arr_size ) /2. ;
        latex->DrawLatex(label_x, (label_y + label_ymax) - (i * 0.04), Form("#scale[0.6]{#bf{%s}}", label_arr[i]));
        std::cout << label_x << " , " << (label_y + label_ymax) - (i * 0.04) << std::endl;
    }
    
    canvas->Print(Form("%s/%s", dir_plots, plot_filename)); // Resets the canvas after printing
    std::cout << "Plotted " << plot_filename << std::endl;
    
    gPad->SetLogy(0);
    
    delete canvas;
    delete histogram;
}

} // End Project_Constants namespace
