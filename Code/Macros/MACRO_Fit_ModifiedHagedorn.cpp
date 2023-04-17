#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TCanvas.h"
#include <cmath>
#include <iostream>
#include <filesystem>
#include <string>


#include "../Scripts_Cpp/jet_ml_constants.h" // Contains default values for many parameters.
using namespace Jet_ML_Constants;

const bool printOut = true;
const bool debug = false;

/*
    THIS FUNCTION IS ~ONLY~ FOR FITTING TO HEP DATA:
        doi:10.1016/j.physletb.2010.12.020
    IT INCLUDES A 1/x TERM TO FIT TO THE HISTOGRAM!
    If you fit data without the other function, it will be wrong!
    If you generate events with this function, it will be wrong!
*/
double ModifiedHagedorn_Func_Fit(
    double* x,
    double* par
    ) {
    double arg1 = par[0] * pow(x[0], 2) / pow( (pow(x[0], 2) + pow(par[1], 2)), 0.5);
    double arg2 = pow( 1 + (x[0] /par[2]), par[3] );
    double arg3 = 1 / x[0];
    double fitval = arg1 * arg2 * arg3;
    return fitval;
}

/*
    THIS FUNCTION IS ~ONLY~ FOR EVENT GENERATION!
 */
double ModifiedHagedorn_Func(
    double* x,
    double* par
    ) {
    double arg1 = par[0] * pow(x[0], 2) / pow( (pow(x[0], 2) + pow(par[1], 2)), 0.5);
    double arg2 = pow( 1 + (x[0] /par[2]), par[3] );
    double fitval = arg1 * arg2;
    return fitval;
}

double* Fit_ModifiedHagedorn(
    string input_file_path,
    string input_hist,
    string input_hist_stat_err,
    string input_hist_syst_err,
    string output_directory
    ) {
    
    TFile* input_file = new TFile(input_file_path.c_str(),"READ");
    
    TH1F* th1f_PbPb_data    = (TH1F*) input_file->Get(input_hist.c_str());
    TH1F* th1f_PbPb_stat_err = (TH1F*) input_file->Get(input_hist_stat_err.c_str());
    TH1F* th1f_PbPb_syst_err  = (TH1F*) input_file->Get(input_hist_syst_err.c_str());
    
    std::cout << "File and hists accessed" << std::endl;

    // Morgan K Code:
    // Adds error to each bin
    for ( int b = 0 ; b < th1f_PbPb_data->GetNbinsX() ; b++ ) {
        double data_val = th1f_PbPb_data->GetBinContent(b);
        double stat_err = th1f_PbPb_stat_err->GetBinContent(b);
        double syst_err = th1f_PbPb_syst_err->GetBinContent(b);
        double total_err = sqrt( stat_err*stat_err + syst_err*syst_err);
        if (total_err / data_val < 0.4) total_err = 0.4 * data_val;
        th1f_PbPb_data->SetBinError( b + 1 , total_err );
    }

    // Declares TF1
    // ("Title", name_of_fitting_function, min_fit_bound, max_fit_bound, number_of_fit_function_parameters)
    TF1* tf1_PbPb_pt_func = new TF1("tf1_PbPb_pt_func", ModifiedHagedorn_Func_Fit, 0.5, 20, 4);
    
    // Set values as initial guesses for fit function
    tf1_PbPb_pt_func->SetParameters(45400., m_pion, .7, -7);

    // Replace with HEP data histogram, (local tf1 with set params, fit args)
    // "R" = fit only within set range
    // "N" = tells plot funcs not to draw fit func
    // "L" = log liklihood fitting good for counts
    th1f_PbPb_data->Fit(tf1_PbPb_pt_func,"RNL");
    
    // TCanvas("name", "title", width (px), height (px))
    TCanvas* canvas = new TCanvas("canvas", "", 1000, 600);
    
    // Turns on ticks, turns off stats
    gPad->SetTicks();
    gStyle->SetOptStat(0);
    gPad->SetLogy(0);
    
    // Turns on Latex formatting relative to canvas
    TLatex* latex = new TLatex();
    latex->SetNDC(kTRUE);
    
    // Build histogram plot
    th1f_PbPb_data->SetLineColor(plot_black);
    th1f_PbPb_data->SetMarkerColor(plot_black);
    th1f_PbPb_data->SetMarkerStyle(mark_circ_open[0]);
    th1f_PbPb_data->SetMarkerSize(mark_circ_open[1]);
    
    tf1_PbPb_pt_func->SetLineColor(plot_red);
    
    // Print plot with log y-axis
    gPad->SetLogy(1);
    th1f_PbPb_data->Draw();
    tf1_PbPb_pt_func->Draw("same");
    canvas->Print((output_directory + "PbPb_pt_distribution_log.pdf").c_str());
    std::cout << "Plotted plot." << std::endl;
    
    // Print plot with linear y-axis
    gPad->SetLogy(0);
    th1f_PbPb_data->Draw();
    tf1_PbPb_pt_func->Draw("same");
    canvas->Print((output_directory + "PbPb_pt_distribution.pdf").c_str());
    std::cout << "Plotted plot." << std::endl;
    
    // Accessing fit param values
    double par_arr[4];
    
    for (int i=0 ; i < 4 ; i++ ) par_arr[i] = double(tf1_PbPb_pt_func->GetParameter(i));
    
    std::cout << par_arr << std::endl;
    
    input_file->Close();
    
    return par_arr;
}

void MACRO_Fit_ModifiedHagedorn() {
    string input_file_path = "../../External/HEP_Data/HEPData-ins879583-v1-root.root";
    string input_hist = "Table 1/Hist1D_y1";
    string input_hist_stat_err = "Table 1/Hist1D_y1_e1";
    string input_hist_syst_err = "Table 1/Hist1D_y1_e2";
    string output_directory = "../../Files/Thesis_Data/Plots_Thesis/";
    
    Fit_ModifiedHagedorn(
        input_file_path,
        input_hist,
        input_hist_stat_err,
        input_hist_syst_err,
        output_directory
    );
}
