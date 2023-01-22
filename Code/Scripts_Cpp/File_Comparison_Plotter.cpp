#include "TFile.h"
#include "TTree.h"

#include "jet_ml_constants.h"
using namespace Jet_ML_Constants;

namespace Comparison_Plotter {

// Script for plotting two files at once
// need to add a macro file template

void File_Comparison_TH1 (
    char  file_A_path[400],             // Full file path to file A
    char  tree_A_name[100],             // Tree to use for comparison for file A
    char  file_A_label[100],            // Plot legend label for file A
    char  tree_A_branches[30][50],      // Array of branches of which to print comparisons
    char  file_B_path[400],             // Full file path to file B
    char  tree_B_name[100],             // Tree to use for comparison for file B
    char  file_B_label[100],            // Plot legend label for file B
    char  tree_B_branches[30][50],     // Array of branches of which to print comparisons
    char  output_directory[400],        //
    const int num_bins,
    int   branch_bins[30],        //
    float branch_range[30][2],    //
    char  branch_x_units[30][50]  //
    ) {
    
    TFile* file_A = new TFile(file_A_path, "READ");
    TTree* tree_A = (TTree*) file_A->Get(tree_A_name);
    std::cout << "File A and Tree A accessed: " << file_A_path << ", " << tree_A_name << std::endl;
    
    TFile* file_B = new TFile(file_B_path, "READ");
    TTree* tree_B = (TTree*) file_B->Get(tree_B_name);
    std::cout << "File B and Tree B accessed: " << file_B_path << ", " << tree_B_name << std::endl;
    
    // Loops through branch names to compare and plots comparisons.
    for ( int i=0 ; i < num_bins ; i++ ) {
        if ( !branch_bins ) break;
        
        char branch_A[100];
        snprintf(branch_A, 100, "%s", tree_A_branches[i]);
        char branch_B[100];
        snprintf(branch_B, 100, "%s", tree_B_branches[i]);
        std::cout << "Plotting Tree A: " << branch_A << ", and Tree B: " << branch_B << std::endl;

        int   n_bins = branch_bins[i];
        float x_min  = branch_range[i][0];
        float x_max  = branch_range[i][1];
        char plot_labels[200];
        snprintf(plot_labels, 200, "Comparison of %s; %s; Probability Density", branch_A, branch_x_units[i]);
        
        gStyle->SetOptStat(0);
        
        std::cout << "    Filling histograms..." << std::endl;
        TH1F* canvas_plot = new TH1F("canvas_plot", plot_labels, n_bins, x_min, x_max);
        TH1F* plot_A = new TH1F(branch_A, branch_A, n_bins, x_min, x_max);
        char tree_A_draw_command[100];
        snprintf(tree_A_draw_command, 100, "%s>>%s", branch_A, branch_A);
        tree_A->Draw(tree_A_draw_command);
        TH1F* plot_B = new TH1F(branch_B, branch_B, n_bins, x_min, x_max);
        char tree_B_draw_command[100];
        snprintf(tree_B_draw_command, 100, "%s>>%s", branch_B, branch_B);
        tree_B->Draw(tree_B_draw_command);
        
        plot_A->Scale( 1. / plot_A->Integral(),"WIDTH");
        plot_B->Scale( 1. / plot_B->Integral(),"WIDTH");
        float y_max = 1;
        if ( plot_A->GetMaximum() > plot_B->GetMaximum() ) y_max = 1.2 * plot_A->GetMaximum();
        else y_max = 1.2 * plot_B->GetMaximum();
        canvas_plot->SetMaximum(y_max);
        
        std::cout << "    Histograms filled. Plotting..." << std::endl;
        
        TCanvas* canvas = new TCanvas("canvas", "", 1000, 600);
        gPad->SetTicks();
        gPad->SetLogy(0);
        TLatex* latex = new TLatex();
        latex->SetNDC(kTRUE);
        
        plot_A->SetLineColor(plot_red);
        plot_A->SetMarkerColor(plot_red);
        plot_A->SetMarkerStyle(mark_circ_open[0]);
        plot_A->SetMarkerSize(mark_circ_open[1]);
        
        plot_B->SetLineColor(plot_blue);
        plot_B->SetMarkerColor(plot_blue);
        plot_B->SetMarkerStyle(mark_diam_open[0]);
        plot_B->SetMarkerSize(mark_diam_open[1]);
        
        TLegend* legend = new TLegend(0.65, 0.75, 0.85, 0.85);
        legend->SetLineWidth(0);
        legend->SetFillStyle(0);
        legend->AddEntry(plot_A, file_A_label, "lp");
        legend->AddEntry(plot_B, file_B_label, "lp");
        
        canvas_plot->Draw();
        plot_A->Draw("same");
        plot_B->Draw("same");
        legend->Draw("same");
        
        char plot_file_name[400];
        snprintf(plot_file_name, 400, "%s/hist_compare_%s.pdf", output_directory, branch_A);
        
        std::cout << "    Histogram plotted!" << std::endl;
        
        canvas->Print(plot_file_name);
        delete plot_A;
        delete plot_B;
        delete canvas_plot;
        delete legend;
        delete canvas;
    } // end of for loop
    
    file_A->Close();
    file_B->Close();
}
    
} // End of namespace
