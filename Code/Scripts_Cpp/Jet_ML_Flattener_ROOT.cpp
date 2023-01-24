#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

#include "jet_ml_constants.h" // Contains default values for many parameters.
using namespace Jet_ML_Constants;

namespace Jet_Flattener {

/*
    Make_Flat_Jet_Distribution()
    --- WORK IN PROGRESS ---
    Should flatten a TTree based on the lowest number of jets in a pT range.
    Need to modify to use random generator to throw out jets
 */
void Make_Flat_Distribution(
    char  mlprep_file_path[500],    // Full name of ML_Prep root file
    char  input_tree_name[100],     // Name of the TTree to be flattened
    char  pt_true_label[100],       // Label for jet_pt_true in the given TTree
    const int n_bins,               // number of bins for flattening
    const float jet_pt_min,         // [GeV] min
    const float jet_pt_max,         // [GeV] max
    int   min_bin_n                 // If the min bin count is below this, uses this value instead
    ) {
    
    TFile* file = new TFile(mlprep_file_path, "UPDATE");
    TTree* input_tree = (TTree*) file->Get(input_tree_name);
    
    float jet_pt_true;
    input_tree->SetBranchAddress(pt_true_label, &jet_pt_true);

    std::cout << "Input file accessed. pt_true branch accessed." << std::endl;
    
    char output_tree_name[100];
    snprintf(output_tree_name, 100, "%s_Flat", input_tree_name);
    char output_tree_desc[100];
    snprintf(output_tree_desc, 100, "Flattened version of %s", input_tree_name);
    TTree* output_tree = input_tree->CloneTree(0);
    output_tree->SetObject(output_tree_name, output_tree_desc);
    std::cout << "Output tree built." << std::endl;
    
    // Makes a histogram of the pt_true bins
    TH1* input_tree_pt = new TH1F("input_tree_pt", "input_tree_pt", n_bins, jet_pt_min, jet_pt_max);
    char tree_draw_command[100];
    snprintf(tree_draw_command, 100, "%s>>input_tree_pt", pt_true_label);
    input_tree->Draw(tree_draw_command);
    
    int min_bin_n_tree = input_tree_pt->GetMinimum();
    if ( min_bin_n < min_bin_n_tree ) min_bin_n = min_bin_n_tree;
    std::cout << "Minimum bin size is: " << min_bin_n << std::endl;
    
    // Sets up arrays to manage bin tracking and counting
    int   bin_n_array[n_bins];                    // Array of the number of jets in each bin
    float bin_keep_ratio[n_bins];                 // Array of the ratio of jets to keep from each bin
    float bin_upper_edge[n_bins];                 // Array of the upper edge of each bin
    float bin_width = input_tree_pt->GetBinWidth(0);    // Sets the width of each bin
    for ( int i = 0 ; i < n_bins ; i++ ) {
        int bin_size = input_tree_pt->GetBinContent(i+1);
        bin_n_array[i] = bin_size;
        float keep_ratio = float(min_bin_n) / float(bin_size);
        if ( keep_ratio >= 1. ) keep_ratio = 1.;
        bin_keep_ratio[i] = keep_ratio;
        bin_upper_edge[i] = jet_pt_min + ((i+1) * bin_width);
        std::cout << "Bin " << i << " has " << bin_size << " entries, upper edge: " << bin_upper_edge[i] << ", and keep ratio: " << bin_keep_ratio[i] << std::endl;
    }
    
    // Copies events over to the output tree
    std::cout << "Copying entries to flat tree..." << std::endl;
    TRandom3* random_gen = new TRandom3();
    for ( int j = 0 ; j < input_tree->GetEntries() ; j++ ) {
        input_tree->GetEntry(j);
        for ( int i = 0 ; i < n_bins ; i++ ) {
            if ( jet_pt_true > bin_upper_edge[i] ) continue;        // Skips if the pt_true isn't within the bin
            if ( random_gen->Rndm() > bin_keep_ratio[i] ) break;    // Ends the loop if not randomly accepted
            output_tree->Fill();
            break; // Ends the loop once filled
        }
    }
    std::cout << "Flattening complete!" << std::endl;
    output_tree->Write("", TObject::kOverwrite);
    file->Write("", TObject::kOverwrite);
    file->Close();
}

} // End of namespace
