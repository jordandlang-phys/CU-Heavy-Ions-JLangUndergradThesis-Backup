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

const int label_arr_size = 2;
//char line2[100];
//int int1 = sprintf(line2, "2.76 TeV, N_{event} = %i", nEvent);
char* label_arr[label_arr_size] = {
    "FastJet, R = 0.4", "p_{T min, jet} > 5.0 GeV"};

const bool printOut = true;
const bool debug = false;

const int max_jets = 100;
const int max_parts = 400;
//const char* output_file_name = "Jet_ML_Prep.root";

void Jet_ML_Prep(char* file_name, char* output_tree_description, float pt_min, float pt_max, int lowest_jet) {
    
    char output_file_name[200];
    sprintf(output_file_name, "ML_Prep_%s.root", file_name);
    
    char output_tree_name[200];
    sprintf(output_tree_name, "Tree_%s", file_name);
    
    std::cout << "----- Preparing " << output_tree_name << " -----" << std::endl;
    
    // Combined file of jets from thermal and PYTHIA data
    char combined_file_path[200];
    sprintf(combined_file_path, "%s/Comb_%s_Trees.root", dir_data, file_name);
    TFile* combined_file = new TFile(combined_file_path, "READ");
    
    TTree* combJet_tree = (TTree*) combined_file->Get("FastJet_Tree");
    
    std::cout << "Reading combined file." << std::endl;
    
    int     c_jet_n;
    float   c_jet_pt[max_jets];
    float   c_jet_y[max_jets];
    float   c_jet_phi[max_jets];
    float   c_jet_E[max_jets];
    float   c_jet_mass[max_jets];
    float   c_jet_area[max_jets];
    float   c_jet_area_err[max_jets];
    int     c_jet_const_n[max_jets];
    float   c_jet_const_pt[max_jets][max_parts];
    float   c_jet_const_eta[max_jets][max_parts];
    float   c_jet_const_phi[max_jets][max_parts];
    float   c_jet_const_E[max_jets][max_parts];
    
    combJet_tree->SetBranchAddress("jet_n",         &c_jet_n);
    combJet_tree->SetBranchAddress("jet_pt",        c_jet_pt);
    combJet_tree->SetBranchAddress("jet_y",         c_jet_y);
    combJet_tree->SetBranchAddress("jet_phi",       c_jet_phi);
    combJet_tree->SetBranchAddress("jet_E",         c_jet_E);
    combJet_tree->SetBranchAddress("jet_area",      c_jet_area);
    combJet_tree->SetBranchAddress("jet_area_err",  c_jet_area_err);
    combJet_tree->SetBranchAddress("jet_mass",      c_jet_mass);
    combJet_tree->SetBranchAddress("jet_const_n",   c_jet_const_n);
    combJet_tree->SetBranchAddress("jet_const_pt",  c_jet_const_pt);
    combJet_tree->SetBranchAddress("jet_const_eta", c_jet_const_eta);
    combJet_tree->SetBranchAddress("jet_const_phi", c_jet_const_phi);
    combJet_tree->SetBranchAddress("jet_const_E",   c_jet_const_E);
    
    // File of jets from PYTHIA
    
    char pythia_file_path[200];
    sprintf(pythia_file_path, "%s/Pyth_%s_Trees.root", dir_data, file_name);
    TFile* pythia_file = new TFile(pythia_file_path, "READ");
    
    TTree* pythJet_tree = (TTree*) pythia_file->Get("FastJet_Tree");
    
    std::cout << "Reading PYTHIA file." << std::endl;
    
    int     p_jet_n;
    float   p_jet_pt[max_jets];
    float   p_jet_y[max_jets];
    float   p_jet_phi[max_jets];
    float   p_jet_E[max_jets];
    float   p_jet_mass[max_jets];
    float   p_jet_area[max_jets];
    float   p_jet_area_err[max_jets];
    int     p_jet_const_n[max_jets];
    float   p_jet_const_pt[max_jets][max_parts];
    float   p_jet_const_eta[max_jets][max_parts];
    float   p_jet_const_phi[max_jets][max_parts];
    float   p_jet_const_E[max_jets][max_parts];
    
    pythJet_tree->SetBranchAddress("jet_n",         &p_jet_n);
    pythJet_tree->SetBranchAddress("jet_pt",        p_jet_pt);
    pythJet_tree->SetBranchAddress("jet_y",         p_jet_y);
    pythJet_tree->SetBranchAddress("jet_phi",       p_jet_phi);
    pythJet_tree->SetBranchAddress("jet_E",         p_jet_E);
    pythJet_tree->SetBranchAddress("jet_area",      p_jet_area);
    pythJet_tree->SetBranchAddress("jet_area_err",  p_jet_area_err);
    pythJet_tree->SetBranchAddress("jet_mass",      p_jet_mass);
    pythJet_tree->SetBranchAddress("jet_const_n",   p_jet_const_n);
    pythJet_tree->SetBranchAddress("jet_const_pt",  p_jet_const_pt);
    pythJet_tree->SetBranchAddress("jet_const_eta", p_jet_const_eta);
    pythJet_tree->SetBranchAddress("jet_const_phi", p_jet_const_phi);
    pythJet_tree->SetBranchAddress("jet_const_E",   p_jet_const_E);
    
    // Sets output file and plot directory
    char output_file_path[200];
    sprintf(output_file_path, "%s/%s", dir_data, output_file_name);
    TFile* output_file = new TFile(output_file_path, "UPDATE");
    TTree* output_tree = new TTree(output_tree_name, output_tree_description);
    
    char subdir_plots[200];
    sprintf(subdir_plots, "%s/MachineLearning", dir_plots);
    std::__fs::filesystem::create_directories(subdir_plots);
    
    // Plotted Data Output Tree
    // NOTE: This assumes jets have already been sorted from highest E to lowest E by FastJet!
    
    int nEvent = pythJet_tree->GetEntries();
    
    int jet_total_n = 0;
    int jet_hard_n = 0;
    int jet_soft_n = 0;
    float  jet_E_raw[max_jets];
    float  jet_E_reco[max_jets];
    float  jet_E_true[max_jets];
    float  background_E_median[nEvent];
    float  background_pt_median[nEvent];
    float  thermal_E_median[nEvent];
    float  jet_const_n_mean[nEvent];
    float  jet_const_n_median[nEvent];
    
    // ML X value variables (features)
    int    jet_n;
    float  jet_pt_raw[max_jets];
    float  jet_pt_corr[max_jets];
    float  jet_pt_true_pythia[max_jets];
    float  jet_pt_true_paper[max_jets];
    float  jet_mass[max_jets];
    float  jet_area[max_jets];
    float  jet_area_err[max_jets];
    int    jet_const_n[max_jets];
    float  const_pt_mean[max_jets];
    float  const_pt_median[max_jets];
    float  const_1_pt[max_jets];
    float  const_2_pt[max_jets];
    float  const_3_pt[max_jets];
    float  const_4_pt[max_jets];
    float  const_5_pt[max_jets];
    float  const_6_pt[max_jets];
    float  const_7_pt[max_jets];
    float  const_8_pt[max_jets];
    float  const_9_pt[max_jets];
    float  const_10_pt[max_jets];
    float  jet_y[max_jets];
    float  jet_phi[max_jets];
    float  jet_rho[max_jets];
    float  background_rho;
    
    output_tree->Branch("jet_n",                &jet_n);
    output_tree->Branch("jet_pt_raw",           jet_pt_raw,         "jet_pt_raw[100]/F");
    output_tree->Branch("jet_pt_corr",          jet_pt_corr,        "jet_pt_corr[100]/F");
    output_tree->Branch("jet_pt_true_pythia",   jet_pt_true_pythia, "jet_pt_true_pythia[100]/F");
    output_tree->Branch("jet_pt_true_paper",    jet_pt_true_paper,  "jet_pt_true_paper[100]/F");
    output_tree->Branch("jet_mass",             jet_mass,           "jet_mass[100]/F");
    output_tree->Branch("jet_area",             jet_area,           "jet_area[100]/F");
    output_tree->Branch("jet_area_err",         jet_area_err,       "jet_area_err[100]/F");
    output_tree->Branch("jet_const_n",          jet_const_n,        "jet_const_n[100]/I");
    output_tree->Branch("const_pt_mean",        const_pt_mean,      "const_pt_mean[100]/F");
    output_tree->Branch("const_pt_median",      const_pt_median,    "const_pt_median[100]/F");
    output_tree->Branch("const_1_pt",           const_1_pt,         "const_1_pt[100]/F");
    output_tree->Branch("const_2_pt",           const_2_pt,         "const_2_pt[100]/F");
    output_tree->Branch("const_3_pt",           const_3_pt,         "const_3_pt[100]/F");
    output_tree->Branch("const_4_pt",           const_4_pt,         "const_4_pt[100]/F");
    output_tree->Branch("const_5_pt",           const_5_pt,         "const_5_pt[100]/F");
    output_tree->Branch("const_6_pt",           const_6_pt,         "const_6_pt[100]/F");
    output_tree->Branch("const_7_pt",           const_7_pt,         "const_7_pt[100]/F");
    output_tree->Branch("const_8_pt",           const_8_pt,         "const_8_pt[100]/F");
    output_tree->Branch("const_9_pt",           const_9_pt,         "const_9_pt[100]/F");
    output_tree->Branch("const_10_pt",          const_10_pt,        "const_10_pt[100]/F");
    output_tree->Branch("jet_y",                jet_y,              "jet_y[100]/F");
    output_tree->Branch("jet_phi",              jet_phi,            "jet_phi[100]/F");
    output_tree->Branch("jet_rho",              jet_rho,            "jet_rho[100]/F");
    output_tree->Branch("background_rho",       &background_rho);
    
    int jet_const_n_total = 0;
    int jet_true_pythia_counter = 0;
    int jet_true_paper_counter = 0;
    
    for ( int e = 0 ; e < combJet_tree->GetEntries() ; e++ ) {
        combJet_tree->GetEntry(e);
        pythJet_tree->GetEntry(e);
        
        // Finds the median pt of the background jets
        
        float  background_jet_pt[max_jets];
        float  background_jet_area[max_jets];
        float  background_rho_arr[max_jets];
        float  background_area_median;
        int    background_jet_n = 0;
        float  jet_const_n_arr[max_parts];
        int    jet_const_n_event = 0;
        
        for ( int j = 0 ; j < c_jet_n ; j++ ) {
            jet_const_n_total++;
            jet_const_n_arr[j] = c_jet_const_n[j];
            jet_const_n_event++;
            if (j > 1) {
                background_jet_pt[j-2]   = c_jet_pt[j];
                background_jet_area[j-2] = c_jet_area[j];
                background_rho_arr[j-2]  = c_jet_pt[j] / c_jet_area[j];
                background_jet_n++;
            }
        }
        
        background_area_median  = TMath::Median(background_jet_n, background_jet_area);
        background_pt_median[e] = TMath::Median(background_jet_n, background_jet_pt) * background_area_median;
        background_rho          = TMath::Median(background_jet_n, background_rho_arr);
        
        jet_const_n_mean[e] = TMath::Mean(jet_const_n_event, jet_const_n_arr);
        jet_const_n_median[e] = TMath::Median(jet_const_n_event, jet_const_n_arr);
        
        // Loops through jets for machine learning predictors
        jet_n = c_jet_n;
        if (debug) std::cout << "Jets for event " << e << ": " << jet_n << std::endl;
        
        float const_total_pt;
        float const_pythia_pt;
        
        if ( lowest_jet == 0 ) lowest_jet = jet_n;
        
        for ( int j = 0 ; j < lowest_jet ; j++ ) {
            
            // Iterate through pythia jets to match jet pT_true
            int pythia_match = -1;
            
            for ( int k = 0 ; k < p_jet_n ; k++ ) {
                if ( ( pow((p_jet_y[k] - c_jet_y[j]), 2) + pow((p_jet_phi[k] - c_jet_phi[j]), 2) ) < fj_rSquared) {
                    pythia_match = k;
                    if ( debug ) std::cout << "Match found!" << std::endl;
                }
            }
            
            if ( pythia_match < 0 ) continue;
            
            // Resets values for each iteration
            jet_pt_true_paper[j]    = 0;
            jet_pt_true_pythia[j]   = 0;
            jet_pt_raw[j]       = 0;
            jet_pt_corr[j]      = 0;
            jet_mass[j]         = 0;
            jet_area[j]         = 0;
            jet_const_n[j]      = 0;
            const_pt_mean[j]    = 0;
            const_pt_median[j]  = 0;
            const_1_pt[j]       = 0;
            const_2_pt[j]       = 0;
            const_3_pt[j]       = 0;
            const_4_pt[j]       = 0;
            const_5_pt[j]       = 0;
            const_6_pt[j]       = 0;
            const_7_pt[j]       = 0;
            const_8_pt[j]       = 0;
            const_9_pt[j]       = 0;
            const_10_pt[j]      = 0;
            jet_y[j]            = 0;
            jet_phi[j]          = 0;
            jet_rho[j]          = 0;
            
            const_total_pt = 0.;
            const_pythia_pt = 0.;
            
            // Iterate through constituent particles to collect their pT for mean and median
            // This checks EVERY pythia particle in the jet, regardless of the jet match status
            float  jet_const_pt_arr[max_parts];
            
            for ( int p = 0 ; p < c_jet_const_n[j] ; p++ ) {
                jet_const_pt_arr[p] = c_jet_const_pt[j][p];
                const_total_pt += c_jet_const_pt[j][p];
                
                // Iterate through pythia particles to match and find Jet pT True by paper method
                for ( int k = 0 ; k < p_jet_n ; k++) {
                    for ( int q = 0 ; q < p_jet_const_n[k] ; q++ ) {
                        if ( pow(c_jet_const_eta[j][p] -  p_jet_const_eta[k][q], 2) + pow(c_jet_const_phi[j][p] -  p_jet_const_phi[k][q], 2) < pow(0.0001, 2) ) {
                            
                            const_pythia_pt = const_pythia_pt + p_jet_const_pt[k][q];
                            
                            if (debug) std::cout << e << "-" << j << ": Match found! --- " << p_jet_const_pt[k][q] << " - " << const_pythia_pt << std::endl;
                        }
                    }
                }
            }
            
            std::sort(jet_const_pt_arr, jet_const_pt_arr + jet_const_n[j], greater<>());
            
            jet_pt_raw[j]       = c_jet_pt[j];
            jet_pt_corr[j]      = c_jet_pt[j] - (background_rho * c_jet_area[j]);
            jet_mass[j]         = c_jet_mass[j];
            jet_area[j]         = c_jet_area[j];
            jet_const_n[j]      = c_jet_const_n[j];
            const_pt_mean[j]    = TMath::Mean(jet_const_n[j], jet_const_pt_arr);
            const_pt_median[j]  = TMath::Median(jet_const_n[j], jet_const_pt_arr);
            const_1_pt[j]       = jet_const_pt_arr[0];
            const_2_pt[j]       = jet_const_pt_arr[1];
            const_3_pt[j]       = jet_const_pt_arr[2];
            const_4_pt[j]       = jet_const_pt_arr[3];
            const_5_pt[j]       = jet_const_pt_arr[4];
            const_6_pt[j]       = jet_const_pt_arr[5];
            const_7_pt[j]       = jet_const_pt_arr[6];
            const_8_pt[j]       = jet_const_pt_arr[7];
            const_9_pt[j]       = jet_const_pt_arr[8];
            const_10_pt[j]      = jet_const_pt_arr[9];
            jet_y[j]            = c_jet_y[j];
            jet_phi[j]          = c_jet_phi[j];
            jet_rho[j]          = background_rho;
            
            jet_pt_true_paper[j]  = jet_pt_raw[j] * (const_pythia_pt / const_total_pt);
            
//            if ( jet_pt_true_paper[j] != 0 && (e % 1000) == 0 ) std::cout << e << "-" << j << ": Truth_Paper Jet: " << jet_pt_true_paper[j] << std::endl;
            
            if ( pythia_match < 0 ) continue;
            
            jet_pt_true_pythia[j] = p_jet_pt[pythia_match];
                        
            jet_true_pythia_counter++;
            if ( jet_pt_true_pythia[j] != 0 && (e % 1000) == 0 ) std::cout << e << "-" << j << ": Truth_PYTHIA Jet: " << jet_pt_true_pythia[j] << " ----- " << std::endl;
        }
        output_tree->Fill();
    }
    
    output_tree->Write("", TObject::kOverwrite);
    output_file->Write("", TObject::kOverwrite);
    
    delete output_tree;
    
    std::cout << "Data written to file." << std::endl;
    
    pythia_file->Close();
    combined_file->Close();
    output_file->Close();
    
    delete pythia_file;
    delete combined_file;
    delete output_file;
    
    std::cout << "Files closed." << std::endl;
    
    std::cout << "----- Completed " << output_tree_name << " -----" << std::endl;
}

void Jet_ML_Prep_ROOT() {
    
//    char output_file_path[200];
//    sprintf(output_file_path, "%s/%s", dir_data, output_file_name);
//    TFile* output_file = new TFile(output_file_path, "UPDATE");
//    output_file->Close();
    
//    Jet_ML_Prep(
//        "40_60_Train",
//        "TRAINING data for machine learning. 40-60 GeV with pT^3 bias.",
//        40., 60.);
//
//    Jet_ML_Prep(
//        "40_60_Test",
//        "TESTING data for machine learning. 40-60 GeV with pT^3 bias.",
//        40., 60.);

    Jet_ML_Prep(
        "10_90_B0_Train",
        "TRAINING data for machine learning. 10-90 GeV with no bias.",
        10., 90.,
        0); // Top n jets to test (NOT index)
    
    Jet_ML_Prep(
        "40_60_B0_Test",
        "TESTING data for machine learning. 40-60 GeV with no bias.",
        40., 60.,
        0); // Top n jets to test (NOT index)
    
    Jet_ML_Prep(
        "10_90_B4_Train",
        "TRAINING data for machine learning. 10-90 GeV with a bias of pT^4.",
        10., 90.,
        0); // Top n jets to test (NOT index)
    
    Jet_ML_Prep(
        "40_60_B4_Test",
        "TESTING data for machine learning. 40-60 GeV with a bias of pT^4.",
        40., 60.,
        0); // Top n jets to test (NOT index)
    
    Jet_ML_Prep(
        "10_90_B8_Train",
        "TRAINING data for machine learning. 10-90 GeV with a bias of pT^8.",
        10., 90.,
        0); // Top n jets to test (NOT index)
    
    Jet_ML_Prep(
        "40_60_B8_Test",
        "TESTING data for machine learning. 40-60 GeV with a bias of pT^8.",
        40., 60.,
        0); // Top n jets to test (NOT index)
    
}
