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

#include "jet_ml_constants.h"
using namespace Jet_ML_Constants;

const int label_arr_size = 2;
//char line2[100];
//int int1 = sprintf(line2, "2.76 TeV, N_{event} = %i", nEvent);
char* label_arr[label_arr_size] = {
    "FastJet, R = 0.4", "p_{T min, jet} > 5.0 GeV"};

const bool printOut = true;
const bool debug    = false;

const int max_jets = 100;
const int max_parts = 400;

//const char* output_file_name = "Jet_ML_Prep.root";

bool Check_File_Compatibility(
    char  file_path_A[500],
    char  file_path_B[500],
    ) {
    
    bool files_ok = true;
    
    TFile* file_A = new TFile(file_path_A, "READ");
    TTree* tree_A;
    if (file1->GetListOfKeys()->Contains("Generator_Parameters")) {
        tree_A = (TTree*) file_A->Get("Generator_Parameters"); }
    else if (file1->GetListOfKeys()->Contains("Pythia_Generator_Parameters")) {
        tree_A = (TTree*) file_A->Get("Pythia_Generator_Parameters"); }
    else if (file1->GetListOfKeys()->Contains("Thermal_Generator_Parameters")) {
        tree_A = (TTree*) file_A->Get("Thermal_Generator_Parameters"); }
    else files_ok = false;
        
    TFile* file_B = new TFile(file_path_B, "READ");
    TTree* tree_B;
    if (file1->GetListOfKeys()->Contains("Generator_Parameters")) {
        tree_B = (TTree*) file_B->Get("Generator_Parameters"); }
    else if (file1->GetListOfKeys()->Contains("Pythia_Generator_Parameters")) {
        tree_B = (TTree*) file_B->Get("Pythia_Generator_Parameters"); }
    else if (file1->GetListOfKeys()->Contains("Thermal_Generator_Parameters")) {
        tree_B = (TTree*) file_B->Get("Thermal_Generator_Parameters"); }
    else files_ok = false;
    
    if (files_ok) {
        int   a_event_count;
        float a_pt_bias_power;
        float a_jet_pt_min;
        float a_jet_pt_max;
        
        tree_A->SetBranchAddress("event_count",       &a_event_count);
        tree_A->SetBranchAddress("pt_bias_power",     &a_pt_bias_power);
        tree_A->SetBranchAddress("jet_pt_min",        &a_jet_pt_min);
        tree_A->SetBranchAddress("jet_pt_max",        &a_jet_pt_max);
        
        int   b_event_count;
        float b_pt_bias_power;
        float b_jet_pt_min;
        float b_jet_pt_max;
        
        tree_B->SetBranchAddress("event_count",       &b_event_count);
        tree_B->SetBranchAddress("pt_bias_power",     &b_pt_bias_power);
        tree_B->SetBranchAddress("jet_pt_min",        &b_jet_pt_min);
        tree_B->SetBranchAddress("jet_pt_max",        &b_jet_pt_max);
        
        tree_A->GetEntry(0);
        tree_B->GetEntry(0);
        
        // Checks that files use the same generator data
        if ((a_event_count   != b_event_count)   ||
            (a_pt_bias_power != b_pt_bias_power) ||
            (a_jet_pt_min    != b_jet_pt_min)    ||
            (a_jet_pt_max    != b_jet_pt_max) ) {
            files_ok = false;
            std::cout << "Input files are incompatible! Generator data is not identical between files." << std::endl;
            }
        }
        
    return files_ok;
    }

void Jet_ML_Prep(
    char  output_base_name[100],        // Base/stem of output file name (no file extensions or prefixes).
    char  output_directory[400],        // Path to desired output directory.
    char  comb_file_name[100],          // Full name of combined output file
    char  pyth_file_name[100],          // Full name of combined output file
    float jet_pt_min = NULL,
    float jet_pt_max = NULL,
    int   lowest_jet = 0,               // Only accepts the top [lowest_jet] jets. If lowest_jet = 0, accepts all jets
    float fastjet_match_radius = 0.3    // Two jets are matched if they are within this radius squared.
    ) {
    
    // Combined file of jets from thermal and PYTHIA data
    char combined_file_path[200];
    sprintf(combined_file_path, "%s/Comb_%s.root", dir_data, file_name);
    TFile* combined_file = new TFile(combined_file_path, "READ");
    
    TTree* combJet_tree = (TTree*) combined_file->Get("FastJet_Tree");
    TTree* pyth_param_tree = (TTree*) combined_file->Get("Pyth_Generator_Parameters");
    
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
    
    if ( !Check_File_Compatibility(combined_file_path, pythia_file_path) ) break;
    
    // Sets output file and plot directory
    char output_file_path[200];
    sprintf(output_file_path, "%s/%s", dir_data, output_file_name);
    TFile* output_file = new TFile(output_file_path, "UPDATE");
    
    TTree* pyth_param_tree_copy = pyth_param_tree->CopyTree("Pyth_Generator_Parameters");
    pyth_param_tree_copy->Write("", TObject::kOverwrite);
    
    int   p_gen_event_count;
    float p_gen_pt_bias_power;
    float p_gen_jet_pt_min;
    float p_gen_jet_pt_max;
    
    pyth_param_tree->SetBranchAddress("event_count",       &p_gen_event_count);
    pyth_param_tree->SetBranchAddress("pt_bias_power",     &p_gen_pt_bias_power);
    pyth_param_tree->SetBranchAddress("jet_pt_min",        &p_gen_jet_pt_min);
    pyth_param_tree->SetBranchAddress("jet_pt_max",        &p_gen_jet_pt_max);
    
    pyth_param_tree->GetEntry(0);
    
    char output_file_path[200];
    sprintf(output_file_path, "%s/ML_Prep_%s_B%.0f_%.0f_%.0f_N%i.root", output_directory, output_base_name, gen_pt_bias_power, gen_jet_pt_min, gen_jet_pt_max, gen_event_count);
    
    char output_tree_name[200];
    sprintf(output_tree_name, "Jet_ML_Tree", output_base_name);
    char output_tree_description[200];
    sprintf(output_tree_description, "Generated from bias of ", output_base_name);
    TTree* output_tree = new TTree(output_tree_name, output_tree_description);
    
    std::cout << "----- Preparing " << output_tree_name << " -----" << std::endl;
    
    char plot_directory[200];
    sprintf(plot_directory, "%s/MachineLearning", output_directory);
    std::__fs::filesystem::create_directories(subdir_plots);
    
    // Plotted Data Output Tree
    // NOTE: This assumes jets have already been sorted from highest E to lowest E by FastJet!
    
    int    nEvent = pythJet_tree->GetEntries();
    
    float  background_pt_median[nEvent];
    float  thermal_E_median[nEvent];
    float  jet_const_n_mean[nEvent];
    float  jet_const_n_median[nEvent];
    
    // ML X value variables (features)
    int    jet_index;
    float  jet_pt_raw;
    float  jet_pt_corr;
    float  jet_pt_true_pythia;
    float  jet_pt_true_paper;
    float  jet_mass;
    float  jet_area;
    float  jet_area_err;
    int    jet_const_n;
    float  const_pt_mean;
    float  const_pt_median;
    float  const_1_pt;
    float  const_2_pt;
    float  const_3_pt;
    float  const_4_pt;
    float  const_5_pt;
    float  const_6_pt;
    float  const_7_pt;
    float  const_8_pt;
    float  const_9_pt;
    float  const_10_pt;
    float  jet_y;
    float  jet_phi;
    float  jet_rho;
    float  background_rho;
    
    output_tree->Branch("jet_index",            &jet_index);
    output_tree->Branch("jet_pt_raw",           &jet_pt_raw,         "jet_pt_raw/F");
    output_tree->Branch("jet_pt_corr",          &jet_pt_corr,        "jet_pt_corr/F");
    output_tree->Branch("jet_pt_true_pythia",   &jet_pt_true_pythia, "jet_pt_true_pythia/F");
    output_tree->Branch("jet_pt_true_paper",    &jet_pt_true_paper,  "jet_pt_true_paper/F");
    output_tree->Branch("jet_mass",             &jet_mass,           "jet_mass/F");
    output_tree->Branch("jet_area",             &jet_area,           "jet_area/F");
    output_tree->Branch("jet_area_err",         &jet_area_err,       "jet_area_err/F");
    output_tree->Branch("jet_const_n",          &jet_const_n,        "jet_const_n/I");
    output_tree->Branch("const_pt_mean",        &const_pt_mean,      "const_pt_mean/F");
    output_tree->Branch("const_pt_median",      &const_pt_median,    "const_pt_median/F");
    output_tree->Branch("const_1_pt",           &const_1_pt,         "const_1_pt/F");
    output_tree->Branch("const_2_pt",           &const_2_pt,         "const_2_pt/F");
    output_tree->Branch("const_3_pt",           &const_3_pt,         "const_3_pt/F");
    output_tree->Branch("const_4_pt",           &const_4_pt,         "const_4_pt/F");
    output_tree->Branch("const_5_pt",           &const_5_pt,         "const_5_pt/F");
    output_tree->Branch("const_6_pt",           &const_6_pt,         "const_6_pt/F");
    output_tree->Branch("const_7_pt",           &const_7_pt,         "const_7_pt/F");
    output_tree->Branch("const_8_pt",           &const_8_pt,         "const_8_pt/F");
    output_tree->Branch("const_9_pt",           &const_9_pt,         "const_9_pt/F");
    output_tree->Branch("const_10_pt",          &const_10_pt,        "const_10_pt/F");
    output_tree->Branch("jet_y",                &jet_y,              "jet_y/F");
    output_tree->Branch("jet_phi",              &jet_phi,            "jet_phi/F");
    output_tree->Branch("jet_rho",              &jet_rho,            "jet_rho/F");
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
        if (debug) std::cout << "Jets for event " << e << ": " << c_jet_n << std::endl;
        
        float const_total_pt;
        float const_pythia_pt;
        int pythia_match = -1;
        
        if ( lowest_jet == 0 ) lowest_jet = c_jet_n;
        
        for ( int j = 0 ; j < lowest_jet ; j++ ) {
            
            // Iterate through pythia jets to match jet pT_true
            pythia_match = -1;
            
            for ( int k = 0 ; k < p_jet_n ; k++ ) {
                if ( pythia_match >= 0 ) continue;
                if ( ( pow((p_jet_y[k] - c_jet_y[j]), 2) + pow((p_jet_phi[k] - c_jet_phi[j]), 2) ) < fj_rSquared) {
                    if ( p_jet_pt[k] > pt_min && p_jet_pt[k] < pt_max ) pythia_match = k;
                    if ( debug ) std::cout << "Match found!" << std::endl;
                }
            }
            
            if ( pythia_match < 0 ) continue;
            if ( p_jet_pt[pythia_match] < pt_min && p_jet_pt[pythia_match] > pt_max ) continue;
            
            // Iterate through constituent particles to collect their pT for mean and median
            float  jet_const_pt_arr[max_parts];
            
            for ( int p = 0 ; p < c_jet_const_n[j] ; p++ ) {
                jet_const_pt_arr[p] = c_jet_const_pt[j][p];
                const_total_pt += c_jet_const_pt[j][p];
            }
            
            std::sort(jet_const_pt_arr, jet_const_pt_arr + jet_const_n, greater<>());
            
            jet_index       = j;
            jet_pt_raw      = c_jet_pt[j];
            jet_pt_corr     = c_jet_pt[j] - (background_rho * c_jet_area[j]);
            jet_mass        = c_jet_mass[j];
            jet_area        = c_jet_area[j];
            jet_const_n     = c_jet_const_n[j];
            const_pt_mean   = TMath::Mean(jet_const_n, jet_const_pt_arr);
            const_pt_median = TMath::Median(jet_const_n, jet_const_pt_arr);
            const_1_pt      = jet_const_pt_arr[0];
            const_2_pt      = jet_const_pt_arr[1];
            const_3_pt      = jet_const_pt_arr[2];
            const_4_pt      = jet_const_pt_arr[3];
            const_5_pt      = jet_const_pt_arr[4];
            const_6_pt      = jet_const_pt_arr[5];
            const_7_pt      = jet_const_pt_arr[6];
            const_8_pt      = jet_const_pt_arr[7];
            const_9_pt      = jet_const_pt_arr[8];
            const_10_pt     = jet_const_pt_arr[9];
            jet_y           = c_jet_y[j];
            jet_phi         = c_jet_phi[j];
            jet_rho         = background_rho;
            
            jet_pt_true_pythia = p_jet_pt[pythia_match];
            jet_pt_true_paper = 0.0;
                        
            jet_true_pythia_counter++;
            if ( jet_pt_true_pythia != 0 && (e % 1000) == 0 ) std::cout << e << "-" << j << ": Truth_PYTHIA Jet: " << jet_pt_true_pythia << " ----- " << std::endl;
            
            output_tree->Fill();
        }
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
    
//    Jet_ML_Prep(
//        "10_90_B8_Test",
//        "TEST data prepped for ML. 10-90 GeV, bias of pT^8, includes all jets.",
//        10., 90.,
//        0); // Top n jets to test (NOT index, use 0 for all jets)

    Jet_ML_Prep(
        "10_90_B8_Train",
        "TRAIN data prepped for ML. 10-90 GeV, bias of pT^8, includes all jets.",
        10., 90.,
        0); // Top n jets to test (NOT index, use 0 for all jets)
//
//    Jet_ML_Prep(
//        "40_60_B8_Test",
//        "TESTING data for machine learning. 40-60 GeV with a bias of pT^8.",
//        40., 60.,
//        0); // Top n jets to test (NOT index, use 0 for all jets)
    
}
