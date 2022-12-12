#include "TFile.h"
#include "TTree.h"
#include "_header.h"
#include <fstream>
#include <vector>
#include <sstream>
using namespace Project_Constants;

void Build_Results_TTree_FromCSV(
    char* file_dir,         // file directory path
    char* input_file_name,  // input file name
    char* output_file_name, // output file name
    char* output_tree_name  // name for TTree in output file
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
    char* input_file_name,      // Full file name with path for input .root file
    char* results_tree_name,    // Name of TTree to get ML data from
    char* weights_tree_name,    // Name of TTree to get weights for ML estimator
    char* plot_name,
    float test_pt_min,
    float test_pt_max
    ) {
    
    TFile* input_file = new TFile(input_file_name, "READ");
    TTree* tree_results = (TTree*) input_file->Get(input_tree_name);
    TTree* tree_weights = (TTree*) input_file->Get(weights_tree_name);
    std::cout << "Input file accessed." << std::endl;
    
    float jet_pt_true;
    float jet_pt_reco;
    float jet_pt_corr;
    
    tree_results->SetBranchAddress("jet_pt_true", &jet_pt_true);
    tree_results->SetBranchAddress("jet_pt_reco", &jet_pt_reco);
    tree_results->SetBranchAddress("jet_pt_reco", &jet_pt_reco);
    
    float x_min = test_pt_min - 20.;
    float x_max = test_pt_max + 20.;
    int  x_bins = int((x_max - x_min) / 1);
    
    TH1F* th1f_jet_pt_true = new TH1F("th1f_jet_pt_reco", "True Jet pT; p_{T, ch jet} [GeV]; Probability Density", x_bins, x_min, x_max);
    TH1F* th1f_jet_pt_reco = new TH1F("th1f_jet_pt_reco", "Jet pT Reconstructed by ML; p_{T, ch jet} [GeV]; Probability Density", x_bins, x_min, x_max);
    TH1F* th1f_jet_pt_corr = new TH1F("th1f_jet_pt_reco", "Jet pT from Simple Correction; p_{T, ch jet} [GeV]; Probability Density", x_bins, x_min, x_max);
    
    th1f_jet_pt_true->Sumw2();
    th1f_jet_pt_reco->Sumw2();
    th1f_jet_pt_corr->Sumw2();
    
    for ( int j=0 ; j < tree_results->GetEntries() ; j++ ) {
        tree_results->GetEntry();
        
        th1f_jet_pt_true->Fill(jet_pt_true);
        th1f_jet_pt_reco->Fill(jet_pt_reco);
        th1f_jet_pt_corr->Fill(jet_pt_corr);
    }
    
    // STOPPED HERE FOR THE NIGHT
    
    
    
}



void Jet_ML_Plotter_ROOT() {
    
    char* file_dir = dir_data;
    char* output_file_name = "ML_Results/Train_B8_10_90.root";
    
    char* ml_weights_array[3][2] = {
        {"ML_Results/Train_B8_F1_10_90_LR_Coeffs.csv", "Tree_Train_B8_F1_10_90_LR_Coeffs"},
        {"ML_Results/Train_B8_F3_10_90_LR_Coeffs.csv", "Tree_Train_B8_F3_10_90_LR_Coeffs"},
        {"ML_Results/Train_B8_F12_10_90_LR_Coeffs.csv", "Tree_Train_B8_F12_10_90_LR_Coeffs"}
    };
    
    for ( int i=0 ; i < 3 ; i++ ) {
        Build_Weights_TTree_FromCSV(
            file_dir,
            ml_weights_array[i][0],
            output_file_name,
            ml_weights_array[i][1]
            );
    }
    
    char* ml_results_array[21][2] = {
        {"ML_Results/Train_B8_F1_10_90_Test_18_22.csv", "Tree_Train_B8_F1_10_90_Test_18_22"},
        {"ML_Results/Train_B8_F1_10_90_Test_28_32.csv", "Tree_Train_B8_F1_10_90_Test_28_32"},
        {"ML_Results/Train_B8_F1_10_90_Test_38_42.csv", "Tree_Train_B8_F1_10_90_Test_38_42"},
        {"ML_Results/Train_B8_F1_10_90_Test_48_52.csv", "Tree_Train_B8_F1_10_90_Test_48_52"},
        {"ML_Results/Train_B8_F1_10_90_Test_58_62.csv", "Tree_Train_B8_F1_10_90_Test_58_62"},
        {"ML_Results/Train_B8_F1_10_90_Test_68_72.csv", "Tree_Train_B8_F1_10_90_Test_68_72"},
        {"ML_Results/Train_B8_F1_10_90_Test_78_82.csv", "Tree_Train_B8_F1_10_90_Test_78_82"},
        {"ML_Results/Train_B8_F3_10_90_Test_18_22.csv", "Tree_Train_B8_F3_10_90_Test_18_22"},
        {"ML_Results/Train_B8_F3_10_90_Test_28_32.csv", "Tree_Train_B8_F3_10_90_Test_28_32"},
        {"ML_Results/Train_B8_F3_10_90_Test_38_42.csv", "Tree_Train_B8_F3_10_90_Test_38_42"},
        {"ML_Results/Train_B8_F3_10_90_Test_48_52.csv", "Tree_Train_B8_F3_10_90_Test_48_52"},
        {"ML_Results/Train_B8_F3_10_90_Test_58_62.csv", "Tree_Train_B8_F3_10_90_Test_58_62"},
        {"ML_Results/Train_B8_F3_10_90_Test_68_72.csv", "Tree_Train_B8_F3_10_90_Test_68_72"},
        {"ML_Results/Train_B8_F3_10_90_Test_78_82.csv", "Tree_Train_B8_F3_10_90_Test_78_82"},
        {"ML_Results/Train_B8_F12_10_90_Test_18_22.csv", "Tree_Train_B8_F12_10_90_Test_18_22"},
        {"ML_Results/Train_B8_F12_10_90_Test_28_32.csv", "Tree_Train_B8_F12_10_90_Test_28_32"},
        {"ML_Results/Train_B8_F12_10_90_Test_38_42.csv", "Tree_Train_B8_F12_10_90_Test_38_42"},
        {"ML_Results/Train_B8_F12_10_90_Test_48_52.csv", "Tree_Train_B8_F12_10_90_Test_48_52"},
        {"ML_Results/Train_B8_F12_10_90_Test_58_62.csv", "Tree_Train_B8_F12_10_90_Test_58_62"},
        {"ML_Results/Train_B8_F12_10_90_Test_68_72.csv", "Tree_Train_B8_F12_10_90_Test_68_72"},
        {"ML_Results/Train_B8_F12_10_90_Test_78_82.csv", "Tree_Train_B8_F12_10_90_Test_78_82"},
    };
    
    for ( int i=0 ; i < 21 ; i++ ) {
        Build_Results_TTree_FromCSV(
            file_dir,
            ml_results_array[i][0],
            output_file_name,
            ml_results_array[i][1]
            );
    }
    
}

