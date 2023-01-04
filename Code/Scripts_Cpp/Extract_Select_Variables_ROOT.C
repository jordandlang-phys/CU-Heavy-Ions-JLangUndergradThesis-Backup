#include "TFile.h"
#include "TTree.h"
#include "_header.h"
using namespace Project_Constants;

void Extract_Select_Variables(char* file_name) {
    
    char input_file_path[200];
    sprintf(input_file_path, "%s/%s.csv", dir_data, file_name);
    
    fstream csv_file;
    csv_file.open( input_file_path, ios::in );
    
    std::cout << "Accessed input file." << std::endl;
    
    char output_file_path[200];
    sprintf(output_file_path, "%s/%s.root", dir_data, file_name);
    TFile* output_file = new TFile(output_file_path, "RECREATE");
    TTree* output_tree = new TTree("Tree_ML", "TTree of selected variables from machine learning: jet area, jet pT_raw, jet pT_true, jet pT_corrected, and jet pT using linear regression (lr), random forest regression (rf), and neural network regression (nn)");
    
    std::cout << "Created output file." << std::endl;
    
    output_tree->ReadFile(
        input_file_path,
        "jet_area/F:jet_pt_raw/F:jet_pt_true/F:jet_pt_corr/F:jet_pt_ml_lr/F:jet_pt_ml_rf/F:jet_pt_ml_nn/F", ',');
    
    output_tree->Write("", TObject::kOverwrite);
    
    delete output_tree;
    
    output_file->Close();
    
    delete output_file;
}

void Extract_Select_Variables_ROOT() {
    Extract_Select_Variables("ML_Results_10_90_Tree_40_60_Test_12feat_ptTruePythia");
}
