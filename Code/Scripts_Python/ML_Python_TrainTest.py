import sklearn
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import accuracy_score
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPRegressor
import numpy as np
from statistics import mean, median
import math
from array import array
from os import path
from array import array
import csv



def Build_FeatureArrays_FromCSV(
    csv_file_path
    ) :
    """
    Simple function to build an array from a csv that has the following columns (number indicates index):
    
    0  jet_pt_raw,      1  jet_pt_corr,     2  jet_mass,        3  jet_area, 
    4  jet_area_err,    5  jet_const_n,     6  const_pt_mean,   7  const_pt_median, 
    8  const_1_pt,      9  const_2_pt,      10 const_3_pt,      11 const_4_pt,
    12 const_5_pt,      13 const_6_pt,      14 const_7_pt,      15 const_8_pt,
    16 const_9_pt,      17 const_10_pt,     18 jet_y,           19 jet_phi,
    20 jet_rho,         21 jet_pt_true,     22 jet_index
    """
    
    input_csv  = open(csv_file_path)
    csv_reader = csv.reader(input_csv)
    
    X_values  = []
    y_values  = []
    sc_values = []
    
    skip_header_toggle = True  # Skips the header row of the csv
    jet_counter = 0            # Counts total number of jets in file
    
    print("Preparing to collect data from csv backup file...")
    
    for row in csv_reader:
        if skip_header_toggle:
            skip_header_toggle = False
            continue
        X_values.append([
            float(row[0]),  float(row[1]),  float(row[2]),  float(row[3]), 
            float(row[4]),  float(row[5]),  float(row[6]),  float(row[7]),
            float(row[8]),  float(row[9]),  float(row[10]), float(row[11]),
            float(row[12]), float(row[13]), float(row[14]), float(row[15]),
            float(row[16]), float(row[17]), float(row[18]), float(row[19]),
            float(row[20])])
        y_values.append(float(row[21]))
        sc_values.append(float(row[1]))
        
        jet_counter   = jet_counter + 1
        
        if jet_counter % 10000 == 0 : print(f"Jet: {jet_counter:2.0f} | pTraw: {float(row[0]):3.3f} | pTcorr: {float(row[1]): 3.3f} | pTtrue: {float(row[21]): 5.3f}")
    
    input_csv.close()
    print("Backup .csv file closed.")

    print(f"All data transferred to array. Testing with {jet_counter} jets.\n")
    print(f"Data set lengths: {len(X_values)} / {len(y_values)} / {len(sc_values)}")
    
    return X_values, y_values, sc_values



def Write_MLResults_ToCSV(
    output_filename,  # File path with name (MUST include '.csv' at end)
    pt_true_array,    # Array of pt_true values (probably 'y_test' array)
    pt_sc_array,      # Array of pt values from rho*area method (probably 'sc_test' array)
    pt_reco_array,    # Array of pt_reco values from ML (probably 'results' array)
    feature_array,    # Array of pt_true values (probably 'X_test_select' array)
    feature_label     # Array of feature labels
    ) :
    """
    Writes ML results to a csv file
    """
    # Create csv file
    output_csv  = open(output_filename, 'w', newline='')
    csv_writer  = csv.writer(output_csv)
    
    # Build and write header
    csv_header = ['jet_pt_true', 'jet_pt_reco']
    for label in feature_label:
        csv_header.append(label)
    csv_writer.writerow(csv_header)
    
    # Add data
    for i in range(len(pt_true_array)):
        csv_row = [pt_true_array[i], pt_reco_array[i]]
        for feature in feature_array[i]:
            csv_row.append(feature)
        csv_writer.writerow(csv_row)
        
    output_csv.close()
    print("ML results .csv file closed.")
    
    return


def Write_MLWeights_ToCSV(
    output_filename,  # File path with name (MUST include '.csv' at end)
    weight_array,     # Array of feature weight values (probably 'lr_coeffs' or 'rf_importance')
    feature_label     # Array of feature labels
    ) :
    """
    Writes ML results to a csv file
    """
    # Create csv file
    output_csv  = open(output_filename, 'w', newline='')
    csv_writer  = csv.writer(output_csv)
    
    # Build and write header
    csv_header = []
    for label in feature_label:
        csv_header.append(label)
    csv_writer.writerow(csv_header)
    
    # Add data
    csv_row = []
    for weight in weight_array:
        csv_row.append(weight)
    csv_writer.writerow(csv_row)
        
    output_csv.close()
    print("ML weights .csv file closed.")
    
    return



def Train_LinearRegression(
    X_train,                     # Array of array of input features
    y_train,                     # Array of target values
    features_labels,             # Array of string labels for each feature
    use_scaler = True    # If True, uses StandardScalar
    ):
    print("\n----- Fitting Linear Regression Estimator -----\n")
    
    # Sets Scaler
    if use_scaler:
        scaler = StandardScaler(with_mean=True, with_std=True)
        print("\nUsing StandardScaler. Data will be recentered and normalized.\n")
    else:
        scaler = StandardScaler(with_mean=False, with_std=False)
        print("\nData will not be rescaled.\n")
        
    # Creates a linear regression model in a pipeline
    lr_estimator = LinearRegression()
    lr_pipeline = make_pipeline(
        scaler,
        lr_estimator )
    lr_coeffs = 0
    
    # Fits the regression model
    output = lr_pipeline.fit(X_train, y_train)
    print("\nLinear Regression Fit:\n", output)

    # Outputs regression coefficients
    lr_coeffs = lr_estimator.coef_
    intercept = lr_estimator.intercept_
    print(intercept)
    print(lr_coeffs)
    lr_coeffs = np.append(lr_coeffs, intercept)
    print(lr_coeffs)

    print("Regression Coefficients:")
    for i in range(len(lr_coeffs)) :
        print(features_labels[i], lr_coeffs[i])
    
        
    print(type(lr_pipeline))
    
    return lr_pipeline, lr_coeffs



def Train_RandomForestRegression(
    X_train,                     # Array of array of input features
    y_train,                     # Array of target values
    features_labels,             # Array of string labels for each feature
    use_scaler = True    # If True, uses StandardScalar
    ):
    print("\n----- Fitting Random Forest Regression Estimator -----\n")    
    
    # Sets Scaler
    if use_scaler:
        scaler = StandardScaler(with_mean=True, with_std=True)
        print("\nUsing StandardScaler. Data will be recentered and normalized.\n")
    else:
        scaler = StandardScaler(with_mean=False, with_std=False)
        print("\nData will not be rescaled.\n")
        
    # Creates a random forest estimator in a pipeline
    rf_estimator = RandomForestRegressor()
    rf_pipeline = make_pipeline(
        scaler,
        rf_estimator )
    rf_features = 0

    # Fits the regression model
    output = rf_pipeline.fit(X_train, y_train)
    print("\nRandom Tree Regression Fit:\n", output)

    # Outputs feature importances
    rf_features_arr = rf_estimator.feature_importances_

    print("Feature Importance:")
    for i in range(len(rf_features)) :
        print(features_labels[i], rf_features[i])
    
    return rf_pipeline, rf_features



def Train_MLPRegression(
    X_train,                     # Array of array of input features
    y_train,                     # Array of target values
    features_labels,             # Array of string labels for each feature
    use_scaler = True,   # If True, uses StandardScalar
    max_iterations = 100         # Maximum number of iterations to use for training
    ):
    print("\n----- Fitting Neural Network Regression Estimator -----\n")   
    
    # Sets Scaler
    if use_scaler:
        scaler = StandardScaler(with_mean=True, with_std=True)
        print("\nUsing StandardScaler. Data will be recentered and normalized.\n")
    else:
        scaler = StandardScaler(with_mean=False, with_std=False)
        print("\nData will not be rescaled.\n")
        
    # Creates a MLP estimator in a pipeline
    mlp_estimator = MLPRegressor(max_iter=max_iterations)
    mlp_pipeline = make_pipeline(
        scaler,
        mlp_estimator)

    # Fits the regression model
    output = mlp_pipeline.fit(X_train, y_train)
    print("\nMultilayer Perceptron Regression Fit:\n", output)
    
    return mlp_pipeline



def Train_All_Estimators(
    X_train,             # Array of array of input features
    y_train,             # Array of target values
    features_labels,     # Array of string labels for each feature
    use_scaler = True,   # If True, uses StandardScalar
    use_lr = True,       # If True, trains a linear regression estimator
    use_rf = True,       # If True, trains a random forest regression estimator
    use_mlp = True       # If True, trains a multilayer perceptron regression estimator
    ):
    """
    Function for training Machine Learning Estimators.
    Takes in feature array (X_train) - this is actually an array of arrays, 
    a target array (y_train), and an array of feature labels (feature_arr_labels).
    
    Note that X_train and y_train MUST be the same length, 
    and feature_arr_lables should be 1 longer than the number of features used.
    
    Returns:
    lr_pipeline      Trained sklearn LinearRegression estimator
    rf_pipeline      Trained sklearn RandomForestRegressor estimator
    mlp_pipeline     Trained sklearn MLPRegressor estimator
    lr_coeffs_arr    Array of coefficients used by lr_pipeline, corresponding to each input feature
    rf_features_arr  Array of feature importances used by rf_pipeline, corresponding to each input feature
    """
    
    print("\n===== Running Train_Estimators Function =====\n")
    
    # --- LINEAR REGRESSION ---
    
    if use_lr:
        Train_LinearRegression(X_train, y_train, features_labels, use_scaler = use_scaler)
    else:
        lr_pipeline = False
        lr_coeffs_arr = False

    
    # --- RANDOM FOREST REGRESSION ---
    
    if use_rf:
        Train_RandomForestRegression(X_train, y_train, features_labels, use_scaler = use_scaler)
    else:
        rf_pipeline = False
        rf_features_arr = False
        
        
    # --- MULTILAYER PERCEPTRON REGRESSION ---
    
    if use_mlp:
        Train_MLPRegression(X_train, y_train, features_labels, use_scaler = use_scaler)
    else:
        mlp_pipeline = False
    
    return lr_pipeline, rf_pipeline, mlp_pipeline, lr_coeffs_arr, rf_features_arr

    
    
def Test_Estimator(
    ml_pipeline, # Trained sklearn machine learning estimator pipeline
    X_test,      # Array of arrays of input features
    y_test       # Array of pT truth values for comparison
    ) :
    """
    Simple function to test a trained sklearn estimator on new data.
    
    Returns:
    results        Array of predicted values from input features
    results_delta  Array of differences between prediction and truth (prediction - truth)
    """
    
    if not ml_pipeline:
        print("No pipeline provided!")
        return False, False
    
    results = ml_pipeline.predict(X_test)
    results_delta = []
    
    for j in range(len(results)):
        results_delta.append(results[j] - y_test[j])
    
    return results, results_delta

    

def Test_All_Estimators(
    X_test,          # Array of arrays of input features
    y_test,          # Array of pT truth values for comparison
    lr_pipeline,     # Trained Linear Regression estimator pipeline
    rf_pipeline,     # Trained Random Forest estimator pipeline
    mlp_pipeline     # Trained Multilayer Perceptron estimator pipeline
    ) :
    """
    Function for testing the Machine Learning estimators on testing data sets.
    
    Note: The estimators used for this MUST be trained with the same features used to test!
    You cannot change (increase, decrease, or swap) the estimators between training and testing.
    
    Note: The number of datapoints for testing CAN be different between training and testing.
    For example, you can train with 100k points and test 500k points.
    
    Returns:
    lr_results         Array of predictions from LinearRegression
    lr_results_delta   Array of differences between predictions from LinearRegression and truth
    rf_results         Array of predictions from RandomForestRegression
    rf_results_delta   Array of differences between predictions from RandomForestRegression and truth
    mlp_results        Array of predictions from MLPRegression
    mlp_results_delta  Array of differences between predictions from MLPRegression and truth
    """
    
    # Initializes ML results variables
    lr_results, lr_results_delta   = False, False
    rf_results, rf_results_delta   = False, False
    mlp_results, mlp_results_delta = False, False
    
    # Tests each estimator and assigns outcomes to results
    if lr_pipeline:
        lr_results, lr_results_delta   = Test_Estimator(lr_pipeline, X_test, y_test)
    
    if rf_pipeline:
        rf_results, rf_results_delta   = Test_Estimator(rf_pipeline, X_test, y_test)
    
    if mlp_pipeline:
        mlp_results, mlp_results_delta = Test_Estimator(mlp_pipeline, X_test, y_test)
    
    return lr_results, lr_results_delta, rf_results, rf_results_delta, mlp_results, mlp_results_delta