import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data
DF = pd.read_csv("~/Desktop/SNP Analysis/Accessions shortlisted/GQ_DP_final.csv")

# Adding polynomial and interaction terms
DF['mean_GQ_squared'] = DF['GQ_mn']**2
DF['mean_DT_squared'] = DF['DP_mn']**2
DF['interaction_term'] = DF['mean_GQ'] * DF['mean_DT']

# Prepare data for XGBoost for R1, R2, and R3
X = DF[['GQ_mn', 'DP_mn', 'GQ_sd', 'DP_sd', 'mean_GQ_squared', 'mean_DT_squared', 'interaction_term']]
y_R1 = DF['R.1']
y_R2 = DF['R.2']
y_R3 = DF['R.3']

# Define the parameter grid
param_grid = {
    'eta': [0.01, 0.05, 0.1],
    'max_depth': [3, 5, 7],
    'subsample': [0.7, 0.8, 1],
    'colsample_bytree': [0.7, 0.8, 1]
}

# Grid search function
def xgb_grid_search(params, X, y, num_rounds, nfold):
    best_rmse = float('inf')
    best_params = None
    
    # Create XGBoost model
    model = xgb.XGBRegressor(objective="reg:squarederror", n_estimators=num_rounds)
    
    grid_search = GridSearchCV(model, param_grid=params, cv=nfold, scoring='neg_root_mean_squared_error')
    grid_search.fit(X, y)
    
    best_rmse = -grid_search.best_score_  # Negative because GridSearchCV minimizes
    best_params = grid_search.best_params_
    
    return best_params, best_rmse

# Running the model for R1, R2, and R3
max_iterations = 10
rmse_threshold = 0.5
iteration = 0

# Function to plot predictions
def plot_predictions(actual, predicted, model_name):
    plt.figure(figsize=(8, 6))
    plt.scatter(actual, predicted, color='blue', alpha=0.6)
    plt.plot([min(actual), max(actual)], [min(actual), max(actual)], color='red', linestyle='dashed')
    plt.title(f"{model_name} Predictions vs Actual")
    plt.xlabel("Actual")
    plt.ylabel("Predicted")
    plt.xlim(min(actual), max(actual))
    plt.ylim(min(actual), max(actual))
    plt.show()

# Running iterations for R1, R2, and R3
for R, y in zip(['R1', 'R2', 'R3'], [y_R1, y_R2, y_R3]):
    iteration = 0
    while iteration < max_iterations:
        iteration += 1
        print(f"Iteration: {iteration} for {R}")
        
        # Grid search for best parameters
        best_params, best_rmse = xgb_grid_search(param_grid, X, y, num_rounds=100, nfold=5)
        
        # Train final XGBoost model
        model = xgb.XGBRegressor(objective="reg:squarederror", **best_params, n_estimators=100)
        model.fit(X, y)
        
        # Make predictions
        predictions = model.predict(X)
        
        # Compute RMSE
        rmse = np.sqrt(mean_squared_error(y, predictions))
        print(f"Tuned XGBoost RMSE for {R}: {rmse}")
        
        # Check if RMSE meets the threshold
        if rmse < rmse_threshold:
            print(f"Super accuracy achieved for {R}. Stopping...\n")
            break
        
        # Stop if max iterations are reached
        if iteration >= max_iterations:
            print(f"Max iterations reached for {R}. Stopping...\n")
            break
        
    # Plot predictions vs actual values
    plot_predictions(y, predictions, f"Tuned XGBoost for {R}")
    
    # Line plots for actual vs predicted
    plt.figure(figsize=(8, 6))
    plt.plot(y, color='green', linewidth=3, label="Actual")
    plt.plot(predictions, color='red', linewidth=1, label="Predicted")
    plt.title(f"{R} Actual vs Predicted")
    plt.xlabel("Index")
    plt.ylabel("Values")
    plt.legend()
    plt.show()

# If needed, repeat similar steps for each of R1, R2, and R3 or perform any other analysis
