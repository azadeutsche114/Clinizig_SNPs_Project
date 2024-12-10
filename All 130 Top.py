import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import RFE
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.linear_model import LassoCV
import xgboost as xgb
import os

# Load the Shortlist.ST dataset
shortlist_st = pd.read_csv("~/Desktop/SNP Analysis/Processed Data/Shortlist ST.csv")

# Create a list to store the file paths for chunks
chunk_paths = [f"/home/hp/Desktop/SNP Analysis/Processed Data/Final_SNP/chunk{i}.csv" for i in range(1, 113)]

# Process each chunk
for i, chunk_path in enumerate(chunk_paths):
    # Load the data chunk
    data = pd.read_csv(chunk_path)
    data.index = ["X_" + str(pos) for pos in data['Pos']]
    data = data.drop(columns=['Pos']).T
    data['no_of_seeds'] = shortlist_st['Average']
    
    # Separate features and response variable
    X = data.iloc[:, :-1]
    y = data['no_of_seeds']
    
    # Remove features with low standard deviation
    X = X.loc[:, X.std(axis=0) >= 0.1]
    
    # 1. Stabilize Feature Importance using Random Forest
    rf_model = RandomForestRegressor(n_estimators=1000, random_state=123)
    rf_model.fit(X, y)
    
    # Extract feature importance
    importance_rf = pd.Series(rf_model.feature_importances_, index=X.columns)
    top_snps_rf = importance_rf.nlargest(100).index.tolist()
    
    # 2. Stability Selection using Recursive Feature Elimination (RFE) with Random Forest
    rfe_model = RFE(rf_model, n_features_to_select=50, step=1)
    rfe_model.fit(X, y)
    best_features = list(X.columns[rfe_model.support_])
    
    # 3. Feature Selection using XGBoost
    dtrain = xgb.DMatrix(data=X.values, label=y.values, feature_names=X.columns.tolist())
    xgb_params = {"objective": "reg:squarederror", "max_depth": 6, "eta": 0.3}
    xgb_model = xgb.train(xgb_params, dtrain, num_boost_round=100)
    importance_xgb = xgb_model.get_score(importance_type='weight')
    top_features_xgb = sorted(importance_xgb, key=importance_xgb.get, reverse=True)[:5]
    
    # 4. Lasso Neural Network for Feature Selection
    lasso_model = LassoCV(cv=5, random_state=123).fit(X, y)
    lasso_selected_features = list(X.columns[lasso_model.coef_ != 0])
    
    # 5. Recursive Feature Elimination (RFE) with Random Forest (additional step sizes)
    rfe_control = RFE(rf_model, n_features_to_select=30, step=1)
    rfe_control.fit(X, y)
    top_snps_rfe = list(X.columns[rfe_control.support_])
    
    # Ensemble of Features
    fa = set(top_snps_rf).intersection(best_features, top_features_xgb, lasso_selected_features, top_snps_rfe)
    if not fa:
        fb = set(top_snps_rf).intersection(best_features, top_features_xgb, top_snps_rfe)
    else:
        fb = fa
    if not fb:
        final_selected_features = set(top_snps_rf).intersection(best_features, top_snps_rfe)
    else:
        final_selected_features = fb
    
    # Combine feature importances
    feature_scores = pd.DataFrame({
        'RandomForest': [importance_rf.get(f, np.nan) for f in final_selected_features],
        'RFE': [rfe_model.ranking_[list(X.columns).index(f)] if f in X.columns else np.nan for f in final_selected_features],
        'XGBoost': [importance_xgb.get(f, np.nan) for f in final_selected_features]
    }, index=final_selected_features)
    
    feature_scores['mean_rank'] = feature_scores.mean(axis=1, skipna=True)
    best_feature = feature_scores['mean_rank'].idxmin()
    
    # Output the important SNPs
    output_df = pd.DataFrame({'Important SNP': [best_feature]})
    output_path = f"/home/hp/Desktop/SNP Analysis/Output/chunk{i+1}.csv"
    output_df.to_csv(output_path, index=False)
    
    # Status update
    print(f"Processed chunk {i+1}")
    print(best_feature)
