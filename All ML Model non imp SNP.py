import pandas as pd
import numpy as np
from sklearn.linear_model import RidgeCV, LassoCV
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.metrics import mean_squared_error, roc_curve, auc
from sklearn.model_selection import cross_val_score
import xgboost as xgb
import matplotlib.pyplot as plt
import seaborn as sns
import vcf

# Reading data
assessions_df = pd.read_excel("Desktop/SNP Analysis/Accessions shortlisted/Assessions.xlsx", header=None, skiprows=1)
access_no = np.unique(assessions_df.iloc[0, :].astype(int))

st_df = pd.read_csv("~/Desktop/SNP Analysis/Accessions shortlisted/Supplementary Table S1.csv")

# Load VCF data
vcf_reader = vcf.Reader(filename='/media/hp/SVM/1001_GENOME_PROJECT/1001 genomes chunk file/Chunk/vcf_chunk_1.vcf')
vcf_meta = pd.DataFrame(vcf_reader.metadata)

# Extract Genotype data
gq_data = []
dp_data = []
for record in vcf_reader:
    gq_data.append([sample.data.GQ for sample in record.samples])
    dp_data.append([sample.data.DP for sample in record.samples])

gq_data = np.array(gq_data).T
dp_data = np.array(dp_data).T

# Filter based on Accessions
gq_data = gq_data[:, np.isin(vcf_reader.samples, access_no)]
dp_data = dp_data[:, np.isin(vcf_reader.samples, access_no)]

# Replace NAs with 0
gq_data[np.isnan(gq_data)] = 0
dp_data[np.isnan(dp_data)] = 0

# Compute means and standard deviations
gq_mn = np.mean(gq_data, axis=1)
gq_sd = np.std(gq_data, axis=1)
dp_mn = np.mean(dp_data, axis=1)
dp_sd = np.std(dp_data, axis=1)

# Prepare DataFrame
df = pd.DataFrame({
    'GQ_mn': gq_mn,
    'GQ_sd': gq_sd,
    'DP_mn': dp_mn,
    'DP_sd': dp_sd,
    'mean_GQ_squared': gq_mn**2,
    'mean_DT_squared': dp_mn**2,
    'interaction_term': gq_mn * dp_mn,
    'R_1': st_df['R.1']  # Assuming 'R.1' is the response variable in ST
})

# Ridge Regression
ridge_model = RidgeCV(store_cv_values=True)
ridge_model.fit(df[['GQ_mn', 'DP_mn', 'GQ_sd', 'DP_sd', 'mean_GQ_squared', 'mean_DT_squared', 'interaction_term']], df['R_1'])
ridge_predictions = ridge_model.predict(df[['GQ_mn', 'DP_mn', 'GQ_sd', 'DP_sd', 'mean_GQ_squared', 'mean_DT_squared', 'interaction_term']])
ridge_rmse = np.sqrt(mean_squared_error(df['R_1'], ridge_predictions))
print(f"Ridge RMSE: {ridge_rmse}")

# Lasso Regression
lasso_model = LassoCV()
lasso_model.fit(df[['GQ_mn', 'DP_mn', 'GQ_sd', 'DP_sd', 'mean_GQ_squared', 'mean_DT_squared', 'interaction_term']], df['R_1'])
lasso_predictions = lasso_model.predict(df[['GQ_mn', 'DP_mn', 'GQ_sd', 'DP_sd', 'mean_GQ_squared', 'mean_DT_squared', 'interaction_term']])
lasso_rmse = np.sqrt(mean_squared_error(df['R_1'], lasso_predictions))
print(f"Lasso RMSE: {lasso_rmse}")

# Support Vector Regression
svr_model = SVR(kernel='rbf')
svr_model.fit(df[['GQ_mn', 'DP_mn', 'GQ_sd', 'DP_sd', 'mean_GQ_squared', 'mean_DT_squared', 'interaction_term']], df['R_1'])
svr_predictions = svr_model.predict(df[['GQ_mn', 'DP_mn', 'GQ_sd', 'DP_sd', 'mean_GQ_squared', 'mean_DT_squared', 'interaction_term']])
svr_rmse = np.sqrt(mean_squared_error(df['R_1'], svr_predictions))
print(f"SVR RMSE: {svr_rmse}")

# K-Nearest Neighbors Regressor
knn_model = KNeighborsRegressor()
knn_model.fit(df[['GQ_mn', 'DP_mn', 'GQ_sd', 'DP_sd', 'mean_GQ_squared', 'mean_DT_squared', 'interaction_term']], df['R_1'])
knn_predictions = knn_model.predict(df[['GQ_mn', 'DP_mn', 'GQ_sd', 'DP_sd', 'mean_GQ_squared', 'mean_DT_squared', 'interaction_term']])
knn_rmse = np.sqrt(mean_squared_error(df['R_1'], knn_predictions))
print(f"KNN RMSE: {knn_rmse}")

# Decision Tree Regressor
dtr_model = DecisionTreeRegressor()
dtr_model.fit(df[['GQ_mn', 'DP_mn', 'GQ_sd', 'DP_sd', 'mean_GQ_squared', 'mean_DT_squared', 'interaction_term']], df['R_1'])
dtr_predictions = dtr_model.predict(df[['GQ_mn', 'DP_mn', 'GQ_sd', 'DP_sd', 'mean_GQ_squared', 'mean_DT_squared', 'interaction_term']])
dtr_rmse = np.sqrt(mean_squared_error(df['R_1'], dtr_predictions))
print(f"DTR RMSE: {dtr_rmse}")

# Random Forest Regressor
rfr_model = RandomForestRegressor()
rfr_model.fit(df[['GQ_mn', 'DP_mn', 'GQ_sd', 'DP_sd', 'mean_GQ_squared', 'mean_DT_squared', 'interaction_term']], df['R_1'])
rfr_predictions = rfr_model.predict(df[['GQ_mn', 'DP_mn', 'GQ_sd', 'DP_sd', 'mean_GQ_squared', 'mean_DT_squared', 'interaction_term']])
rfr_rmse = np.sqrt(mean_squared_error(df['R_1'], rfr_predictions))
print(f"Random Forest RMSE: {rfr_rmse}")

# Gradient Boosting Regressor
gbm_model = GradientBoostingRegressor(n_estimators=100, max_depth=3, learning_rate=0.1)
gbm_model.fit(df[['GQ_mn', 'DP_mn', 'GQ_sd', 'DP_sd', 'mean_GQ_squared', 'mean_DT_squared', 'interaction_term']], df['R_1'])
gbm_predictions = gbm_model.predict(df[['GQ_mn', 'DP_mn', 'GQ_sd', 'DP_sd', 'mean_GQ_squared', 'mean_DT_squared', 'interaction_term']])
gbm_rmse = np.sqrt(mean_squared_error(df['R_1'], gbm_predictions))
print(f"GBM RMSE: {gbm_rmse}")

# XGBoost Regressor
xgb_model = xgb.XGBRegressor(n_estimators=100, max_depth=3, learning_rate=0.1)
xgb_model.fit(df[['GQ_mn', 'DP_mn', 'GQ_sd', 'DP_sd', 'mean_GQ_squared', 'mean_DT_squared', 'interaction_term']], df['R_1'])
xgb_predictions = xgb_model.predict(df[['GQ_mn', 'DP_mn', 'GQ_sd', 'DP_sd', 'mean_GQ_squared', 'mean_DT_squared', 'interaction_term']])
xgb_rmse = np.sqrt(mean_squared_error(df['R_1'], xgb_predictions))
print(f"XGBoost RMSE: {xgb_rmse}")

# ROC Curves (for binary classification)
threshold = np.median(df['R_1'])
binary_r1 = (df['R_1'] >= threshold).astype(int)

# Calculate ROC curves for each model
def plot_roc(model_predictions, model_name):
    fpr, tpr, _ = roc_curve(binary_r1, model_predictions)
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, label=f'{model_name} (AUC = {roc_auc:.2f})')

plt.figure()
plot_roc(ridge_predictions, 'Ridge')
plot_roc(lasso_predictions, 'Lasso')
plot_roc(svr_predictions, 'SVR')
plot_roc(knn_predictions, 'KNN')
plot_roc(dtr_predictions, 'Decision Tree')
plot_roc(rfr_predictions, 'Random Forest')
plot_roc(gbm_predictions, 'Gradient Boosting')
plot_roc(xgb_predictions, 'XGBoost')
plt.plot([0, 1], [0, 1], linestyle='--', color='gray')
plt.legend(loc='lower right')
plt.title("ROC Curves for Regression Models")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.show()
