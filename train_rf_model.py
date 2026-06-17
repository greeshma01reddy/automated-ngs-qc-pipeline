import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import os

print(" Initializing RNA Structural Motif Machine Learning Engine...")

# Create results directory if it doesn't exist
os.makedirs("results", exist_ok=True)

# 1. Generate synthetic structural data for training
# Simulating features of RNA motifs (e.g., Hairpins, Internal Loops, Bulges)
np.random.seed(42)
n_samples = 1200

data = {
    'Loop_Length': np.random.randint(3, 15, size=n_samples),
    'Stem_Length': np.random.randint(4, 20, size=n_samples),
    'Free_Energy_DeltaG': np.random.uniform(-35.0, -5.0, size=n_samples),
    'GC_Content': np.random.uniform(0.3, 0.8, size=n_samples),
    'MFE_Probability': np.random.uniform(0.4, 0.99, size=n_samples),
    'Bulge_Count': np.random.choice([0, 1, 2], size=n_samples, p=[0.7, 0.2, 0.1])
}

df = pd.DataFrame(data)

# Simulate a target label: 1 for Stable Functional Motif, 0 for Unstable/Non-functional
# Making the classification dependent on Free Energy, GC Content, and Stem Length
df['Target_Motif'] = ((df['Free_Energy_DeltaG'] < -20) & (df['GC_Content'] > 0.5) & (df['Stem_Length'] > 8)).astype(int)

# Split features and target
X = df.drop(columns=['Target_Motif'])
y = df['Target_Motif']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 2. Train the Tuned Random Forest Classifier
rf = RandomForestClassifier(n_estimators=150, max_depth=8, random_state=42)
rf.fit(X_train, y_train)

accuracy = rf.score(X_test, y_test)
print(f"✅ Model Training Complete. Test Set Accuracy: {accuracy * 100:.1f}%")

# 3. Extract Feature Importances
importances = rf.feature_importances_
feature_imp_df = pd.DataFrame({
    'Feature': X.columns,
    'Importance': importances
}).sort_values(by='Importance', ascending=False)

# Save feature importances to a CSV file for the front-end dashboard to read
feature_imp_df.to_csv("results/rf_feature_importances.csv", index=False)
print("💾 Top structural feature mappings exported safely to 'results/rf_feature_importances.csv'")
