## Data prep

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
import matplotlib.pyplot as plt
import seaborn as sns

# Load your data
df = pd.read_csv('Trimmed_Aerial_Data_draft.csv')

# Create presence/absence column (assuming 'Number of Mantas' indicates presence)
df['presence'] = (df['Number of Mantas'] > 0).astype(int)

# Select features - these are the environmental variables
features = [
    'tavg_c',              # Average temperature
    'prcp_mm_prev_day',    # Previous day's precipitation
    'tide_height_m',       # Tide height
    'tide_stage',          # Tide stage (Ebb/Flood/Slack)
    'lunar_illum',         # Moon illumination
    'nearest_inlet_km'     # Distance to nearest inlet
]

# Convert categorical variables to numerical
df['tide_stage'] = df['tide_stage'].map({
    'Ebb (Outgoing)': 0,
    'Flood (Incoming)': 1,
    'Slack': 0.5
})

# Handle missing values
df = df.fillna(df.median())

# Prepare X (features) and y (target)
X = df[features]
y = df['presence']

##Train test split
# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, random_state=42, stratify=y
)


#Train a random forest model

# Initialize the model
model = RandomForestClassifier(
    n_estimators=100,
    random_state=42,
    class_weight='balanced'  # Important for imbalanced data
)

# Train the model
model.fit(X_train, y_train)

# Make predictions
y_pred = model.predict(X_test)

# Evaluate the model
print(classification_report(y_test, y_pred))

## feature importance

# Plot feature importance
feature_importance = pd.Series(
    model.feature_importances_,
    index=features
).sort_values(ascending=False)

plt.figure(figsize=(10, 6))
sns.barplot(x=feature_importance, y=feature_importance.index)
plt.title('Feature Importance')
plt.tight_layout()
plt.show()


##create prediction map

# Predict probabilities for all points
df['predicted_probability'] = model.predict_proba(X)[:, 1]

# Create a simple map of presence probability
plt.figure(figsize=(12, 8))
scatter = plt.scatter(
    df['Longitude'], 
    df['Latitude'], 
    c=df['predicted_probability'],
    cmap='viridis',
    alpha=0.6
)
plt.colorbar(scatter, label='Predicted Probability of Presence')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Manta Ray Presence Probability')
plt.show()


#model interpretation

# Partial dependence plots
from sklearn.inspection import PartialDependenceDisplay

fig, ax = plt.subplots(figsize=(15, 10))
display = PartialDependenceDisplay.from_estimator(
    model, X, features, 
    kind='both', 
    centered=True,
    pd_linewidth=2,
    ax=ax
)
_ = display.figure_.suptitle('Partial Dependence Plots', fontsize=16)
