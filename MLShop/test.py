#!/usr/bin/env python3 
# Pengxuan Zhu 
# zhupx99@icloud.com

import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score
import numpy as np
import pandas as pd 


def get_data():
    bkg_list = {
        "WW": {
            "path": "../Events/BKG_WW.csv",
            "xsect": 195.2
        },
        "TauTau": {
            "path": "../Events/BKG_TauTau.csv",
            "xsect": 50.2 
        },
        "Wmuv": {
            "path": "../Events/BKG_Wmuv.csv",
            "xsect": 4.9
        },
        "MuMuvv": {
            "path": "../Events/BKG_MuMuvv.csv",
            "xsect": 56.9
        },
        "ZZ": {
            "path": "../Events/BKG_ZZ.csv",
            "xsect": 22.7
        },
        "Zll": {
            "path": "../Events/BKG_Zll.csv",
            "xsect": 22.5
        },
        "Zh": {
            "path": "../Events/BKG_Zh.csv",
            "xsect": 3.5
        }
    }
    sg_list = {
        "path": "../Events/SGN_122_45.csv",
        "xsect": 1.5
    }

    labelwmRC = [
        "pVa.px", "pVa.py", "pVa.pz", "pVa.E", 
        "pVb.px", "pVb.py", "pVb.pz", "pVb.E", 
        "pMiss.px", "pMiss.py", "pMiss.pz", "pMiss.E", 
        "pIa_O.px", "pIa_O.py", "pIa_O.pz", "pIa_O.E", 
        "pIb_O.px", "pIb_O.py", "pIb_O.pz", "pIb_O.E", 
        "pIa_A.px", "pIa_A.py", "pIa_A.pz", "pIa_A.E", 
        "pIb_A.px", "pIb_A.py", "pIb_A.pz", "pIb_A.E", 
        "pIa_B.px", "pIa_B.py", "pIa_B.pz", "pIa_B.E", 
        "pIb_B.px", "pIb_B.py", "pIb_B.pz", "pIb_B.E", 
        "pIa_C.px", "pIa_C.py", "pIa_C.pz", "pIa_C.E", 
        "pIb_C.px", "pIb_C.py", "pIb_C.pz", "pIb_C.E", 
        "mRCmin", "mRC_B", "mRC_C", "mRCmax", "mRCLSP", "mLSPmax", 
        "mRecoil", "mVV", "mInv", "dRVaVb", 
        "dRVaIa_O", "dRVaIa_A", "dRVaIa_B", "dRVaIa_C", 
        "dRVbIb_O", "dRVbIb_A", "dRVbIb_B", "dRVbIb_C", 
        "dRIaIb_O", "dRIaIb_A", "dRIaIb_B", "dRIaIb_C", 
        "ctheta_pO", "ctheta_pA", "chteta_pB", "chteta_pC", 
        "chteta_pMax", "ctheta_mO", "ctheta_mA", "chteta_mB", "chteta_mC", "chteta_mMax", ]

    labelwmRC1 = [
        "pVa.px", "pVa.py", "pVa.pz", "pVa.E", 
        "pVb.px", "pVb.py", "pVb.pz", "pVb.E", 
        "mRCmin", "mRCmax", "mRCLSP", "mLSPmax", 
        "mRecoil", "mVV", "mInv", "dRVaVb", ]


    labelwomRC = [
        "pVa.px", "pVa.py", "pVa.pz", "pVa.E", 
        "pVb.px", "pVb.py", "pVb.pz", "pVb.E", 
        "pMiss.px", "pMiss.py", "pMiss.pz", "pMiss.E", 
    ]

    all_features = []
    all_labels = []
    all_sample_weights = []
    
    # Loop through each CSV file and weight
    for proc, bkg in bkg_list.items():
        print(f"Loading {bkg['path']}")
        df = pd.read_csv(bkg['path'])
        df = df.dropna()

        df['label'] = 0

        # Extract features and labels
        features = df[labelwmRC1].values
        labels = df["label"].values
        # features = df.drop(columns=[label_column]).values
        # labels = df[label_column].values
        
        # Apply the weight to all samples in the current file
        sample_weights = np.full(labels.shape, bkg['xsect'])
        
        # Append to the full dataset
        all_features.append(features)
        all_labels.append(labels)
        all_sample_weights.append(sample_weights)
    
    df = pd.read_csv(sg_list['path'])
    df = df.dropna()
    df['label'] = 1
    
    features    = df[labelwmRC1].values
    labels      = df['label'].values
    sample_weights = np.full(labels.shape, sg_list['xsect'])
    
    all_features.append(features)
    all_labels.append(labels)
    all_sample_weights.append(sample_weights) 
       
    # Combine all data
    all_features = np.vstack(all_features)
    all_labels = np.hstack(all_labels)
    all_sample_weights = np.hstack(all_sample_weights)
    
    
    # Normalize features to the range [0, 1]
    min_vals = all_features.min(axis=0)
    max_vals = all_features.max(axis=0)
    print(min_vals, max_vals)
    normalized_features = (all_features - min_vals) / (max_vals - min_vals)

    # Split into training and testing sets
    X_train, X_test, y_train, y_test, train_weights, test_weights = train_test_split(
        normalized_features, all_labels, all_sample_weights,
        test_size=0.2, random_state=42
    )
    
    return X_train, X_test, y_train, y_test, train_weights, test_weights


# Step 1: Generate synthetic data
def generate_data(n_samples=10000):
    np.random.seed(42)
    # Background: Random data with noise
    bg_features = np.random.normal(0, 1, size=(n_samples, 4))  # 4 features
    bg_labels = np.zeros((n_samples, 1))
    
    # Signal: Slightly different distribution
    signal_features = np.random.normal(1, 1, size=(n_samples, 4))  # 4 features
    signal_labels = np.ones((n_samples, 1))
    
    # Combine and shuffle
    features = np.vstack((bg_features, signal_features))
    labels = np.vstack((bg_labels, signal_labels))
    
    indices = np.arange(len(labels))
    np.random.shuffle(indices)
    return features[indices], labels[indices]


# Step 2: Prepare data
# features, labels = get_data()
# print(features, labels)
X_train, X_test, y_train, y_test, train_weights, test_weights = get_data()

# Convert to PyTorch tensors
X_train = torch.tensor(X_train, dtype=torch.float32)
# y_train = torch.tensor(y_train, dtype=torch.float32)
X_test = torch.tensor(X_test, dtype=torch.float32)
# y_test = torch.tensor(y_test, dtype=torch.float32)
y_train = torch.tensor(y_train, dtype=torch.float32).view(-1, 1)
y_test = torch.tensor(y_test, dtype=torch.float32).view(-1, 1)

# Step 3: Define the model
# class SimpleNN(nn.Module):

#     def __init__(self, input_size):
#         super(SimpleNN, self).__init__()
#         self.fc1 = nn.Linear(input_size, 64)
#         self.relu = nn.ReLU()
#         self.fc2 = nn.Linear(64, 32)
#         self.fc3 = nn.Linear(32, 1)
#         self.sigmoid = nn.Sigmoid()
    
#     def forward(self, x):
#         x = self.relu(self.fc1(x))
#         x = self.relu(self.fc2(x))
#         x = self.sigmoid(self.fc3(x))
#         return x

class SimpleNN(nn.Module):
    def __init__(self, input_size):
        super(SimpleNN, self).__init__()
        self.fc1 = nn.Linear(input_size, 64)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(64, 32)
        self.fc3 = nn.Linear(32, 1)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = self.relu(self.fc1(x))
        x = self.relu(self.fc2(x))
        x = self.sigmoid(self.fc3(x))  # Ensures outputs are in [0, 1]
        return x

# Instantiate the model
input_size = X_train.shape[1]
model = SimpleNN(input_size=input_size)
# model = SimpleNN(input_size=4)

# Step 4: Define loss function and optimizer
criterion = nn.BCELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

# Step 5: Train the model
epochs = 20
batch_size = 64
# for epoch in range(epochs):
#     model.train()
#     permutation = torch.randperm(X_train.size()[0])
#     for i in range(0, X_train.size()[0], batch_size):
#         indices = permutation[i:i + batch_size]
#         batch_X, batch_y = X_train[indices], y_train[indices]
        
#         # Forward pass
#         outputs = model(batch_X).squeeze()
#         loss = criterion(outputs, batch_y.squeeze())
        
#         # Backward pass and optimization
#         optimizer.zero_grad()
#         loss.backward()
#         optimizer.step()
    
#     # Print loss for the epoch
#     print(f"Epoch {epoch+1}/{epochs}, Loss: {loss.item():.4f}")

# Training loop
for epoch in range(epochs):
    model.train()
    permutation = torch.randperm(X_train.size()[0])
    for i in range(0, X_train.size()[0], batch_size):
        indices = permutation[i:i + batch_size]
        batch_X, batch_y = X_train[indices], y_train[indices]

        # Forward pass
        outputs = model(batch_X)  # Already includes sigmoid
        loss = criterion(outputs, batch_y)

        # Backward pass and optimization
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    # Print loss for the epoch
    print(f"Epoch {epoch+1}/{epochs}, Loss: {loss.item():.4f}")

# Step 6: Evaluate the model
model.eval()
with torch.no_grad():
    y_pred = model(X_test).squeeze()
    y_pred_labels = (y_pred > 0.5).float()
    accuracy = accuracy_score(y_test, y_pred_labels)
    roc_auc = roc_auc_score(y_test, y_pred)
    print(f"Accuracy: {accuracy:.4f}, ROC-AUC: {roc_auc:.4f}")
