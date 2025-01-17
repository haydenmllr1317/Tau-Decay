import pickle
# from fast_histogram import histogram1d
import pandas as pd
# import matplotlib.pyplot as plt
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from sklearn.model_selection import train_test_split
from torch.optim.lr_scheduler import ReduceLROnPlateau



pkl_JanIdea = '/isilon/export/home/hdmiller/cms_work/Tau-Decay/Pickled-Data/RECO-Data/nomMass_JanIdea.pkl'
pkl_START_Upsilon = '/isilon/export/home/hdmiller/cms_work/Tau-Decay/Pickled-Data/RECO-Data/nomMass_START_Upsilon.pkl'
df_JanIdea = pd.read_pickle(pkl_JanIdea)
# print(df_JanIdea.keys())
df_START_Upsilon = pd.read_pickle(pkl_START_Upsilon)
# print(df_START_Upsilon.keys())

df = pd.concat([df_JanIdea, df_START_Upsilon])
# print(df.keys())

# df['local_pi_lv1_pt'] = pd.concat([df['local_pi_m_lv1_pt'], df['local_pi_p_lv1_pt']], axis = 1)
# df['local_pi_lv2_pt'] = pd.concat([df['local_pi_m_lv2_pt'], df['local_pi_p_lv2_pt']], axis = 1)
# df['local_pi_lv3_pt'] = pd.concat([df['local_pi_m_lv3_pt'], df['local_pi_p_lv3_pt']], axis = 1)
# df['local_pi_lv1_phi'] = pd.concat([df['local_pi_m_lv1_phi'], df['local_pi_p_lv1_phi']], axis = 1)
# df['local_pi_lv2_phi'] = pd.concat([df['local_pi_m_lv2_phi'], df['local_pi_p_lv2_phi']], axis = 1)
# df['local_pi_lv3_phi'] = pd.concat([df['local_pi_m_lv3_phi'], df['local_pi_p_lv3_phi']], axis = 1)
# df['local_pi_lv1_theta'] = pd.concat([df['local_pi_m_lv1_theta'], df['local_pi_p_lv1_theta']], axis = 1)
# df['local_pi_lv2_theta'] = pd.concat([df['local_pi_m_lv2_theta'], df['local_pi_p_lv2_theta']], axis = 1)
# df['local_pi_lv3_theta'] = pd.concat([df['local_pi_m_lv3_theta'], df['local_pi_p_lv3_theta']], axis = 1)

# print(df.keys())
def convert_lists(df):
    for column in df.keys():
        # print(column)
        df[column] = df[column].apply(lambda x: x[0] if len(x) > 0 else None)
    
    return df

df = convert_lists(df)
print(len(df))
print(df.shape)

print(df['local_pi_m_lv1_pt'].shape)


data_pi1_pt = pd.concat([df['local_pi_m_lv1_pt'], df['local_pi_p_lv1_pt']], ignore_index=True)
data_pi2_pt = pd.concat([df['local_pi_m_lv2_pt'], df['local_pi_p_lv2_pt']], ignore_index=True)
data_pi3_pt = pd.concat([df['local_pi_m_lv3_pt'], df['local_pi_p_lv3_pt']], ignore_index=True)
data_pi1_phi = pd.concat([df['local_pi_m_lv1_phi'], df['local_pi_p_lv1_phi']], ignore_index=True)
data_pi2_phi = pd.concat([df['local_pi_m_lv2_phi'], df['local_pi_p_lv2_phi']], ignore_index=True)
data_pi3_phi = pd.concat([df['local_pi_m_lv3_phi'], df['local_pi_p_lv3_phi']], ignore_index=True)
data_pi1_theta = pd.concat([df['local_pi_m_lv1_theta'], df['local_pi_p_lv1_theta']], ignore_index=True)
data_pi2_theta = pd.concat([df['local_pi_m_lv2_theta'], df['local_pi_p_lv2_theta']], ignore_index=True)
data_pi3_theta = pd.concat([df['local_pi_m_lv3_theta'], df['local_pi_p_lv3_theta']], ignore_index=True)

data_n_pt = pd.concat([df['local_neu_lv_pt'], df['local_antineu_lv_pt']], ignore_index=True)
data_n_phi = pd.concat([df['local_neu_lv_phi'], df['local_antineu_lv_phi']], ignore_index=True)
data_n_theta = pd.concat([df['local_neu_lv_theta'], df['local_antineu_lv_theta']], ignore_index=True)


# df['local_pi_lv1_pt'] = pd.concat([df['local_pi_m_lv1_pt'], df['local_pi_p_lv1_pt']], ignore_index=True)
# print(df['local_pi_lv1_pt'].shape)
# df['local_pi_lv2_pt'] = pd.concat([df['local_pi_m_lv2_pt'], df['local_pi_p_lv2_pt']], ignore_index=True)
# df['local_pi_lv3_pt'] = pd.concat([df['local_pi_m_lv3_pt'], df['local_pi_p_lv3_pt']], ignore_index=True)
# df['local_pi_lv1_phi'] = pd.concat([df['local_pi_m_lv1_phi'], df['local_pi_p_lv1_phi']], ignore_index=True)
# df['local_pi_lv2_phi'] = pd.concat([df['local_pi_m_lv2_phi'], df['local_pi_p_lv2_phi']], ignore_index=True)
# df['local_pi_lv3_phi'] = pd.concat([df['local_pi_m_lv3_phi'], df['local_pi_p_lv3_phi']], ignore_index=True)
# df['local_pi_lv1_theta'] = pd.concat([df['local_pi_m_lv1_theta'], df['local_pi_p_lv1_theta']], ignore_index=True)
# df['local_pi_lv2_theta'] = pd.concat([df['local_pi_m_lv2_theta'], df['local_pi_p_lv2_theta']], ignore_index=True)
# df['local_pi_lv3_theta'] = pd.concat([df['local_pi_m_lv3_theta'], df['local_pi_p_lv3_theta']], ignore_index=True)

# df['local_n_lv_pt'] = pd.concat([df['local_neu_lv_pt'], df['local_antineu_lv_pt']], ignore_index=True)
# df['local_n_lv_theta'] = pd.concat([df['local_neu_lv_theta'], df['local_antineu_lv_theta']], ignore_index=True)
# df['local_n_lv_phi'] = pd.concat([df['local_neu_lv_phi'], df['local_antineu_lv_phi']], ignore_index=True)


# print(df['local_n_lv_pt'].shape)

data_cos1 = np.cos(data_pi1_phi.values)
data_sin1 = np.sin(data_pi1_phi.values)
data_cos2 = np.cos(data_pi2_phi.values)
data_sin2 = np.sin(data_pi2_phi.values)
data_cos3 = np.cos(data_pi3_phi.values)
data_sin3 = np.sin(data_pi3_phi.values)

all_data = {'pi_lv1_pt': data_pi1_pt, 'pi_lv2_pt': data_pi2_pt, 'pi_lv3_pt': data_pi3_pt, 
          'pi_lv1_cosphi': data_cos1, 'pi_lv2_cosphi': data_cos2, 'pi_lv3_cosphi': data_cos3,
          'pi_lv1_sinphi': data_sin1, 'pi_lv2_sinphi': data_sin2, 'pi_lv3_cosphi': data_sin3,
          'pi_lv1_theta': data_pi1_theta, 'pi_lv2_theta': data_pi2_theta, 'pi_lv3_pt': data_pi3_theta,
          'n_pt': data_n_pt, 'n_theta': data_n_theta}

new_df = pd.DataFrame(all_data)

print(new_df.shape)
print(new_df.keys())

# df['cos(phi)1'] = np.cos(df[['local_pi_m_lv1_phi']].values)
# df['sin(phi)1'] = np.sin(df[['local_pi_m_lv1_phi']].values)
# df['cos(phi)2'] = np.cos(df[['local_pi_m_lv2_phi']].values)
# df['sin(phi)2'] = np.sin(df[['local_pi_m_lv2_phi']].values)
# df['cos(phi)3'] = np.cos(df[['local_pi_m_lv3_phi']].values)
# df['sin(phi)3'] = np.sin(df[['local_pi_m_lv3_phi']].values)

# df['cos(phi)'] = np.cos(df[['local_neu_lv_phi']].values)
# df['sin(phi)'] = np.sin(df[['local_neu_lv_phi']].values)

#print(df['local_taum_lv_mass'])
# print(df.keys())

# X = df[['local_pi_lv1_pt', 'cos(phi)1', 'sin(phi)1', 'local_pi_lv1_theta',
#         'local_pi_lv2_pt', 'cos(phi)2', 'sin(phi)2', 'local_pi_lv2_theta',
#         'local_pi_lv3_pt', 'cos(phi)3', 'sin(phi)3', 'local_pi_lv3_theta']].values
# Y = df[['local_neu_lv_pt', 'local_neu_lv_theta']].values
# # Y = df[['local_neu_lv_pt', 'local_neu_lv_phi', 'local_neu_lv_theta']].values

# # X = df[['local_pi_m_lv1_pt', 'cos(phi)1', 'sin(phi)1', 'local_pi_m_lv1_theta',
# #         'local_pi_m_lv2_pt', 'cos(phi)2', 'sin(phi)2', 'local_pi_m_lv2_theta',
# #         'local_pi_m_lv3_pt', 'cos(phi)3', 'sin(phi)3', 'local_pi_m_lv3_theta']].values

# # X = df[['local_pi_m_lv1_pt', 'local_pi_m_lv1_phi', 'local_pi_m_lv1_theta',
# #         'local_pi_m_lv2_pt', 'local_pi_m_lv2_phi', 'local_pi_m_lv2_theta',
# #         'local_pi_m_lv3_pt', 'local_pi_m_lv3_phi', 'local_pi_m_lv3_theta',
# #         'local_pi_p_lv1_pt', 'local_pi_p_lv1_phi', 'local_pi_p_lv1_theta',
# #         'local_pi_p_lv2_pt', 'local_pi_p_lv2_phi', 'local_pi_p_lv2_theta',
# #         'local_pi_p_lv3_pt', 'local_pi_p_lv3_phi', 'local_pi_p_lv3_theta']].values
# # Y = df[['local_neu_lv_pt', 'local_neu_lv_phi', 'local_neu_lv_theta',
# #        'local_antineu_lv_pt', 'local_antineu_lv_phi', 'local_antineu_lv_theta']].values

# X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

# class CustomDataset(Dataset):
#     def __init__(self, features, labels):
#         self.features = torch.tensor(features, dtype=torch.float32)
#         self.labels = torch.tensor(labels, dtype=torch.float32)
        
#     def __len__(self):
#         return len(self.features)
    
#     def __getitem__(self, idx):
#         return self.features[idx], self.labels[idx]

# train_dataset = CustomDataset(X_train, Y_train)
# test_dataset = CustomDataset(X_test, Y_test)

# train_loader = DataLoader(train_dataset, batch_size=64, shuffle=True)
# test_loader = DataLoader(test_dataset, batch_size=64, shuffle=False)

# class SimpleDNN(nn.Module):
#     def __init__(self, input_dim, output_dim):
#         super(SimpleDNN, self).__init__()
#         self.fc1 = nn.Linear(input_dim, 128)
#         # self.fc2 = nn.Linear(128, 256)
#         # self.fc3 = nn.Linear(256, 256)
#         # self.fc4 = nn.Linear(256, 128)
#         self.fc5 = nn.Linear(128, 64)
#         self.fc6 = nn.Linear(64, 48)
#         self.fc7 = nn.Linear(48, 32)
#         self.fc8 = nn.Linear(32, 24)
#         self.fc9 = nn.Linear(24, 16)
#         self.fc10 = nn.Linear(16, 12)
#         self.fc11 = nn.Linear(12, 8)
#         self.fc12 = nn.Linear(8, 6)
#         self.fc13 = nn.Linear(6, output_dim)
#         # self.dropout_a = nn.Dropout(p=0.5)
#         # self.dropout_b = nn.Dropout(p=0.2)
#         # self.dropout_c = nn.Dropout(p=0.1)
#         # she had 12 layers, up to 2560 neurons, all relu with droppout from 0.3, up to 0.5, and then progressively down to 0.05

#     def forward(self, x):
#         x = F.relu(self.fc1(x))
#         # x = self.dropout_a(x)
#         # x = F.relu(self.fc2(x))
#         # x = self.dropout_a(x)
#         # x = F.relu(self.fc3(x))
#         # x = self.dropout_a(x)
#         # x = F.relu(self.fc4(x))
#         # x = self.dropout_b(x)
#         x = F.relu(self.fc5(x))
#         # x = self.dropout_c(x)
#         x = F.relu(self.fc6(x))
#         x = F.relu(self.fc7(x))
#         x = F.relu(self.fc8(x))
#         x = F.relu(self.fc9(x))
#         x = F.relu(self.fc10(x))
#         x = F.relu(self.fc11(x))
#         x = F.relu(self.fc12(x))
#         x = F.relu(self.fc13(x))

#         return x

# # Instantiate the model
# input_dim = X_train.shape[1]
# output_dim = Y_train.shape[1]
# print(output_dim)
# model = SimpleDNN(input_dim, output_dim)

# # Loss and optimizer
# #criterion = nn.BCELoss()
# criterion = nn.MSELoss()
# optimizer = optim.Adam(model.parameters(), lr=0.001)
# # scheduler = ReduceLROnPlateau(optimizer, mode='min', factor=0.1, patience=5, verbose=True)


# # Training loop
# num_epochs = 50
# model.train()
# for epoch in range(num_epochs):
#     for features, labels in train_loader:
#         # Forward pass
#         outputs = model(features)
#         loss = criterion(outputs.squeeze(), labels)
        
#         # Backward pass and optimization
#         optimizer.zero_grad()
#         loss.backward()
#         optimizer.step()
    
#     print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}')
#     # scheduler.step(loss)


# # Evaluation
# model.eval()
# count = 0
# for name, param in model.named_parameters():
#     if param.requires_grad:
#         #print(f'Layer: {name} | Size: {param.size()} | Number of parameters: {param.numel()}')
#         count = count + param.numel()
#         print("Total param count:" + str(count))
# with torch.no_grad():
#     # correct = 0
#     # total = 0
#     # for features, labels in test_loader:
#     #     outputs = model(features)
#     #     predicted = (outputs.squeeze() > 0.5).float()
#     #     total += labels.size(0)
#     #     correct += (predicted == labels).sum().item()

#     # accuracy = correct / total
#     # print(f'Accuracy of the model on the test set: {accuracy * 100:.2f}%')

#     for features, labels in test_loader:
#         outputs = model(features)
#         loss = criterion(outputs.squeeze(), labels)

#     print(f'Eval Loss: {loss.item():.4f}')
