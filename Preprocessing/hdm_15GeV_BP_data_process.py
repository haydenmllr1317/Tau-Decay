import pickle
from fast_histogram import histogram1d
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df_taup = {}
df_taum = {}
df_pis = {}

df_taum[1] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP1_taum_data.pkl')
df_taup[1] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP1_taup_data.pkl')
df_pis[1] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP1_pis_data.pkl')
df_taum[2] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP2_taum_data.pkl')
df_taup[2] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP2_taup_data.pkl')
df_pis[2] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP2_pis_data.pkl')
df_taum[3] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP3_taum_data.pkl')
df_taup[3] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP3_taup_data.pkl')
df_pis[3] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP3_pis_data.pkl')
df_taum[4] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP4_taum_data.pkl')
df_taup[4] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP4_taup_data.pkl')
df_pis[4] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP4_pis_data.pkl')
df_taum[5] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP5_taum_data.pkl')
df_taup[5] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP5_taup_data.pkl')
df_pis[5] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP5_pis_data.pkl')
df_taum[6] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP6_taum_data.pkl')
df_taup[6] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP6_taup_data.pkl')
df_pis[6] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP6_pis_data.pkl')
df_taum[7] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP7_taum_data.pkl')
df_taup[7] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP7_taup_data.pkl')
df_pis[7] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP7_pis_data.pkl')
df_taum[8] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP8_taum_data.pkl')
df_taup[8] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP8_taup_data.pkl')
df_pis[8] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP8_pis_data.pkl')
df_taum[9] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP9_taum_data.pkl')
df_taup[9] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP9_taup_data.pkl')
df_pis[9] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP9_pis_data.pkl')
df_taum[10] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP10_taum_data.pkl')
df_taup[10] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP10_taup_data.pkl')
df_pis[10] = pd.read_pickle('/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data/15GeV_BP10_pis_data.pkl')




df_combined_taup = pd.concat([df_taup[1], df_taup[2], df_taup[3], df_taup[4], df_taup[5], df_taup[6], df_taup[7], df_taup[8], df_taup[9], df_taup[10]], axis=0)
df_combined_taum = pd.concat([df_taum[1], df_taum[2], df_taum[3], df_taum[4], df_taum[5], df_taum[6], df_taum[7], df_taum[8], df_taum[9], df_taum[10]], axis=0)
df_combined_pis = pd.concat([df_pis[1], df_pis[2], df_pis[3], df_pis[4], df_pis[5], df_pis[6], df_pis[7], df_pis[8], df_pis[9], df_pis[10]], axis=0)

def convert_lists(df):
    for column in df.keys():
        print(column)
        df[column] = df[column].apply(lambda x: x[0] if len(x) > 0 else None)
    
    return df

# def create_momentums(df):
#     df['px'] = df['pt']*np.cos(df['phi'])
#     df['py'] = df['pt']*np.sin(df['phi'])
#     df['pz'] = df['pt']*np.sinh(df['eta'])
    
#     return df


df_final_taup = convert_lists(df_combined_taup)
df_final_taum = convert_lists(df_combined_taum)
df_final_pis = convert_lists(df_combined_pis)
# df_neutrino = convert_lists(df_neutrino)
# df_antineutrino = convert_lists(df_antineutrino)
# df_upsilon = convert_lists(df_upsilon)

# df_minus1 = df_final_pis[['minus1_pt', 'minus1_eta', 'minus1_phi']]
# df_minus2 = df_final_pis[['minus2_pt', 'minus2_eta', 'minus2_phi']]
# df_minus3 = df_final_pis[['minus3_pt', 'minus3_eta', 'minus3_phi']]
# df_plus1 = df_final_pis[['plus1_pt', 'plus1_eta', 'plus1_phi']]
# df_plus2 = df_final_pis[['plus2_pt', 'plus2_eta', 'plus2_phi']]
# df_plus3 = df_final_pis[['plus3_pt', 'plus3_eta', 'plus3_phi']]

# df_minus1 = df_minus1.rename(columns={"minus1_pt": "pt", "minus1_eta": "eta", "minus1_phi": "phi"})
# df_minus2 = df_minus2.rename(columns={"minus2_pt": "pt", "minus2_eta": "eta", "minus2_phi": "phi"})
# df_minus3 = df_minus3.rename(columns={"minus3_pt": "pt", "minus3_eta": "eta", "minus3_phi": "phi"})
# df_plus1 = df_plus1.rename(columns={"plus1_pt": "pt", "plus1_eta": "eta", "plus1_phi": "phi"})
# df_plus2 = df_plus2.rename(columns={"plus2_pt": "pt", "plus2_eta": "eta", "plus2_phi": "phi"})
# df_plus3 = df_plus3.rename(columns={"plus3_pt": "pt", "plus3_eta": "eta", "plus3_phi": "phi"})


# df_taup_rotated = create_momentums(df_final_taup)
# df_taum_rotated = create_momentums(df_final_taum)
# df_neutrino = create_momentums(df_neutrino)
# df_antineutrino = create_momentums(df_antineutrino)
# df_upsilon = create_momentums(df_upsilon)

# df_minus1 = create_momentums(df_minus1)
# df_minus2 = create_momentums(df_minus2)
# df_minus3 = create_momentums(df_minus3)
# df_plus1 = create_momentums(df_plus1)
# df_plus2 = create_momentums(df_plus2)
# df_plus3 = create_momentums(df_plus3)

# pi_mass_df = pd.DataFrame({'mass': [0.13957] * len(df_minus1)})

# tau_inv = []

# tau_inv =  np.square(np.sqrt(np.square(df_taup_rotated['px'])+np.square(df_taup_rotated['py'])+np.square(df_taup_rotated['pz'])+np.square(df_taup_rotated['mass'])) + \
#     np.sqrt(np.square(df_taum_rotated['px'])+np.square(df_taum_rotated['py'])+np.square(df_taum_rotated['pz'])+np.square(df_taum_rotated['mass']))) - \
#     (np.square(df_taup_rotated['px']+df_taum_rotated['px'])+np.square(df_taup_rotated['py']+df_taum_rotated['py'])+np.square(df_taup_rotated['pz']+df_taum_rotated['pz']))

# pi_neutrino_inv = []

# pi_neutrino_inv = np.square(np.sqrt(np.square(df_minus1['px'])+np.square(df_minus1['py'])+np.square(df_minus1['pz'])+np.square(pi_mass_df['mass'])) + \
#     np.sqrt(np.square(df_minus2['px'])+np.square(df_minus2['py'])+np.square(df_minus2['pz'])+np.square(pi_mass_df['mass'])) + \
#     np.sqrt(np.square(df_minus3['px'])+np.square(df_minus3['py'])+np.square(df_minus3['pz'])+np.square(pi_mass_df['mass'])) + \
#     np.sqrt(np.square(df_plus1['px'])+np.square(df_plus1['py'])+np.square(df_plus1['pz'])+np.square(pi_mass_df['mass'])) + \
#     np.sqrt(np.square(df_plus2['px'])+np.square(df_plus2['py'])+np.square(df_plus2['pz'])+np.square(pi_mass_df['mass'])) + \
#     np.sqrt(np.square(df_plus3['px'])+np.square(df_plus3['py'])+np.square(df_plus3['pz'])+np.square(pi_mass_df['mass'])) + \
#     np.sqrt(np.square(df_neutrino['px'])+np.square(df_neutrino['py'])+np.square(df_neutrino['pz'])) + \
#     np.sqrt(np.square(df_antineutrino['px'])+np.square(df_antineutrino['py'])+np.square(df_antineutrino['pz']))) - \
#     np.square(df_minus1['px']+df_minus2['px']+df_minus3['px']+df_plus1['px']+df_plus2['px']+df_plus3['px']+df_neutrino['px']+df_antineutrino['px']) - \
#     np.square(df_minus1['py']+df_minus2['py']+df_minus3['py']+df_plus1['py']+df_plus2['py']+df_plus3['py']+df_neutrino['py']+df_antineutrino['py']) - \
#     np.square(df_minus1['pz']+df_minus2['pz']+df_minus3['pz']+df_plus1['pz']+df_plus2['pz']+df_plus3['pz']+df_neutrino['pz']+df_antineutrino['pz'])
    
# pi_inv = []

# pi_inv = np.square(np.sqrt(np.square(df_minus1['px'])+np.square(df_minus1['py'])+np.square(df_minus1['pz'])+np.square(pi_mass_df['mass'])) + \
#     np.sqrt(np.square(df_minus2['px'])+np.square(df_minus2['py'])+np.square(df_minus2['pz'])+np.square(pi_mass_df['mass'])) + \
#     np.sqrt(np.square(df_minus3['px'])+np.square(df_minus3['py'])+np.square(df_minus3['pz'])+np.square(pi_mass_df['mass'])) + \
#     np.sqrt(np.square(df_plus1['px'])+np.square(df_plus1['py'])+np.square(df_plus1['pz'])+np.square(pi_mass_df['mass'])) + \
#     np.sqrt(np.square(df_plus2['px'])+np.square(df_plus2['py'])+np.square(df_plus2['pz'])+np.square(pi_mass_df['mass'])) + \
#     np.sqrt(np.square(df_plus3['px'])+np.square(df_plus3['py'])+np.square(df_plus3['pz'])+np.square(pi_mass_df['mass']))) - \
#     np.square(df_minus1['px']+df_minus2['px']+df_minus3['px']+df_plus1['px']+df_plus2['px']+df_plus3['px']) - \
#     np.square(df_minus1['py']+df_minus2['py']+df_minus3['py']+df_plus1['py']+df_plus2['py']+df_plus3['py']) - \
#     np.square(df_minus1['pz']+df_minus2['pz']+df_minus3['pz']+df_plus1['pz']+df_plus2['pz']+df_plus3['pz'])

# upsilon_inv = []

# upsilon_inv = np.square(df_upsilon['mass'])







# Parameters for the histogram
# bins = 80  # Number of bins
bins = 100  # Number of bins
# range_min = 0
range_min = -200
# range_max = 180
range_max = 100



# Compute the histogram
print(df_final_pis['mass'])
# print(len(df_final_taup['mass']))
# print(len(df_final_taum['mass']))
print(df_final_taum['neu_pt'])

hist1 = histogram1d(df_final_pis['mass'], bins=bins, range=[range_min, range_max])
hist2 = histogram1d(df_final_taup['mass'], bins=bins, range=[range_min, range_max])
hist3 = histogram1d(df_final_taum['mass'], bins=bins, range=[range_min, range_max])

# Create bin edges
bin_edges = np.linspace(range_min, range_max, bins + 1)

# Plotting the histogram
plt.figure(figsize=(10, 6))
plt.bar(bin_edges[:-1], hist1, width=np.diff(bin_edges), edgecolor='black')
plt.bar(bin_edges[:-1], hist2, width=np.diff(bin_edges), edgecolor='red')
plt.bar(bin_edges[:-1], hist3, width=np.diff(bin_edges), edgecolor='blue')
# plt.title('Invariant Mass of Lepton Decays w/o Neutrino Momentum')
plt.title('Invariant Mass of tt Pairs (Blue) vs Six Pions (Black)')
plt.xlabel('Invariant Mass')
plt.ylabel('Frequency')
# plt.show()
plt.savefig("taus_vs_pions_15GeV_BP.png")

# print(df_taup)