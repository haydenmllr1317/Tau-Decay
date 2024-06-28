import pickle
from fast_histogram import histogram1d
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# pkl_taup = '/Users/gavinpitt/Desktop/taup_file.pkl'
pkl_taup = '/isilon/export/home/hdmiller/fold_empty_attempt/taup_data.pkl'
# pkl_taum = '/Users/gavinpitt/Desktop/taum_file.pkl'
pkl_taum = '/isilon/export/home/hdmiller/fold_empty_attempt/taum_data.pkl'
# pkl_pi = '/Users/gavinpitt/Desktop/pi_file.pkl'
pkl_pi = '/isilon/export/home/hdmiller/fold_empty_attempt/pi_data.pkl'
# pkl_neutrino = '/Users/gavinpitt/Desktop/neutrino_file.pkl'
# pkl_antineutrino = '/Users/gavinpitt/Desktop/antineutrino_file.pkl'
# pkl_upsilon = '/Users/gavinpitt/Desktop/upsilon_file.pkl'

df_taup = pd.read_pickle(pkl_taup)
df_taum = pd.read_pickle(pkl_taum)
df_pi = pd.read_pickle(pkl_pi)
# df_neutrino = pd.read_pickle(pkl_neutrino)
# df_antineutrino = pd.read_pickle(pkl_antineutrino)
# df_upsilon = pd.read_pickle(pkl_upsilon)

def convert_lists(df):
    for column in df.keys():
        print(column)
        df[column] = df[column].apply(lambda x: x[0] if len(x) > 0 else None)
    
    return df

def create_momentums(df):
    df['px'] = df['pt']*np.cos(df['phi'])
    df['py'] = df['pt']*np.sin(df['phi'])
    df['pz'] = df['pt']*np.sinh(df['eta'])
    
    return df


df_taup = convert_lists(df_taup)
df_taum = convert_lists(df_taum)
# df_pi = convert_lists(df_pi)
# df_neutrino = convert_lists(df_neutrino)
# df_antineutrino = convert_lists(df_antineutrino)
# df_upsilon = convert_lists(df_upsilon)

# df_minus1 = df_pi[['minus1_pt', 'minus1_eta', 'minus1_phi']]
# df_minus2 = df_pi[['minus2_pt', 'minus2_eta', 'minus2_phi']]
# df_minus3 = df_pi[['minus3_pt', 'minus3_eta', 'minus3_phi']]
# df_plus1 = df_pi[['plus1_pt', 'plus1_eta', 'plus1_phi']]
# df_plus2 = df_pi[['plus2_pt', 'plus2_eta', 'plus2_phi']]
# df_plus3 = df_pi[['plus3_pt', 'plus3_eta', 'plus3_phi']]

# df_minus1 = df_minus1.rename(columns={"minus1_pt": "pt", "minus1_eta": "eta", "minus1_phi": "phi"})
# df_minus2 = df_minus2.rename(columns={"minus2_pt": "pt", "minus2_eta": "eta", "minus2_phi": "phi"})
# df_minus3 = df_minus3.rename(columns={"minus3_pt": "pt", "minus3_eta": "eta", "minus3_phi": "phi"})
# df_plus1 = df_plus1.rename(columns={"plus1_pt": "pt", "plus1_eta": "eta", "plus1_phi": "phi"})
# df_plus2 = df_plus2.rename(columns={"plus2_pt": "pt", "plus2_eta": "eta", "plus2_phi": "phi"})
# df_plus3 = df_plus3.rename(columns={"plus3_pt": "pt", "plus3_eta": "eta", "plus3_phi": "phi"})


df_taup = create_momentums(df_taup)
df_taum = create_momentums(df_taum)
# df_neutrino = create_momentums(df_neutrino)
# df_antineutrino = create_momentums(df_antineutrino)
# df_upsilon = create_momentums(df_upsilon)

# df_minus1 = create_momentums(df_minus1)
# df_minus2 = create_momentums(df_minus2)
# df_minus3 = create_momentums(df_minus3)
# df_plus1 = create_momentums(df_plus1)
# df_plus2 = create_momentums(df_plus2)
# df_plus3 = create_momentums(df_plus3)

# pi_mass_df = pd.DataFrame({'mass': [0.13957] * 2612})

tau_inv = []

tau_inv =  np.square(np.sqrt(np.square(df_taup['px'])+np.square(df_taup['py'])+np.square(df_taup['pz'])+np.square(df_taup['mass'])) + \
    np.sqrt(np.square(df_taum['px'])+np.square(df_taum['py'])+np.square(df_taum['pz'])+np.square(df_taum['mass']))) - \
    (np.square(df_taup['px']+df_taum['px'])+np.square(df_taup['py']+df_taum['py'])+np.square(df_taup['pz']+df_taum['pz']))

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
bins = 20  # Number of bins
# range_min = 0
range_min = 5
# range_max = 180
range_max = 15



# Compute the histogram
# hist1 = histogram1d(pi_inv, bins=bins, range=[range_min, range_max])
# hist2 = histogram1d(pi_neutrino_inv, bins=bins, range=[range_min, range_max])
hist3 = histogram1d(tau_inv, bins=bins, range=[range_min, range_max])

# Create bin edges
bin_edges = np.linspace(range_min, range_max, bins + 1)

# Plotting the histogram
plt.figure(figsize=(10, 6))
# plt.bar(bin_edges[:-1], hist1, width=np.diff(bin_edges), edgecolor='black')
# plt.bar(bin_edges[:-1], hist2, width=np.diff(bin_edges), edgecolor='red')
plt.bar(bin_edges[:-1], hist3, width=np.diff(bin_edges), edgecolor='blue')
# plt.title('Invariant Mass of Lepton Decays w/o Neutrino Momentum')
plt.title('Invariant Mass of tt Pairs')
plt.xlabel('Invariant Mass')
plt.ylabel('Frequency')
# plt.show()
plt.savefig("taus_plot_try1.png")

# print(df_taup)