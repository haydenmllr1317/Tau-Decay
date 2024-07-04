import pickle
from fast_histogram import histogram1d
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
"""
pkl_taup = '/Users/gavinpitt/Desktop/data pickles/taup_file_2.pkl'
pkl_taum = '/Users/gavinpitt/Desktop/data pickles/taum_file_2.pkl'
pkl_pi = '/Users/gavinpitt/Desktop/data pickles/pi_file_2.pkl'
pkl_neutrino = '/Users/gavinpitt/Desktop/data pickles/neutrino_file_2.pkl'
pkl_antineutrino = '/Users/gavinpitt/Desktop/data pickles/antineutrino_file_2.pkl'
pkl_upsilon = '/Users/gavinpitt/Desktop/data pickles/upsilon_file_2.pkl'
pkl_tau = '/Users/gavinpitt/Desktop/data pickles/tau_file_2.pkl'

df_taup = pd.read_pickle(pkl_taup)
df_taum = pd.read_pickle(pkl_taum)
df_pi = pd.read_pickle(pkl_pi)
df_neutrino = pd.read_pickle(pkl_neutrino)
df_antineutrino = pd.read_pickle(pkl_antineutrino)
df_upsilon = pd.read_pickle(pkl_upsilon)
df_tau = pd.read_pickle(pkl_tau)
"""
pkl_toUse = '/Users/gavinpitt/Desktop/data pickles/toUse_file_otto.pkl'

df_toUse = pd.read_pickle(pkl_toUse)
pkl_otto_taup  = pd.read_pickle('/Users/gavinpitt/Desktop/data pickles/toUse_file_otto_taup.pkl')

print(df_toUse)
['local_taup_lv_mass','local_pi_p_lv1_phi', 'local_neu_lv_pt', 'local_pi_p_lv1_pt','local_pi_p_lv2_pt', 'local_pi_p_lv3_pt', 
'local_neu_lv_theta', 'local_pi_p_lv1_theta', 'local_pi_p_lv2_theta', 'local_pi_p_lv3_theta', 'local_pi_p_lv2_phi', 
'local_pi_p_lv3_phi','local_neu_lv_phi']

def convert_lists(df):
    for column in df.keys():
        df[column] = df[column].apply(lambda x: x[0] if len(x) > 0 else None)
    
    return df

def create_momentums(df):
    df['px'] = df['pt']*np.cos(df['phi'])
    df['py'] = df['pt']*np.sin(df['phi'])
    df['pz'] = df['pt']*np.sinh(df['eta'])
    
    return df

df_toUse = convert_lists(df_toUse)

px_pi_1 = df_toUse['local_pi_p_lv1_pt']*np.cos(df_toUse['local_pi_p_lv1_phi'])
px_pi_2 = df_toUse['local_pi_p_lv2_pt']*np.cos(df_toUse['local_pi_p_lv2_phi'])
px_pi_3 = df_toUse['local_pi_p_lv3_pt']*np.cos(df_toUse['local_pi_p_lv3_phi'])

py_pi_1 = df_toUse['local_pi_p_lv1_pt']*np.sin(df_toUse['local_pi_p_lv1_phi'])
py_pi_2 = df_toUse['local_pi_p_lv2_pt']*np.sin(df_toUse['local_pi_p_lv2_phi'])
py_pi_3 = df_toUse['local_pi_p_lv3_pt']*np.sin(df_toUse['local_pi_p_lv3_phi'])

eta_pi_1 = np.negative(np.log(np.tan(df_toUse['local_pi_p_lv1_theta'] / 2)))
eta_pi_2 = np.negative(np.log(np.tan(df_toUse['local_pi_p_lv2_theta'] / 2)))
eta_pi_3 = np.negative(np.log(np.tan(df_toUse['local_pi_p_lv3_theta'] / 2)))

pz_pi_1 = df_toUse['local_pi_p_lv1_pt']*np.sinh(eta_pi_1)
pz_pi_2 = df_toUse['local_pi_p_lv2_pt']*np.sinh(eta_pi_2)
pz_pi_3 = df_toUse['local_pi_p_lv3_pt']*np.sinh(eta_pi_3)

eta_antineu = np.negative(np.log(np.tan(df_toUse['local_antineu_lv_theta'] / 2)))

px_antineu = df_toUse['local_antineu_lv_pt']*np.cos(df_toUse['local_antineu_lv_phi'])
py_antineu = df_toUse['local_antineu_lv_pt']*np.sin(df_toUse['local_antineu_lv_phi'])
pz_antineu = df_toUse['local_antineu_lv_pt']*np.sinh(eta_antineu)

#########

px_pi_1m = df_toUse['local_pi_m_lv1_pt']*np.cos(df_toUse['local_pi_m_lv1_phi'])
px_pi_2m = df_toUse['local_pi_m_lv2_pt']*np.cos(df_toUse['local_pi_m_lv2_phi'])
px_pi_3m = df_toUse['local_pi_m_lv3_pt']*np.cos(df_toUse['local_pi_m_lv3_phi'])

py_pi_1m = df_toUse['local_pi_m_lv1_pt']*np.sin(df_toUse['local_pi_m_lv1_phi'])
py_pi_2m = df_toUse['local_pi_m_lv2_pt']*np.sin(df_toUse['local_pi_m_lv2_phi'])
py_pi_3m = df_toUse['local_pi_m_lv3_pt']*np.sin(df_toUse['local_pi_m_lv3_phi'])

eta_pi_1m = np.negative(np.log(np.tan(df_toUse['local_pi_m_lv1_theta'] / 2)))
eta_pi_2m = np.negative(np.log(np.tan(df_toUse['local_pi_m_lv2_theta'] / 2)))
eta_pi_3m = np.negative(np.log(np.tan(df_toUse['local_pi_m_lv3_theta'] / 2)))

pz_pi_1m = df_toUse['local_pi_m_lv1_pt']*np.sinh(eta_pi_1m)
pz_pi_2m = df_toUse['local_pi_m_lv2_pt']*np.sinh(eta_pi_2m)
pz_pi_3m = df_toUse['local_pi_m_lv3_pt']*np.sinh(eta_pi_3m)

eta_neu = np.negative(np.log(np.tan(df_toUse['local_neu_lv_theta'] / 2)))

px_neu = df_toUse['local_neu_lv_pt']*np.cos(df_toUse['local_neu_lv_phi'])
py_neu = df_toUse['local_neu_lv_pt']*np.sin(df_toUse['local_neu_lv_phi'])
pz_neu = df_toUse['local_neu_lv_pt']*np.sinh(eta_neu)


mass_pi = 0.13957



inv_plus = np.square(np.sqrt(np.square(px_pi_1)+np.square(py_pi_1)+np.square(pz_pi_1)+np.square(mass_pi)) + \
                 np.sqrt(np.square(px_pi_2)+np.square(py_pi_2)+np.square(pz_pi_2)+np.square(mass_pi)) + \
                np.sqrt(np.square(px_pi_3)+np.square(py_pi_3)+np.square(pz_pi_3)+np.square(mass_pi)) + \
                np.sqrt(np.square(px_antineu)+np.square(py_antineu)+np.square(pz_antineu))) - \
                np.square(px_pi_1+px_pi_2+px_pi_3+px_antineu) - \
                np.square(py_pi_1+py_pi_2+py_pi_3+py_antineu) - \
                np.square(pz_pi_1+pz_pi_2+pz_pi_3+pz_antineu)

inv_plus_no_neu = np.square(np.sqrt(np.square(px_pi_1)+np.square(py_pi_1)+np.square(pz_pi_1)+np.square(mass_pi)) + \
                 np.sqrt(np.square(px_pi_2)+np.square(py_pi_2)+np.square(pz_pi_2)+np.square(mass_pi)) + \
                np.sqrt(np.square(px_pi_3)+np.square(py_pi_3)+np.square(pz_pi_3)+np.square(mass_pi))) - \
                np.square(px_pi_1+px_pi_2+px_pi_3) - \
                np.square(py_pi_1+py_pi_2+py_pi_3) - \
                np.square(pz_pi_1+pz_pi_2+pz_pi_3)

inv_minus = np.square(np.sqrt(np.square(px_pi_1m)+np.square(py_pi_1m)+np.square(pz_pi_1m)+np.square(mass_pi)) + \
                 np.sqrt(np.square(px_pi_2m)+np.square(py_pi_2m)+np.square(pz_pi_2m)+np.square(mass_pi)) + \
                np.sqrt(np.square(px_pi_3m)+np.square(py_pi_3m)+np.square(pz_pi_3m)+np.square(mass_pi)) + \
                np.sqrt(np.square(px_neu)+np.square(py_neu)+np.square(pz_neu))) - \
                np.square(px_pi_1m+px_pi_2m+px_pi_3m+px_neu) - \
                np.square(py_pi_1m+py_pi_2m+py_pi_3m+py_neu) - \
                np.square(pz_pi_1m+pz_pi_2m+pz_pi_3m+pz_neu)

inv_minus_no_neu = np.square(np.sqrt(np.square(px_pi_1m)+np.square(py_pi_1m)+np.square(pz_pi_1m)+np.square(mass_pi)) + \
                 np.sqrt(np.square(px_pi_2m)+np.square(py_pi_2m)+np.square(pz_pi_2m)+np.square(mass_pi)) + \
                np.sqrt(np.square(px_pi_3m)+np.square(py_pi_3m)+np.square(pz_pi_3m)+np.square(mass_pi))) - \
                np.square(px_pi_1m+px_pi_2m+px_pi_3m) - \
                np.square(py_pi_1m+py_pi_2m+py_pi_3m) - \
                np.square(pz_pi_1m+pz_pi_2m+pz_pi_3m)

inv_taum = np.square(df_toUse['local_taum_lv_mass'])

inv_taup = np.square(df_toUse['local_taup_lv_mass'])

# Parameters for the histogram
bins = 80  # Number of bins
range_min = 0
range_max = 5



# Compute the histogram
hist1 = histogram1d(inv_plus, bins=bins, range=[range_min, range_max])
hist2 = histogram1d(inv_minus, bins=bins, range=[range_min, range_max])
hist3 = histogram1d(inv_taum, bins=bins, range=[range_min, range_max])
hist4 = histogram1d(inv_taup, bins=bins, range=[range_min, range_max])
hist5 = histogram1d(inv_plus_no_neu, bins=bins, range=[range_min, range_max])
hist6 = histogram1d(inv_minus_no_neu, bins=bins, range=[range_min, range_max])

# Create bin edges
bin_edges = np.linspace(range_min, range_max, bins + 1)

# Plotting the histogram

plt.figure(figsize=(10, 6))
plt.bar(bin_edges[:-1], hist1, width=np.diff(bin_edges), edgecolor='black')
plt.title('Invariant Mass (squared) of Antitau Decay Products')
plt.xlabel('Invariant Mass (squared GeV)')
plt.ylabel('Frequency')
plt.savefig('antitau_decay_sep_otto.png')
#plt.show()

plt.figure(figsize=(10, 6))
plt.bar(bin_edges[:-1], hist2, width=np.diff(bin_edges), edgecolor='black')
plt.title('Invariant Mass (squared) of Tau Decay Products')
plt.xlabel('Invariant Mass (squared GeV)')
plt.ylabel('Frequency')
plt.savefig('tau_decay_sep_otto.png')
#plt.show()

plt.figure(figsize=(10, 6))
plt.bar(bin_edges[:-1], hist3, width=np.diff(bin_edges), edgecolor='black')
plt.title('Invariant Mass (squared) of Tau')
plt.xlabel('Invariant Mass (squared GeV)')
plt.ylabel('Frequency')
plt.savefig('tau_sep_otto.png')
#plt.show()

plt.figure(figsize=(10, 6))
plt.bar(bin_edges[:-1], hist4, width=np.diff(bin_edges), edgecolor='black')
plt.title('Invariant Mass (squared) of Antitau')
plt.xlabel('Invariant Mass (squared GeV)')
plt.ylabel('Frequency')
plt.savefig('antitau_sep_otto.png')
#plt.show()

plt.figure(figsize=(10, 6))
plt.bar(bin_edges[:-1], hist5, width=np.diff(bin_edges), edgecolor='black')
plt.title('Invariant Mass (squared) of Antitau Decay Products w/o Antineutrino')
plt.xlabel('Invariant Mass (squared GeV)')
plt.ylabel('Frequency')
plt.savefig('sep_no_anti_otto.png')
#plt.show()


plt.figure(figsize=(10, 6))
plt.bar(bin_edges[:-1], hist6, width=np.diff(bin_edges), edgecolor='black')
plt.title('Invariant Mass (squared) of Tau Decay Products w/o Neutrino')
plt.xlabel('Invariant Mass (squared GeV)')
plt.ylabel('Frequency')
plt.savefig('sep_no_neutrino_otto.png')
#plt.show()
"""

def convert_lists(df):
    for column in df.keys():
        df[column] = df[column].apply(lambda x: x[0] if len(x) > 0 else None)
    
    return df

def create_momentums(df):
    df['px'] = df['pt']*np.cos(df['phi'])
    df['py'] = df['pt']*np.sin(df['phi'])
    df['pz'] = df['pt']*np.sinh(df['eta'])
    
    return df


df_taup = convert_lists(df_taup)
df_taum = convert_lists(df_taum)
df_pi = convert_lists(df_pi)
df_neutrino = convert_lists(df_neutrino)
df_antineutrino = convert_lists(df_antineutrino)
df_upsilon = convert_lists(df_upsilon)

df_minus1 = df_pi[['minus1_pt', 'minus1_eta', 'minus1_phi']]
df_minus2 = df_pi[['minus2_pt', 'minus2_eta', 'minus2_phi']]
df_minus3 = df_pi[['minus3_pt', 'minus3_eta', 'minus3_phi']]
df_plus1 = df_pi[['plus1_pt', 'plus1_eta', 'plus1_phi']]
df_plus2 = df_pi[['plus2_pt', 'plus2_eta', 'plus2_phi']]
df_plus3 = df_pi[['plus3_pt', 'plus3_eta', 'plus3_phi']]

df_minus1 = df_minus1.rename(columns={"minus1_pt": "pt", "minus1_eta": "eta", "minus1_phi": "phi"})
df_minus2 = df_minus2.rename(columns={"minus2_pt": "pt", "minus2_eta": "eta", "minus2_phi": "phi"})
df_minus3 = df_minus3.rename(columns={"minus3_pt": "pt", "minus3_eta": "eta", "minus3_phi": "phi"})
df_plus1 = df_plus1.rename(columns={"plus1_pt": "pt", "plus1_eta": "eta", "plus1_phi": "phi"})
df_plus2 = df_plus2.rename(columns={"plus2_pt": "pt", "plus2_eta": "eta", "plus2_phi": "phi"})
df_plus3 = df_plus3.rename(columns={"plus3_pt": "pt", "plus3_eta": "eta", "plus3_phi": "phi"})


df_taup = create_momentums(df_taup)
df_taum = create_momentums(df_taum)
df_neutrino = create_momentums(df_neutrino)
df_antineutrino = create_momentums(df_antineutrino)
df_upsilon = create_momentums(df_upsilon)

df_minus1 = create_momentums(df_minus1)
df_minus2 = create_momentums(df_minus2)
df_minus3 = create_momentums(df_minus3)
df_plus1 = create_momentums(df_plus1)
df_plus2 = create_momentums(df_plus2)
df_plus3 = create_momentums(df_plus3)

pi_mass_df = pd.DataFrame({'mass': [0.13957] * 2612})

tau_inv = []

tau_inv =  np.square(np.sqrt(np.square(df_taup['px'])+np.square(df_taup['py'])+np.square(df_taup['pz'])+np.square(df_taup['mass'])) + \
    np.sqrt(np.square(df_taum['px'])+np.square(df_taum['py'])+np.square(df_taum['pz'])+np.square(df_taum['mass']))) - \
    (np.square(df_taup['px']+df_taum['px'])+np.square(df_taup['py']+df_taum['py'])+np.square(df_taup['pz']+df_taum['pz']))

pi_neutrino_inv = []

pi_neutrino_inv = np.square(np.sqrt(np.square(df_minus1['px'])+np.square(df_minus1['py'])+np.square(df_minus1['pz'])+np.square(pi_mass_df['mass'])) + \
    np.sqrt(np.square(df_minus2['px'])+np.square(df_minus2['py'])+np.square(df_minus2['pz'])+np.square(pi_mass_df['mass'])) + \
    np.sqrt(np.square(df_minus3['px'])+np.square(df_minus3['py'])+np.square(df_minus3['pz'])+np.square(pi_mass_df['mass'])) + \
    np.sqrt(np.square(df_plus1['px'])+np.square(df_plus1['py'])+np.square(df_plus1['pz'])+np.square(pi_mass_df['mass'])) + \
    np.sqrt(np.square(df_plus2['px'])+np.square(df_plus2['py'])+np.square(df_plus2['pz'])+np.square(pi_mass_df['mass'])) + \
    np.sqrt(np.square(df_plus3['px'])+np.square(df_plus3['py'])+np.square(df_plus3['pz'])+np.square(pi_mass_df['mass'])) + \
    np.sqrt(np.square(df_neutrino['px'])+np.square(df_neutrino['py'])+np.square(df_neutrino['pz'])) + \
    np.sqrt(np.square(df_antineutrino['px'])+np.square(df_antineutrino['py'])+np.square(df_antineutrino['pz']))) - \
    np.square(df_minus1['px']+df_minus2['px']+df_minus3['px']+df_plus1['px']+df_plus2['px']+df_plus3['px']+df_neutrino['px']+df_antineutrino['px']) - \
    np.square(df_minus1['py']+df_minus2['py']+df_minus3['py']+df_plus1['py']+df_plus2['py']+df_plus3['py']+df_neutrino['py']+df_antineutrino['py']) - \
    np.square(df_minus1['pz']+df_minus2['pz']+df_minus3['pz']+df_plus1['pz']+df_plus2['pz']+df_plus3['pz']+df_neutrino['pz']+df_antineutrino['pz'])
    
pi_inv = []

pi_inv = np.square(np.sqrt(np.square(df_minus1['px'])+np.square(df_minus1['py'])+np.square(df_minus1['pz'])+np.square(pi_mass_df['mass'])) + \
    np.sqrt(np.square(df_minus2['px'])+np.square(df_minus2['py'])+np.square(df_minus2['pz'])+np.square(pi_mass_df['mass'])) + \
    np.sqrt(np.square(df_minus3['px'])+np.square(df_minus3['py'])+np.square(df_minus3['pz'])+np.square(pi_mass_df['mass'])) + \
    np.sqrt(np.square(df_plus1['px'])+np.square(df_plus1['py'])+np.square(df_plus1['pz'])+np.square(pi_mass_df['mass'])) + \
    np.sqrt(np.square(df_plus2['px'])+np.square(df_plus2['py'])+np.square(df_plus2['pz'])+np.square(pi_mass_df['mass'])) + \
    np.sqrt(np.square(df_plus3['px'])+np.square(df_plus3['py'])+np.square(df_plus3['pz'])+np.square(pi_mass_df['mass']))) - \
    np.square(df_minus1['px']+df_minus2['px']+df_minus3['px']+df_plus1['px']+df_plus2['px']+df_plus3['px']) - \
    np.square(df_minus1['py']+df_minus2['py']+df_minus3['py']+df_plus1['py']+df_plus2['py']+df_plus3['py']) - \
    np.square(df_minus1['pz']+df_minus2['pz']+df_minus3['pz']+df_plus1['pz']+df_plus2['pz']+df_plus3['pz'])

upsilon_inv = []

upsilon_inv = np.square(df_upsilon['mass'])



"""




