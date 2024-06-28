import os
import pandas as pd
import numpy as np
import warnings
from coffea.nanoevents import NanoEventsFactory, PFNanoAODSchema
from tqdm import tqdm
import argparse
import awkward as ak
import uproot


# Path to your ROOT file
filepath = '/isilon/export/home/hdmiller/Tau-Decay/Decay-Data/cartesian_upsilon_taus_GENSIM_withBP_10.root'

# Open the ROOT file
# with uproot.open(filepath) as file:
#     # List all keys in the ROOT file
#     keys = file.keys()
#     print(keys)

def main(datapath, savepath, datatype_taup, datatype_taum, datatype_pi, filetype):
    # TODO: should we actually care about this warning? 
    warnings.filterwarnings("ignore", message="Found duplicate branch")
    if filetype == '.h5':
        load_h5()
    else:
        # print("trying to load root file")
        df_taup, df_taum, df_pi = load_root(datapath)
        # load_root(datapath)
        # print("went through root file load process")

    # print("trying to save output files")
    print("saving taup file")
    df_taup.to_pickle(savepath + "/" + datatype_taup + ".pkl")
    print("saving taum file")
    df_taup.to_pickle(savepath + "/" + datatype_taum + ".pkl")
    print("saving pi file")
    df_taup.to_pickle(savepath + "/" + datatype_pi + ".pkl")
    print("end of main")
    return

def load_h5():
    return

def load_root(filepath):
    data_taup = None
    data_taum = None
    data_pi = None
    # print(filepath)
    if os.path.splitext((filepath))[-1] == ".root" and os.path.isfile(filepath):

        warnings.filterwarnings("ignore")

        ### this is the line that I am having trouble with ###
        #events = NanoEventsFactory.from_root(filepath, schemaclass = PFNanoAODSchema).events()
        events = NanoEventsFactory.from_root(filepath, treepath = "tree", schemaclass = PFNanoAODSchema).events()
        # everything_there = NanoEventsFactory.from_root(filepath, treepath = "tree", schemaclass = PFNanoAODSchema)
        # all_keys = everything_there.schema.all_keys()
        # the_tree = NanoEventsFactory.from_root(filepath, schemaclass = PFNanoAODSchema).tree
        # all_keys = the_tree.all_keys()
        print(len(events))

        # random_indices = np.random.choice(len(events), size=200, replace=False)
        # print(random_indices)

    # j = 0
    # for index in tqdm(random_indices):
    for index in range(len(events)):
    # for index in range(20):
        # j += 1
        properties_taup, property_names_taup, properties_taum, property_names_taum, properties_pi, property_names_pi = process_event_root(events[index:index+1])
        if properties_taup != -1 and data_taup == None: 
            data_taup = {property_name: [] for property_name in property_names_taup}
            for i, prop in enumerate(properties_taup): 
                data_taup[property_names_taup[i]].append(prop)
        if properties_taum != -1 and data_taum == None: 
            data_taum = {property_name: [] for property_name in property_names_taum}
            for k, prop in enumerate(properties_taum): 
                data_taum[property_names_taum[k]].append(prop)
        if properties_pi != -1 and data_pi == None: 
            data_pi = {property_name: [] for property_name in property_names_pi}
            for l, prop in enumerate(properties_pi): 
                data_pi[property_names_pi[l]].append(prop)
        elif properties_taup != -1 and properties_taum != -1 and properties_pi != -1: 
            for i, prop in enumerate(properties_taup): 
                data_taup[property_names_taup[i]].append(prop)
            for k, prop in enumerate(properties_taum): 
                data_taum[property_names_taum[k]].append(prop)
            for l, prop in enumerate(properties_pi): 
                data_pi[property_names_pi[l]].append(prop)
            # if len(data['pt']) == 20: 
            # if j == 200: 
            #     print("200 events reached")
            #     # print("taup data:")
            #     # print(pd.DataFrame.from_dict(data_taup))
            #     # print("taum data:")
            #     # print(pd.DataFrame.from_dict(data_taum))
            #     # print("pi data:")
            #     # print(pd.DataFrame.from_dict(data_pi))
            #     df_taup = pd.DataFrame.from_dict(data_taup)
            #     df_taum = pd.DataFrame.from_dict(data_taum)
            #     df_pi = pd.DataFrame.from_dict(data_pi)
            #     return df_taup, df_taum, df_pi
                        # if len(data['pt']) == 20: 

    df_taup = pd.DataFrame.from_dict(data_taup)
    df_taum = pd.DataFrame.from_dict(data_taum)
    df_pi = pd.DataFrame.from_dict(data_pi)

    # return pd.DataFrame.from_dict(data_taup)
    return df_taup, df_taum, df_pi

def process_event_root(events):
    
    # fields = dir(events)
    # print("Available fields:", fields)
    # Available fields: ['Mask', '__class__', '__delattr__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__get__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__self__', '__self_class__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__thisclass__', '_behavior', '_caches', '_layout', '_numbaview', 'add_kind', 'add_systematic', 'antineutrino', 'behavior', 'caches', 'candMatchPi1Info', 'candMatchPi2Info', 'check', 'check1', 'check2', 'describe_variations', 'explodes_how', 'fields', 'initial', 'layout', 'local', 'mask', 'metadata', 'nbytes', 'ndim', 'neutrino', 'numba_type', 'orig', 'pi', 'slot0', 'slot1', 'slot2', 'slot3', 'slot4', 'slot5', 'slot6', 'slot7', 'slot8', 'slot9', 'systematics', 'taum', 'taup', 'toUse', 'to_list', 'to_numpy', 'tolist', 'type', 'upsilon']

    taup = events.taup
    taum = events.taum
    pi = events.pi

    # eta = [ak.to_numpy(pfs['phi']).flatten()[i] for i in range(len(pfs))]
    # phi = [ak.to_numpy(pfs['eta']).flatten()[i] for i in range(len(pfs))]
    # pt = [ak.to_numpy(pfs['pt']).flatten()[i] for i in range(len(pfs))]

    # properties = [pt, eta, phi]
    properties_taup = []
    properties_taum = []
    properties_pi = []
    # property_names = ["pt", "eta", "phi"]
    property_names_taup = []
    property_names_taum = []
    property_names_pi = []

    # add all other fields
    fields_taup = taup.fields
    fields_taum = taum.fields
    fields_pi = pi.fields

    for field in fields_taup: 
        # if field != "pt" and field != "eta" and field != "phi": 
            properties_taup.append([ak.to_numpy(taup[field]).flatten()[i] for i in range(len(taup))])
            property_names_taup.append(field)
    for field in fields_taum: 
        # if field != "pt" and field != "eta" and field != "phi": 
            properties_taum.append([ak.to_numpy(taum[field]).flatten()[i] for i in range(len(taum))])
            property_names_taum.append(field)
    for field in fields_pi: 
        # if field != "pt" and field != "eta" and field != "phi": 
            properties_pi.append([ak.to_numpy(pi[field]).flatten()[i] for i in range(len(pi))])
            property_names_pi.append(field)
    # print(property_names_taup)
    # print(property_names_taum)
    # print(property_names_pi)
    return properties_taup, property_names_taup, properties_taum, property_names_taum, properties_pi, property_names_pi


print("run began okay")
main(filepath, "/isilon/export/home/hdmiller/Tau-Decay/Pickled-Data", "15GeV_BP10_taup_data", "15GeV_BP10_taum_data", "15GeV_BP10_pis_data", ".root")