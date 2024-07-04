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
filepath = '/Users/gavinpitt/Desktop/CERN Workspace/Tau-Decay/Decay-Data/CUT_15GeVMass_JanIdea.root'
# Open the ROOT file
# with uproot.open(filepath) as file:
#     # List all keys in the ROOT file
#     keys = file.keys()
#     print(keys)

def main(datapath, savepath, datatype, filetype):
    # TODO: should we actually care about this warning? 
    warnings.filterwarnings("ignore", message="Found duplicate branch")
    if filetype == '.h5':
        load_h5()
    else:
        # print("trying to load root file")
        df = load_root(datapath)
        # load_root(datapath)
        # print("went through root file load process")

    print("trying to save output file")
    df.to_pickle(savepath + "/" + datatype + ".pkl")
    print("end of main")
    return

def load_h5():
    return

def load_root(filepath):
    data = None
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

    for i in tqdm(range(len(events))):
        properties, property_names = process_event_root(events[i:i+1])
        if properties != -1 and data == None: 
            data = {property_name: [] for property_name in property_names}
            for i, prop in enumerate(properties): 
                data[property_names[i]].append(prop)
        elif properties != -1: 
            for i, prop in enumerate(properties): 
                data[property_names[i]].append(prop)

    return pd.DataFrame.from_dict(data)
            



def process_event_root(events):
    
    fields = dir(events)

 
    available_columns = ['antineutrino', 'neutrino', 'pi', 'taum', 'taup', 'upsilon']
    

    pfs = events.toUse

    """
    pfs.append(events.fields)
    pfs.append(events.initial)
    pfs.append(events.neutrino)
    pfs.append(events.pi)
    pfs.append(events.taum)
    pfs.append(events.taup)
    pfs.append(events.upsilon)
    """

    

    #eta = [ak.to_numpy(pfs['phi']).flatten()[i] for i in range(len(pfs))]
    #phi = [ak.to_numpy(pfs['eta']).flatten()[i] for i in range(len(pfs))]
    #pt = [ak.to_numpy(pfs['pt']).flatten()[i] for i in range(len(pfs))]

    #properties = [pt, eta, phi]
    #property_names = ["pt", "eta", "phi"]

    properties = []
    property_names = []

    # add all other fields
    fields = pfs.fields
   

    for field in fields: 
            #if field != "pt" and field != "eta" and field != "phi": 
        properties.append([ak.to_numpy(pfs[field]).flatten()[i] for i in range(len(pfs))])
        property_names.append(field)

    return properties, property_names


print("run began okay")
main(filepath, '/Users/gavinpitt/Desktop/data pickles', "JanIdea", ".root")
    