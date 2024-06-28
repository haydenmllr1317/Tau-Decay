THIS FILE EXPLAIN EACH DATASET IN DECAY-DATA FOLDER
//
Data Set File Name:
Tau-Decay/Decay-Data/FOM_Studies_genApproach_30May2020_total.root
Looked At?:
Yes
Info:
Likely same data as other sets, 2612 decays at nominal mass. Contained information that allowed us to work in the Upsilon rest frame, and thus reconstruct invariant mass of the combined tau leptons.
//
Data Set File Name:
CUT_15GeVMass_JanIdea.root
Looked At?:
Yes, Gavin had some success here?
Success here with reconstruction data in the tau refereance frames, however could not to for upsilon masses since some data was 
put in as -99.
Info:
This file has 4300+ events with the “upsilon” mass is set to 15 GeV. All the pions have been rotated into the desired frame and you can access that information via the toUse branches.
//
Data Set File Name:
CUT_START_Upsilon_nomMass_Total_5May2020.root
Looked At?:
This is the original set of data we looked at (many variations, highlighted in docs).
Info:
This file has 2612 events with (we believe, given the title), nominal "upsilon" mass. Neutral pions in the decay have been excluded.
//
Data Set File Name:
CUT_nomMass_total_JanIdea.root
Looked At?:
Hayden looking at now.
Info:
This file has 1726 events at seemingly nominal upsilon mass. Might contain vectors rotated into desired frame which would be labeled as "toUse".
//
Data Set File Name:
CUT_nomMass_from_Otto.root
Looked At?:
Gavin looking at now
Info:
This file has around 1000 nominal upsilon mass events. Might contain vectors rotated into desired frame which would be labeled as "toUse".
//
Data Set File Name:
cartesian_upsilon_taus_GENSIM_withBP_*.root, where * \in {1-10}
Looked At?:
Yes, Hayden added them into one file and they all seemed junk, but might be worth another look more closely.
Info:
These are ten files that have 15Gev upsilons (and decay products), with no neutral pions.
File 1: 300+ events
File 2: 339 events
File 3: 356 events
File 4: 320 events
File 5: 317 events
File 6: 333 events
File 7: 334 events
File 8: 352 events
File 9: 377 events
File 10: 311 events
I think it's 3381 total
These are full of -99 fields rahhhhhh so don't use them.
//
Data Set File Name:
GOOD_UpsilonToTauTau_3prong_m15_part2_cartesian_upsilon_taus_UpsilonToTauTau_PUPoissonAve20_102X_upgrade2018_realistic_v18_3prong_m15_miniaod_TOTAL.root
Looked At?:
No
Info:
This is a flat file with 15GeV upsilons with no neutral pions, it contains roughly 2600 events. I would check that this is different data from the other data with 2600 events. There is no evidence that they are similar, the numbers are just the same so worth checking.
