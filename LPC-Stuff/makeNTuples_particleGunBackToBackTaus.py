#! /usr/bin/env python

import ROOT
import sys
import math
from DataFormats.FWLite import Events, Handle
from collections import OrderedDict
from array import array
from ROOT import TLorentzVector

def isAncestor(a,p) :
    if a == p :
        return True
    for i in xrange(0,p.numberOfMothers()) :
        if isAncestor(a,p.mother(i)) :
            return True
    return False


def computeCosine (Vx, Vy, Vz,
                   Wx, Wy, Wz):

    Vnorm = math.sqrt(Vx*Vx + Vy*Vy + Vz*Vz)
    Wnorm = math.sqrt(Wx*Wx + Wy*Wy + Wz*Wz)
    VdotW = Vx*Wx + Vy*Wy + Vz*Wz

    if Vnorm > 0. and Wnorm > 0.:
        cosAlpha = VdotW / (Vnorm * Wnorm)
    else:
        cosAlpha = -99.

    return cosAlpha

def has_mcEvent_match(g_elec, reco_ele, delta_r_cut, charge_label):
    min_delta_r = 9.9
    delta_r = 10.

    for ii, reco in enumerate(reco_ele):
        if reco.pt() > 0.1:
            reco_vec = ROOT.TLorentzVector(reco.px(), reco.py(), reco.pz(), reco.energy())
            delta_r = g_elec.DeltaR(reco_vec)
            if delta_r < min_delta_r:
                min_delta_r = delta_r
                if delta_r < delta_r_cut:
                    #return True
                    if not abs(charge_label):
                        return (delta_r, ii)
                    elif charge_label*reco.charge() >0:
                        return (delta_r, ii)
    
    return (False, 100)


# Make VarParsing object
# https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideAboutPythonConfigFile#VarParsing_Example
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.parseArguments()
print options

# Events takes either
# - single file name
# - list of file names
# - VarParsing options

# use Varparsing object
events = Events (options)

# Generated stuff
handlePruned  = Handle ("std::vector<reco::GenParticle>")
#labelPruned = ("prunedGenParticles")
labelPruned = ("genParticles")

# Create histograms, etc.
ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.gROOT.SetStyle('Plain') # white background

# ntupla
branches = [
# Generated
#             'g_upsilon_M',
#             'g_upsilon_pt',
#             'g_upsilon_eta',
#             'g_upsilon_phi',
            'g_tau_lead_pt',
            'g_tau_lead_eta',
            'g_tau_lead_phi',
            'g_tau_sublead_pt',
            'g_tau_sublead_eta',
            'g_tau_sublead_phi',
            ]


file_out = ROOT.TFile('nTuples_particleGunBacktoBackTaus.root', 'recreate')
file_out.cd()

ntuple   = ROOT.TNtuple('tree', 'tree', ':'.join(branches))

nTot = 0

verb = False
#verb = True

print events.size()
# loop over events
for event in events:
    leps_mydecay = []
    print 'Processing event: %i...'%(nTot)

    tofill   = OrderedDict(zip(branches, [-99.]*len(branches)))
    # Generated stuff
    event.getByLabel (labelPruned, handlePruned)
    pruned_gen_particles = handlePruned.product()
 

    gen_taum = [pp for pp in pruned_gen_particles if pp.pdgId() == -15]
    
        
    gen_taup = [pp for pp in pruned_gen_particles if pp.pdgId() == 15]
   
    foundtaum = False
    for p in gen_taum:
        gtaum = p
        foundtaum = True
        print "p.px() is:", p.px()
        print "p.py() is:", p.py()
        print "p.pz() is:", p.pz()
        pxSq = p.px()**2
        pySq = p.py()**2
        pzSq = p.pz()**2
        pMagSq = pxSq = pySq + pzSq
        pMag  = math.sqrt(pMagSq)
        print "pxSq is:", pxSq
        print "pySq is:", pySq
        print "pzSq is:", pzSq
        print "pMagSq is:", pMagSq
        print "pMag is:", pMag
        print "p.energy() is:", p.energy()
        print "p.pt() is:", p.pt()
        print "p.eta() is:", p.eta()
        print "p.phi() is:", p.phi()
        print "p.mass() is:", p.mass()
        print '######'
    foundtaup = False    
    for p in gen_taup:
        gtaup = p
        foundtaup = True
        print "p.px() is:", p.px()
        print "p.py() is:", p.py()
        print "p.pz() is:", p.pz()
        pxSq = p.px()**2
        pySq = p.py()**2
        pzSq = p.pz()**2
        pMagSq = pxSq = pySq + pzSq
        pMag  = math.sqrt(pMagSq)
        print "pxSq is:", pxSq
        print "pySq is:", pySq
        print "pzSq is:", pzSq
        print "pMagSq is:", pMagSq
        print "pMag is:", pMag
        print "p.energy() is:", p.energy()
        print "p.pt() is:", p.pt()
        print "p.eta() is:", p.eta()
        print "p.phi() is:", p.phi()
        print "p.mass() is:", p.mass()
        print '######'
    if (foundtaum and foundtaup):
        leps_mydecay.append(gtaum)
        leps_mydecay.append(gtaup)
        
        gtau_lead = gtaum if gtaum.pt() > gtaup.pt() else gtaup
        gtau_sublead = gtaup if gtaum.pt() > gtaup.pt() else gtaum
        

        tau_lead_lv  = TLorentzVector(gtau_lead.px(), gtau_lead.py(), gtau_lead.pz(), gtau_lead.energy())
        tau_sublead_lv  = TLorentzVector(gtau_sublead.px(), gtau_sublead.py(), gtau_sublead.pz(), gtau_sublead.energy())
   
       # fill ntuple with gen quantities
      

        tofill['g_tau_lead_pt'] =tau_lead_lv.Pt()
        tofill['g_tau_lead_eta']=tau_lead_lv.Eta()
        tofill['g_tau_lead_phi']=tau_lead_lv.Phi()

        tofill['g_tau_sublead_pt'] =tau_sublead_lv.Pt()
        tofill['g_tau_sublead_eta']=tau_sublead_lv.Eta()
        tofill['g_tau_sublead_phi']=tau_sublead_lv.Phi()

    
        ntuple.Fill(array('f',tofill.values()))
 
        nTot += 1

# make a canvas, draw, and save it
file_out.cd()
ntuple.Write()



file_out.Close()
