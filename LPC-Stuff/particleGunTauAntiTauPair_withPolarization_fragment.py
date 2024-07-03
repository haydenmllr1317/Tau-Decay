import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *
from Configuration.Generator.PSweightsPythia.PythiaPSweightsSettings_cfi import *

generator = cms.EDFilter("Pythia8PtGun", #looked here to get the name right: https://github.com/cms-sw/cmssw/blob/master/GeneratorInterface/Pythia8Interface/plugins/Py8PtGun.cc
#                                 comEnergy = cms.double(13000.0),
#                                 filterEfficiency = cms.untracked.double(1.0),
                                 maxEventsToPrint = cms.untracked.int32(0),
                                 pythiaPylistVerbosity = cms.untracked.int32(0),
                                 pythiaHepMCVerbosity = cms.untracked.bool(False),
                                
  
    PGunParameters = cms.PSet(
       ParticleID = cms.vint32(15), #tau
       AddAntiParticle = cms.bool(True), #make an anti tau too, this tau will be made back to back, I verified this is so
       MinPhi = cms.double(-3.14159265359),
       MaxPhi = cms.double(3.14159265359),
       MinPt = cms.double(2.0), #these numbers from an email from Greg from 9 October 2019
       MaxPt = cms.double(10.0),
       MinEta = cms.double(-2.5), #to be within tracker acceptance 
       MaxEta = cms.double(2.5)
    ),#close PGunParameters  stuff
    
     PythiaParameters = cms.PSet(

    pythia8CommonSettingsBlock,
    pythia8CP5SettingsBlock,
    pythia8PSweightsSettingsBlock,
    processParameters = cms.vstring(
      'Main:timesAllowErrors    = 10000',
      #'HiggsSM:all=true',
      #'25:m0 = 125.0',
      #'25:onMode = off',
      #'25:addChannel = 1  1.00   103   22   553',
      #'Bottomonium:all = on',
      # 'Bottomonium:gg2bbbar(3S1)[3S1(1)]g    = on,off,off',
#       'Bottomonium:gg2bbbar(3S1)[3S1(1)]gm   = on,off,off',
#       'Bottomonium:gg2bbbar(3S1)[3S1(8)]g    = on,off,off',
#       'Bottomonium:qg2bbbar(3S1)[3S1(8)]q    = on,off,off',
#       'Bottomonium:qqbar2bbbar(3S1)[3S1(8)]g = on,off,off',
#       'Bottomonium:gg2bbbar(3S1)[1S0(8)]g    = on,off,off',
#       'Bottomonium:qg2bbbar(3S1)[1S0(8)]q    = on,off,off',
#       'Bottomonium:qqbar2bbbar(3S1)[1S0(8)]g = on,off,off',
#       'Bottomonium:gg2bbbar(3S1)[3PJ(8)]g    = on,off,off',
#       'Bottomonium:qg2bbbar(3S1)[3PJ(8)]q    = on,off,off',
#       'Bottomonium:qqbar2bbbar(3S1)[3PJ(8)]g = on,off,off',
#       '553:m0 = 15.0',
#       '553:mMin = 14.99995',
#       '553:mMax = 15.00005',
#       '553:onMode = off',
#       '553:onIfMatch = 15 -15', # doubt I need this stuff that I commented out, not making any upsilons
      '15:onMode = off',
      '15:onIfAll = 211 211 211',
      '15:onIfAll = 211 211 321',
      '15:onIfAll = 211 321 321',
      '15:onIfAll = 321 321 321',
      '15:onIfAll = 211 211 211',
      '15:onIfAll = 211 211 321',
      '15:onIfAll = 211 321 321',
      '15:onIfAll = 321 321 321',
      'TauDecays:mode = 3', #these TauDecays:etc lines come with help from Steve Mrena and other gen experts
      'TauDecays:tauPolarization = -1' #polarized sample, this will create a tau with polarization -1 (left handed, because particle) and will also make sure the anti tau has polarization + 1 (right handed, because anti particle)
      ), # close processParameters 

    parameterSets = cms.vstring(
      'pythia8CommonSettings',
      'pythia8CP5Settings',
      'pythia8PSweightsSettings',
      'processParameters')#close parameterSets
    ) #close PythiaParameters
) #close pythia8 gun



# Do Not need the filter for this
#upsilonfilter = cms.EDFilter("PythiaFilter",
#    Status = cms.untracked.int32(2),
#    MaxEta = cms.untracked.double(1000.0),
#    MinEta = cms.untracked.double(-1000.0),
#    MinPt = cms.untracked.double(2),
#    ParticleID = cms.untracked.int32(553)
#)

ProductionFilterSequence = cms.Sequence(generator) #*upsilonfilter
