#  _                _              ______  __  __ _      _                 _           
# | |              | |            |  ____|/ _|/ _(_)    (_)               (_)          
# | |     ___ _ __ | |_ ___  _ __ | |__  | |_| |_ _  ___ _  ___ _ __   ___ _  ___  ___ 
# | |    / _ \ '_ \| __/ _ \| '_ \|  __| |  _|  _| |/ __| |/ _ \ '_ \ / __| |/ _ \/ __|
# | |___|  __/ |_) | || (_) | | | | |____| | | | | | (__| |  __/ | | | (__| |  __/\__ \
# |______\___| .__/ \__\___/|_| |_|______|_| |_| |_|\___|_|\___|_| |_|\___|_|\___||___/
#            | |                                                                       
#            |_|                                                                       
#
# Description:  Using the TagProbeFitTreeAnalyzer to calculate
#               muon/electron reco, id and HLT efficiencies
# Date:         08/May/2014
# Author:       Janos Karancsi
# E-mail:       janos.karancsi@cern.ch

import FWCore.ParameterSet.Config as cms

import sys
import os
import ROOT

process = cms.Process("TagProbe")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

isMC = True
runOnTestFile = True
doEleGsf = True
doEleId  = True
doEleHLT = True
doMuTrk  = False
doMuId   = False
doMuHLT  = False

suffix = ""

if isMC:
    data_suffix = "_MC"
else:
    data_suffix = "_Data"
if runOnTestFile:
    suffix = "_Test"

outputEleGsf = "EffPlot"+data_suffix+"_Ele_Gsf"+suffix+".root"
outputEleId  = "EffPlot"+data_suffix+"_Ele_ID" +suffix+".root"
outputEleHLT = "EffPlot"+data_suffix+"_Ele_HLT"+suffix+".root"
outputMuTrk  = "EffPlot"+data_suffix+"_Mu_Trk"+suffix+".root"
outputMuId   = "EffPlot"+data_suffix+"_Mu_ID" +suffix+".root"
outputMuHLT  = "EffPlot"+data_suffix+"_Mu_HLT"+suffix+".root"

testFileEle = "TNP_Ntuple_MC_Ele.root"
testFileMu  = "TNP_Ntuple_MC_Mu.root"

#################################################################

#  _____          _                     _        
# |  __ \        | |                   | |       
# | |  | |  __ _ | |_  __ _  ___   ___ | |_  ___ 
# | |  | | / _` || __|/ _` |/ __| / _ \| __|/ __|
# | |__| || (_| || |_| (_| |\__ \|  __/| |_ \__ \
# |_____/  \__,_| \__|\__,_||___/ \___| \__||___/

mainDir  = "/data/gridout/jkarancs/SusyAnalysis/TagAndProbe_140512"

pathNames_Data_SingleEle = [
    mainDir+"/Ele_Data/SingleElectron_Run2012A-22Jan2013-v1/",
    mainDir+"/Ele_Data/SingleElectron_Run2012B-22Jan2013-v1/",
    mainDir+"/Ele_Data/SingleElectron_Run2012C-22Jan2013-v1/",
    mainDir+"/Ele_Data/SingleElectron_Run2012D-22Jan2013-v1/",
]
pathNames_Data_EleHad = [
    mainDir+"/Ele_Data/ElectronHad_Run2012A-22Jan2013-v1/",
    mainDir+"/Ele_Data/ElectronHad_Run2012B-22Jan2013-v1/",
    mainDir+"/Ele_Data/ElectronHad_Run2012C-22Jan2013-v1/",
    mainDir+"/Ele_Data/ElectronHad_Run2012D-22Jan2013-v1/",
]
pathNames_Data_SingleMu = [
    mainDir+"/Mu_Data/SingleMu_Run2012A-22Jan2013-v1/",
    mainDir+"/Mu_Data/SingleMu_Run2012B-22Jan2013-v1/",
    mainDir+"/Mu_Data/SingleMu_Run2012C-22Jan2013-v1/",
    mainDir+"/Mu_Data/SingleMu_Run2012D-22Jan2013-v1/",
]
pathNames_Data_MuHad = [
    mainDir+"/Mu_Data/MuHad_Run2012A-22Jan2013-v1/",
    mainDir+"/Mu_Data/MuHad_Run2012B-22Jan2013-v1/",
    mainDir+"/Mu_Data/MuHad_Run2012C-22Jan2013-v1/",
    mainDir+"/Mu_Data/MuHad_Run2012D-22Jan2013-v1/",
]
pathNames_MC_Ele = [
    mainDir+"/Ele_MC/TT_Summer12_DR53X-DynIneff-FlatPUmax50_START53_V7A-v1/",
    #mainDir+"/Ele_MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/",
    #mainDir+"/Ele_MC/QCD_HT-250To500_TuneZ2star_8TeV-madgraph-pythia6/",
    #mainDir+"/Ele_MC/QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia6/",
    #mainDir+"/Ele_MC/QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6/",
    #mainDir+"/Ele_MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/",
    #mainDir+"/Ele_MC/T_s-channel_TuneZ2star_8TeV-powheg-tauola/",
    #mainDir+"/Ele_MC/T_t-channel_TuneZ2star_8TeV-powheg-tauola/",
    #mainDir+"/Ele_MC/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/",
    #mainDir+"/Ele_MC/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola/",
    #mainDir+"/Ele_MC/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola/",
    #mainDir+"/Ele_MC/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/",
    #mainDir+"/Ele_MC/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/",
    #mainDir+"/Ele_MC/WW_TuneZ2star_8TeV_pythia6_tauola/",
    #mainDir+"/Ele_MC/WZ_TuneZ2star_8TeV_pythia6_tauola/",
    #mainDir+"/Ele_MC/ZZ_TuneZ2star_8TeV_pythia6_tauola/",
]
pathNames_MC_Mu = [
    mainDir+"/Mu_MC/TT_Summer12_DR53X-DynIneff-FlatPUmax50_START53_V7A-v1/",
    #mainDir+"/Mu_MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/",
    #mainDir+"/Mu_MC/QCD_HT-250To500_TuneZ2star_8TeV-madgraph-pythia6/",
    #mainDir+"/Mu_MC/QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia6/",
    #mainDir+"/Mu_MC/QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6/",
    #mainDir+"/Mu_MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/",
    #mainDir+"/Mu_MC/T_s-channel_TuneZ2star_8TeV-powheg-tauola/",
    #mainDir+"/Mu_MC/T_t-channel_TuneZ2star_8TeV-powheg-tauola/",
    #mainDir+"/Mu_MC/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/",
    #mainDir+"/Mu_MC/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola/",
    #mainDir+"/Mu_MC/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola/",
    #mainDir+"/Mu_MC/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/",
    #mainDir+"/Mu_MC/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/",
    #mainDir+"/Mu_MC/WW_TuneZ2star_8TeV_pythia6_tauola/",
    #mainDir+"/Mu_MC/WZ_TuneZ2star_8TeV_pythia6_tauola/",
    #mainDir+"/Mu_MC/ZZ_TuneZ2star_8TeV_pythia6_tauola/",
]

if isMC:
    weightVar = "weight"
    theUnbinned = cms.vstring("mass",weightVar)
    pathNames_Ele_All = [ pathNames_MC_Ele ]
    pathNames_Mu_All = [ pathNames_MC_Mu ]
else:
    weightVar = ""
    theUnbinned = cms.vstring("mass")
    pathNames_Ele_All = [ pathNames_Data_SingleEle, pathNames_Data_EleHad ]
    pathNames_Mu_All = [ pathNames_Data_SingleMu, pathNames_Data_MuHad ]

selectedEleAll = cms.vstring()
selectedEleSingle = cms.vstring()
selectedMuAll = cms.vstring()
selectedMuSingle = cms.vstring()
if runOnTestFile:
    selectedEleAll.append(testFileEle)
    selectedEleSingle.append(testFileEle)
    selectedMuAll.append(testFileMu)
    selectedMuSingle.append(testFileMu)
else:
    if doEleGsf or doEleId or doEleHLT:
        for (i,dir) in enumerate(pathNames_Ele_All):
            if i == 0 or doEleHLT:
                for (j,dir2) in enumerate(dir):
                    print "Adding files in " + dir2 + ":"
                    filenames = os.listdir( dir2 )
                    # Select files to be added
                    selectedFiles = []
                    addedFiles = []
                    for file in filenames:
                        # discard non-root files
                        if (file.rpartition(".")[2] != "root"): continue
                        # discard duplicates
                        fname = file
                        #fname = file.rpartition("_")[0].rpartition("_")[0] # crab format
                        if (addedFiles.count(fname) > 0): 
                            print "   omitting duplicate file: " + file
                            continue
                        else:
                            addedFiles.append( fname )
                            selectedFiles.append( dir2 + file )
                    for file in selectedFiles:
                        selectedEleAll.append(file)
                        if i == 0:
                            selectedEleSingle.append(file)
                        print "      " + str(len(selectedFiles)) + " selected files"
        if doEleHLT and not isMC:
            print "SingleEle (for Gsf and Id): " + str(len(selectedEleSingle)) + " files added"
            print "SingleEle + EleHad (for HLT):" + str(len(selectedEleAll)) + " files added"
        else:
            print "All Electrons: " + str(len(selectedEleAll)) + " files added"
    
    if doMuTrk or doMuId or doMuHLT:
        for (i,dir) in enumerate(pathNames_Mu_All):
            if i == 0 or doMuHLT:
                for (j,dir2) in enumerate(dir):
                    print "Adding files in " + dir2 + ":"
                    filenames = os.listdir( dir2 )
                    # Select files to be added
                    selectedFiles = []
                    addedFiles = []
                    for file in filenames:
                        # discard non-root files
                        if (file.rpartition(".")[2] != "root"): continue
                        # discard duplicates
                        fname = file
                        #fname = file.rpartition("_")[0].rpartition("_")[0] # for crab format
                        if (addedFiles.count(fname) > 0): 
                            print "   omitting duplicate file: " + file
                            continue
                        else:
                            addedFiles.append( fname )
                            selectedFiles.append( dir2 + file )
                    for file in selectedFiles:
                        selectedMuAll.append(file)
                        if i == 0:
                            selectedMuSingle.append(file)
                        print "      " + str(len(selectedFiles)) + " selected files"
        if doMuHLT and not isMC:
            print "SingleMu (for Trk and Id): " + str(len(selectedMuSingle)) + " files added"
            print "SingleMu + MuHad (for HLT):" + str(len(selectedMuAll)) + " files added"
        else:
            print "All Muons: " + str(len(selectedMuAll)) + " files added"



#################################################################

#  _______                                         _                             
# |__   __|                   /\                  | |                            
#    | | _ __  ___   ___     /  \    _ __    __ _ | | _   _  ____ ___  _ __  ___ 
#    | || '__|/ _ \ / _ \   / /\ \  | '_ \  / _` || || | | ||_  // _ \| '__|/ __|
#    | || |  |  __/|  __/  / ____ \ | | | || (_| || || |_| | / /|  __/| |   \__ \
#    |_||_|   \___| \___| /_/    \_\|_| |_| \__,_||_| \__, |/___|\___||_|   |___/
#                                                      __/ |                     
#                                                     |___/                      

if isMC:
    weightVar = "weight"
    theUnbinned = cms.vstring("mass",weightVar)
else:
    weightVar = ""
    theUnbinned = cms.vstring("mass")

process.TPFTA_template = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    InputTreeName = cms.string("fitter_tree"),
    NumCPU = cms.uint32(6),
    SaveWorkspace = cms.bool(False),
    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(40),
    WeightVariable = cms.string(weightVar),
    PDFs = cms.PSet(
        vpvPlusExpo = cms.vstring(
            "Voigtian::signal1p(mass, mean1[90,80,100], width[2.495], sigma1[2,1,3])",
            "Voigtian::signal2p(mass, mean2[90,80,100], width,        sigma2[4,2,10])",
            "SUM::signalPass(vFrac[0.8,0,1]*signal1p, signal2p)",
            "Voigtian::signal1f(mass, mean1[90,80,100], width[2.495], sigma1[2,1,3])",
            "Voigtian::signal2f(mass, mean2[90,80,100], width,        sigma2[4,2,10])",
            "SUM::signalFail(vFrac[0.8,0,1]*signal1f, signal2f)",
            "Exponential::backgroundPass(mass, lp[-0.1,-1,0.1])",
            "Exponential::backgroundFail(mass, lf[-0.1,-1,0.1])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.90]"
        ),
    ),
    Cuts = cms.PSet(),
)

#  ______  _              _                           
# |  ____|| |            | |                          
# | |__   | |  ___   ___ | |_  _ __  ___   _ __   ___ 
# |  __|  | | / _ \ / __|| __|| '__|/ _ \ | '_ \ / __|
# | |____ | ||  __/| (__ | |_ | |  | (_) || | | |\__ \
# |______||_| \___| \___| \__||_|   \___/ |_| |_||___/
#

# Electrons - Sc -> Gsf
process.FitEleGsfEff = process.TPFTA_template.clone(
    InputFileNames = selectedEleSingle,
    InputDirectoryName = cms.string("SuperClusterToGsfElectronPATTag"),
    OutputFileName = cms.string(outputEleGsf),
        
    Variables = cms.PSet(
        mass     = cms.vstring("Tag-Probe Mass", "40", "140", "GeV/c^{2}"),
        pt       = cms.vstring("SuperCluster p_{T}", "0", "1000", "GeV"), # exactly same as et
        et       = cms.vstring("SuperCluster E_{T}", "0", "1000", "GeV"),
        phi      = cms.vstring("SuperCluster #phi at Vertex", "-3.1416", "3.1416", ""),
        eta      = cms.vstring("SuperCluster #eta", "-2.5", "2.5", ""),
        abs_eta  = cms.vstring("SuperCluster |#eta|", "0", "2.5", ""),
        nvtx     = cms.vstring("Number of Vertices", "0", "999", ""),
        pileup   = cms.vstring("Number of Inelastic Collisions", "0", "60", ""),
        instlumi = cms.vstring("Instantaneous Luminosity", "0", "99999", "#mub^{-1}s^{-1}"),
        weight   = cms.vstring("MC Event Weight", "0", "100000000", ""),
        drjet    = cms.vstring("SuperCluster-Jet_{pt>40} #DeltaR", "0", "10000000000", ""),
        njet     = cms.vstring("Number of Jets_{pt>40}", "0", "100", ""),
        ht       = cms.vstring("Event HT", "0", "600", "GeV/c"),
        met      = cms.vstring("Event #slash{E}_{T}", "0", "600", "GeV"),
    ),
    
    Categories = cms.PSet(
        #mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        passing_gsf = cms.vstring("SuperCluster Passing as Gsf Electron", "dummy[pass=1,fail=0]"),
        passing_pat = cms.vstring("SuperCluster Passing as Pat Electron (with cuts)", "dummy[pass=1,fail=0]"),
    ),
)

# Electrons - Gsf -> Iso,Id
process.FitEleIdEff = process.TPFTA_template.clone(
    InputFileNames = selectedEleSingle,
    InputDirectoryName = cms.string("GsfElectronToIdPATTag"),
    OutputFileName = cms.string(outputEleId),
    
    Variables = cms.PSet(
        mass        = cms.vstring("Tag-Probe Mass", "40", "140", "GeV/c^{2}"),
        pt          = cms.vstring("Gsf Electron p_{T}", "0", "1000", "GeV/c"),
        et          = cms.vstring("Gsf Electron E_{T}", "0", "1000", "GeV"),
        phi         = cms.vstring("Gsf Electron #phi at Vertex", "-3.1416", "3.1416", ""),
        eta         = cms.vstring("Gsf Electron #eta", "-2.5", "2.5", ""),
        abs_eta     = cms.vstring("Gsf Electron |#eta|", "0", "2.5", ""),
        #sc_et       = cms.vstring("SuperCluster E_{T}", "0", "1000", "GeV"),
        #sc_phi      = cms.vstring("SuperCluster #phi at Vertex", "-3.1416", "3.1416", ""),
        #sc_eta      = cms.vstring("SuperCluster #eta", "-2.5", "2.5", ""),
        #track_pt    = cms.vstring("Gsf Track p_{T}", "0", "1000", "GeV/c"),
        #track_phi   = cms.vstring("Gsf Track #phi at Vertex", "-3.1416", "3.1416", ""),
        #track_eta   = cms.vstring("Gsf Track #eta", "-2.5", "2.5", ""),        
        nvtx        = cms.vstring("Number of Vertices", "0", "999", ""),
        pileup      = cms.vstring("Number of Inelastic Collisions", "0", "60", ""),
        instlumi    = cms.vstring("Instantaneous Luminosity", "0", "99999", "#mub^{-1}s^{-1}"),
        weight      = cms.vstring("MC Event Weight", "0", "100000000", ""),
        #d0_b        = cms.vstring("Gsf Electron d0_{Beamspot}"  , "-20", "20", "cm"),
        #dz_b        = cms.vstring("Gsf Electron dz_{Beamspot}"  , "-20", "20", "cm"),
        d0_v        = cms.vstring("Gsf Electron d0_{Vertex}"  , "-20", "20", "cm"),
        dz_v        = cms.vstring("Gsf Electron dz_{Vertex}"  , "-20", "20", "cm"),
        drjet       = cms.vstring("Gsf Electron-Jet_{pt>40} #DeltaR", "0", "10000000000", ""),
        njet        = cms.vstring("Number of Jets_{pt>40}", "0", "100", ""),
        ht          = cms.vstring("Event H_{T}", "0", "600", "GeV/c"),
        met         = cms.vstring("Event #slash{E}_{T}", "0", "600", "GeV"),
        st          = cms.vstring("Gsf Electron S_{T}", "0", "600", "GeV"),
        absdeltapt  = cms.vstring("Gsf Electron |p_{T,Reco] - p_{T,PF}|", "-9999", "9999", "GeV/c"), # -9999: under 10 GeV
        reliso      = cms.vstring("Gsf Electron Relative Isolation", "0", "9999999", ""),
        #passConvRej = cms.vstring("Gsf Electron Passing Conversion Rejection", "0", "1", ""), # passing_pat already has it
    ),
    
    Categories = cms.PSet(
        #mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        passing_pat = cms.vstring("Gsf Electron Passing as Pat Electron (with cuts)", "dummy[pass=1,fail=0]"),
    ),
    
    Cuts = cms.PSet(
        d0_v_m002         = cms.vstring("d0 > -0.02",                "d0_v",       "-0.02"),
        d0_v_p002         = cms.vstring("d0 < 0.02",                 "d0_v",       "0.02"),
        dz_v_m01          = cms.vstring("dz > -0.1",                 "dz_v",       "-0.1"),
        dz_v_p01          = cms.vstring("dz < 0.1",                  "dz_v",       "0.1"),
        drjet_03          = cms.vstring("deltaR (nearest Jet > 0.3", "drjet",      "0.3"),
        reliso_015        = cms.vstring("RelIso  < 0.15",            "reliso",     "0.15"),
        absdeltapt_p10    = cms.vstring("|delta pt| < 10",           "absdeltapt", "10.",),
        absdeltapt_m10000 = cms.vstring("|delta pt| > -10000",       "absdeltapt", "-10000.",), # -9999 for pt<10 GeV
    ),
    
)

# Electrons - Id -> HLT
process.FitEleHLTEff = process.TPFTA_template.clone(
    InputFileNames = selectedEleAll,
    InputDirectoryName = cms.string("goodPATEleToHLT"),
    OutputFileName = cms.string(outputEleHLT),
    
    Variables = cms.PSet(
        mass       = cms.vstring("Tag-Probe Mass", "40", "140", "GeV/c^{2}"),
        pt         = cms.vstring("Pat Electron p_{T}", "0", "1000", "GeV/c"),
        et         = cms.vstring("Pat Electron E_{T}", "0", "1000", "GeV"),
        phi        = cms.vstring("Pat Electron #phi at Vertex", "-3.1416", "3.1416", ""),
        eta        = cms.vstring("Pat Electron #eta", "-2.5", "2.5", ""),
        abs_eta    = cms.vstring("Pat Electron |#eta|", "0", "2.5", ""),
        #sc_et      = cms.vstring("SuperCluster E_{T}", "0", "1000", "GeV"),
        #sc_phi     = cms.vstring("SuperCluster #phi at Vertex", "-3.1416", "3.1416", ""),
        #sc_eta     = cms.vstring("SuperCluster #eta", "-2.5", "2.5", ""),
        #track_pt   = cms.vstring("Gsf Track p_{T}", "0", "1000", "GeV/c"),
        #track_phi  = cms.vstring("Gsf Track #phi at Vertex", "-3.1416", "3.1416", ""),
        #track_eta  = cms.vstring("Gsf Track #eta", "-2.5", "2.5", ""),        
        nvtx       = cms.vstring("Number of Vertices", "0", "999", ""),
        pileup     = cms.vstring("Number of Inelastic Collisions", "0", "60", ""),
        instlumi   = cms.vstring("Instantaneous Luminosity", "0", "99999", "#mub^{-1}s^{-1}"),
        weight     = cms.vstring("MC Event Weight", "0", "100000000", ""),
        d0_b       = cms.vstring("Pat Electron d0_{Beamspot}", "-20", "20", "cm"),
        dz_b       = cms.vstring("Pat Electron dz_{Beamspot}", "-20", "20", "cm"),
        d0_v       = cms.vstring("Pat Electron d0_{Vertex}"  , "-20", "20", "cm"),
        dz_v       = cms.vstring("Pat Electron dz_{Vertex}"  , "-20", "20", "cm"),
        drjet      = cms.vstring("Pat Electron-Jet_{pt>40} #DeltaR", "0", "10000000000", ""),
        njet       = cms.vstring("Number of Jets_{pt>40}", "0", "100", ""),
        ht         = cms.vstring("Event H_{T}", "0", "600", "GeV/c"),
        met        = cms.vstring("Event #slash{E}_{T}", "0", "600", "GeV"),
        st         = cms.vstring("Pat Electron S_{T}", "0", "600", "GeV"),
        absdeltapt = cms.vstring("Pat Electron |p_{T,Reco] - p_{T,PF}|", "-9999", "9999", "GeV/c"), # -9999: under 10 GeV
        reliso     = cms.vstring("Pat Electron Relative Isolation", "0", "9999999", ""),
        # 2012 Cross Triggers
        event_passing_HLT_Ele8_Jet30                     = cms.vstring("Event passing HLT_Ele8_CaloIdT_TrkIdVL_Jet30_v*",                                      "0", "2", ""),
        event_passing_HLT_CleanPFHT350_Ele5_PFMET45      = cms.vstring("Event passing HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*",     "0", "2", ""),
        event_passing_HLT_CleanPFHT350_Ele5_PFMET50      = cms.vstring("Event passing HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*",     "0", "2", ""),
        event_passing_HLT_CleanPFHT300_Ele15_PFMET45     = cms.vstring("Event passing HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*",    "0", "2", ""),
        event_passing_HLT_CleanPFHT300_Ele15_PFMET50     = cms.vstring("Event passing HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*",    "0", "2", ""),
        event_passing_HLT_CleanPFHT300_Ele40             = cms.vstring("Event passing HLT_CleanPFHT300_Ele40_CaloIdVT_TrkIdT_v*",                              "0", "2", ""),
        event_passing_HLT_CleanPFHT300_Ele60             = cms.vstring("Event passing HLT_CleanPFHT300_Ele60_CaloIdVT_TrkIdT_v*",                              "0", "2", ""),
        event_passing_HLT_CleanPFNoPUHT350_Ele5_PFMET45  = cms.vstring("Event passing HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*", "0", "2", ""),
        event_passing_HLT_CleanPFNoPUHT350_Ele5_PFMET50  = cms.vstring("Event passing HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*", "0", "2", ""),
        event_passing_HLT_CleanPFNoPUHT300_Ele15_PFMET45 = cms.vstring("Event passing HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*","0", "2", ""),
        event_passing_HLT_CleanPFNoPUHT300_Ele15_PFMET50 = cms.vstring("Event passing HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*","0", "2", ""),
        event_passing_HLT_CleanPFNoPUHT300_Ele40         = cms.vstring("Event passing HLT_CleanPFNoPUHT300_Ele40_CaloIdVT_TrkIdT_v*",                          "0", "2", ""),
        event_passing_HLT_CleanPFNoPUHT300_Ele60         = cms.vstring("Event passing HLT_CleanPFNoPUHT300_Ele60_CaloIdVT_TrkIdT_v*",                          "0", "2", ""),
    ),
    
    Categories = cms.PSet(
        #mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        passing_HLT_Ele22_CaloIdL_CaloIsoVL                  = cms.vstring("Probe passing HLT_Ele22_CaloIdL_CaloIsoVL_v*",                                        "dummy[pass=1,fail=0]"),
        passing_HLT_Ele27_CaloIdL_CaloIsoVL                  = cms.vstring("Probe passing HLT_Ele27_CaloIdL_CaloIsoVL_v*",                                        "dummy[pass=1,fail=0]"),
        passing_HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL = cms.vstring("Probe passing HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",                       "dummy[pass=1,fail=0]"),
        passing_HLT_Ele27_WP80                               = cms.vstring("Probe passing HLT_Ele27_WP80_v*",                                                     "dummy[pass=1,fail=0]"),
        passing_HLT_Ele30_CaloIdVT_TrkIdT                    = cms.vstring("Probe passing HLT_Ele30_CaloIdVT_TrkIdT_v*",                                          "dummy[pass=1,fail=0]"),
        passing_HLT_Ele32_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL = cms.vstring("Probe passing HLT_Ele32_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",                       "dummy[pass=1,fail=0]"),
        passing_HLT_Ele65_CaloIdVT_TrkIdT                    = cms.vstring("Probe passing HLT_Ele65_CaloIdVT_TrkIdT_v*",                                          "dummy[pass=1,fail=0]"),
        passing_HLT_Ele80_CaloIdVT_TrkIdT                    = cms.vstring("Probe passing HLT_Ele80_CaloIdVT_TrkIdT_v*",                                          "dummy[pass=1,fail=0]"),
        passing_HLT_Ele100_CaloIdVT_TrkIdT                   = cms.vstring("Probe passing HLT_Ele100_CaloIdVT_TrkIdT_v*",                                         "dummy[pass=1,fail=0]"),
        passing_HLT_Ele80_CaloIdVT_GsfTrkIdT                 = cms.vstring("Probe passing HLT_Ele80_CaloIdVT_GsfTrkIdT_v*",                                       "dummy[pass=1,fail=0]"),
        passing_HLT_Ele90_CaloIdVT_GsfTrkIdT                 = cms.vstring("Probe passing HLT_Ele90_CaloIdVT_GsfTrkIdT_v*",                                       "dummy[pass=1,fail=0]"),
        passing_HLT_Ele8_CaloIdT_TrkIdVL_Jet30               = cms.vstring("Probe passing HLT_Ele8_CaloIdT_TrkIdVL_Jet30_v*",                                     "dummy[pass=1,fail=0]"),
        passing_HLT_CleanPFHT350_Ele5_PFMET45                = cms.vstring("Probe passing HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*",    "dummy[pass=1,fail=0]"),
        passing_HLT_CleanPFHT350_Ele5_PFMET50                = cms.vstring("Probe passing HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*",    "dummy[pass=1,fail=0]"),
        passing_HLT_CleanPFNoPUHT350_Ele5_PFMET45            = cms.vstring("Probe passing HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*","dummy[pass=1,fail=0]"),
        passing_HLT_CleanPFNoPUHT350_Ele5_PFMET50            = cms.vstring("Probe passing HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*","dummy[pass=1,fail=0]"),
        passing_HLT_CleanPFHT300_Ele15_PFMET45               = cms.vstring("Probe passing HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*",   "dummy[pass=1,fail=0]"),
        passing_HLT_CleanPFHT300_Ele15_PFMET50               = cms.vstring("Probe passing HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*",   "dummy[pass=1,fail=0]"),
        passing_HLT_CleanPFNoPUHT300_Ele15_PFMET45           = cms.vstring("Probe passing HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET4_v*","dummy[pass=1,fail=0]"),
        passing_HLT_CleanPFNoPUHT300_Ele15_PFMET50           = cms.vstring("Probe passing HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET5_v*","dummy[pass=1,fail=0]"),
        passing_HLT_CleanPFHT300_Ele40                       = cms.vstring("Probe passing HLT_CleanPFHT300_Ele40_CaloIdVT_TrkIdT_v*",                             "dummy[pass=1,fail=0]"),
        passing_HLT_CleanPFHT300_Ele60                       = cms.vstring("Probe passing HLT_CleanPFHT300_Ele60_CaloIdVT_TrkIdT_v*",                             "dummy[pass=1,fail=0]"),
        passing_HLT_CleanPFNoPUHT300_Ele40                   = cms.vstring("Probe passing HLT_CleanPFNoPUHT300_Ele40_CaloIdVT_TrkIdT_v*",                         "dummy[pass=1,fail=0]"),
        passing_HLT_CleanPFNoPUHT300_Ele60                   = cms.vstring("Probe passing HLT_CleanPFNoPUHT300_Ele60_CaloIdVT_TrkIdT_v*",                         "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Ele22_CaloIdL_CaloIsoVL                  = cms.vstring("Tag passing HLT_Ele22_CaloIdL_CaloIsoVL_v*",                                        "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Ele27_CaloIdL_CaloIsoVL                  = cms.vstring("Tag passing HLT_Ele27_CaloIdL_CaloIsoVL_v*",                                        "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL = cms.vstring("Tag passing HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",                       "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Ele27_WP80                               = cms.vstring("Tag passing HLT_Ele27_WP80_v*",                                                     "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Ele30_CaloIdVT_TrkIdT                    = cms.vstring("Tag passing HLT_Ele30_CaloIdVT_TrkIdT_v*",                                          "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Ele32_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL = cms.vstring("Tag passing HLT_Ele32_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",                       "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Ele65_CaloIdVT_TrkIdT                    = cms.vstring("Tag passing HLT_Ele65_CaloIdVT_TrkIdT_v*",                                          "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Ele80_CaloIdVT_TrkIdT                    = cms.vstring("Tag passing HLT_Ele80_CaloIdVT_TrkIdT_v*",                                          "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Ele100_CaloIdVT_TrkIdT                   = cms.vstring("Tag passing HLT_Ele100_CaloIdVT_TrkIdT_v*",                                         "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Ele80_CaloIdVT_GsfTrkIdT                 = cms.vstring("Tag passing HLT_Ele80_CaloIdVT_GsfTrkIdT_v*",                                       "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Ele90_CaloIdVT_GsfTrkIdT                 = cms.vstring("Tag passing HLT_Ele90_CaloIdVT_GsfTrkIdT_v*",                                       "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Ele8_CaloIdT_TrkIdVL_Jet30               = cms.vstring("Tag passing HLT_Ele8_CaloIdT_TrkIdVL_Jet30_v*",                                     "dummy[pass=1,fail=0]"),
        tag_passing_HLT_CleanPFHT350_Ele5_PFMET45                = cms.vstring("Tag passing HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*",    "dummy[pass=1,fail=0]"),
        tag_passing_HLT_CleanPFHT350_Ele5_PFMET50                = cms.vstring("Tag passing HLT_CleanPFHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*",    "dummy[pass=1,fail=0]"),
        tag_passing_HLT_CleanPFNoPUHT350_Ele5_PFMET45            = cms.vstring("Tag passing HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*","dummy[pass=1,fail=0]"),
        tag_passing_HLT_CleanPFNoPUHT350_Ele5_PFMET50            = cms.vstring("Tag passing HLT_CleanPFNoPUHT350_Ele5_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*","dummy[pass=1,fail=0]"),
        tag_passing_HLT_CleanPFHT300_Ele15_PFMET45               = cms.vstring("Tag passing HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET45_v*",   "dummy[pass=1,fail=0]"),
        tag_passing_HLT_CleanPFHT300_Ele15_PFMET50               = cms.vstring("Tag passing HLT_CleanPFHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET50_v*",   "dummy[pass=1,fail=0]"),
        tag_passing_HLT_CleanPFNoPUHT300_Ele15_PFMET45           = cms.vstring("Tag passing HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET4_v*","dummy[pass=1,fail=0]"),
        tag_passing_HLT_CleanPFNoPUHT300_Ele15_PFMET50           = cms.vstring("Tag passing HLT_CleanPFNoPUHT300_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_PFMET5_v*","dummy[pass=1,fail=0]"),
        tag_passing_HLT_CleanPFHT300_Ele40                       = cms.vstring("Tag passing HLT_CleanPFHT300_Ele40_CaloIdVT_TrkIdT_v*",                             "dummy[pass=1,fail=0]"),
        tag_passing_HLT_CleanPFHT300_Ele60                       = cms.vstring("Tag passing HLT_CleanPFHT300_Ele60_CaloIdVT_TrkIdT_v*",                             "dummy[pass=1,fail=0]"),
        tag_passing_HLT_CleanPFNoPUHT300_Ele40                   = cms.vstring("Tag passing HLT_CleanPFNoPUHT300_Ele40_CaloIdVT_TrkIdT_v*",                         "dummy[pass=1,fail=0]"),
        tag_passing_HLT_CleanPFNoPUHT300_Ele60                   = cms.vstring("Tag passing HLT_CleanPFNoPUHT300_Ele60_CaloIdVT_TrkIdT_v*",                         "dummy[pass=1,fail=0]"),
    ),
)

#  __  __  _    _   ____   _   _   _____ 
# |  \/  || |  | | / __ \ | \ | | / ____|
# | \  / || |  | || |  | ||  \| || (___  
# | |\/| || |  | || |  | || . ` | \___ \ 
# | |  | || |__| || |__| || |\  | ____) |
# |_|  |_| \____/  \____/ |_| \_||_____/ 

# Muons - Sta -> Trk
process.FitMuTrkEff = process.TPFTA_template.clone(
    InputFileNames = selectedMuSingle,
    InputDirectoryName = cms.string("fitTkFromSta"),
    OutputFileName = cms.string(outputMuTrk),
    
    Categories = cms.PSet(
        #mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        passing_trk = cms.vstring("Sta Muon Passing as Track", "dummy[pass=1,fail=0]"),
    ),
    
    Variables = cms.PSet(
        mass        = cms.vstring("Tag-Probe Mass", "40", "140", "GeV/c^{2}"),
        pt          = cms.vstring("Sta Muon p_{T}", "0", "9999", "GeV/c"),
        phi         = cms.vstring("Sta Muon #phi at Vertex", "-3.1416", "3.1416", ""),
        eta         = cms.vstring("Sta Muon #eta", "-2.4", "2.4", ""),
        abs_eta     = cms.vstring("Sta Muon |#eta|", "0", "2.4", ""),
        nvtx        = cms.vstring("Number of Vertices", "0", "9999", ""),
        pileup      = cms.vstring("Number of Inelastic Collisions", "0", "9999", ""),
        instlumi    = cms.vstring("Instantaneous Luminosity", "0", "99999", "#mub^{-1}s^{-1}"),
        weight      = cms.vstring("MC Event Weight", "0", "9999", ""),
        #d0_b        = cms.vstring("Sta Muon d0_{Beamspot}"  , "-20", "20", "cm"),
        #dz_b        = cms.vstring("Sta Muon dz_{Beamspot}"  , "-20", "20", "cm"),
        d0_v        = cms.vstring("Sta Muon d0_{Vertex}"  , "-2", "2", "cm"),
        dz_v        = cms.vstring("Sta Muon dz_{Vertex}"  , "-25", "25", "cm"),
        drjet       = cms.vstring("Sta Muon-Jet_{pt>40} #DeltaR", "0", "9999", ""),
        njet        = cms.vstring("Number of Jets_{pt>40}", "0", "9999", ""),
        ht          = cms.vstring("Event H_{T}", "0", "9999", "GeV/c"),
        met         = cms.vstring("Event #slash{E}_{T}", "0", "9999", "GeV"),
        st          = cms.vstring("Sta Muon S_{T}", "0", "9999", "GeV"),
    ),
)

# Muons - Trk -> Iso,Id
process.FitMuIdEff = process.TPFTA_template.clone(
    InputFileNames = selectedMuSingle,
    InputDirectoryName = cms.string("fitGlbFromTk"),
    OutputFileName = cms.string(outputMuId),
    
    Variables = cms.PSet(
        mass        = cms.vstring("Tag-Probe Mass", "40", "140", "GeV/c^{2}"),
        pt          = cms.vstring("Tracker Muon p_{T}", "0", "9999", "GeV/c"),
        phi         = cms.vstring("Tracker Muon #phi at Vertex", "-3.1416", "3.1416", ""),
        eta         = cms.vstring("Tracker Muon #eta", "-2.4", "2.4", ""),
        abs_eta     = cms.vstring("Tracker Muon |#eta|", "0", "2.4", ""),
        nvtx        = cms.vstring("Number of Vertices", "0", "9999", ""),
        pileup      = cms.vstring("Number of Inelastic Collisions", "0", "9999", ""),
        instlumi    = cms.vstring("Instantaneous Luminosity", "0", "99999", "#mub^{-1}s^{-1}"),
        weight      = cms.vstring("MC Event Weight", "0", "9999", ""),
        #d0_b        = cms.vstring("Tracker Muon d0_{Beamspot}"  , "-20", "20", "cm"),
        #dz_b        = cms.vstring("Tracker Muon dz_{Beamspot}"  , "-20", "20", "cm"),
        d0_v        = cms.vstring("Tracker Muon d0_{Vertex}"  , "-2", "2", "cm"),
        dz_v        = cms.vstring("Tracker Muon dz_{Vertex}"  , "-25", "25", "cm"),
        drjet       = cms.vstring("Tracker Muon-Jet_{pt>40} #DeltaR", "0", "9999", ""),
        njet        = cms.vstring("Number of Jets_{pt>40}", "0", "9999", ""),
        ht          = cms.vstring("Event H_{T}", "0", "9999", "GeV/c"),
        met         = cms.vstring("Event #slash{E}_{T}", "0", "9999", "GeV"),
        st          = cms.vstring("Tracker Muon S_{T}", "0", "9999", "GeV"),
        absdeltapt  = cms.vstring("Global (Pat) Muon |p_{T,Reco] - p_{T,PF}|", "0", "999999", "GeV/c"), # -9999: under 10 GeV, 999999: no global match
        deltapt     = cms.vstring("Global (Pat) Muon |(p_{T,Reco] - p_{T,PF})/p_{T,Reco]|", "0", "999999", "GeV/c"), # -9999: under 10 GeV, 999999: no global match
        reliso      = cms.vstring("Global (Pat) Muon Relative Isolation (PF)", "0", "9999", ""), # 999999: no global match
    ),
    
    Categories = cms.PSet(
        passing_glb = cms.vstring("Tracker Muon Passing as Global (Pat) Muon (With cuts)", "dummy[pass=1,fail=0]"),
    ),
    
    Cuts = cms.PSet(
        d0_v_m002         = cms.vstring("d0 > -0.02",                "d0_v",       "-0.02"),
        d0_v_p002         = cms.vstring("d0 < 0.02",                 "d0_v",       "0.02"),
        dz_v_m05          = cms.vstring("dz > -0.5",                 "dz_v",       "-0.5"),
        dz_v_p05          = cms.vstring("dz < 0.5",                  "dz_v",       "0.5"),
        drjet_03          = cms.vstring("deltaR (nearest Jet > 0.3", "drjet",      "0.3"),
        reliso_012        = cms.vstring("ReliIso (PF) < 0.12",       "reliso",     "0.12"),
        absdeltapt_p5     = cms.vstring("|delta pt| < 5",            "absdeltapt", "5.",),
        absdeltapt_m10000 = cms.vstring("|delta pt| > -10000",       "absdeltapt", "-10000.",),
    ),
)

# Muons - Id -> HLT
process.FitMuHLTEff = process.TPFTA_template.clone(
    InputFileNames = selectedMuAll,
    InputDirectoryName = cms.string("fitHltFromGlb"),
    OutputFileName = cms.string(outputMuHLT),
    
    Variables = cms.PSet(
        mass        = cms.vstring("Tag-Probe Mass", "40", "140", "GeV/c^{2}"),
        pt          = cms.vstring("Global (Pat) Muon p_{T}", "0", "9999", "GeV/c"),
        phi         = cms.vstring("Global (Pat) Muon #phi at Vertex", "-3.1416", "3.1416", ""),
        eta         = cms.vstring("Global (Pat) Muon #eta", "-2.4", "2.4", ""),
        abs_eta     = cms.vstring("Global (Pat) Muon |#eta|", "0", "2.4", ""),
        nvtx        = cms.vstring("Number of Vertices", "0", "9999", ""),
        pileup      = cms.vstring("Number of Inelastic Collisions", "0", "9999", ""),
        instlumi    = cms.vstring("Instantaneous Luminosity", "0", "99999", "#mub^{-1}s^{-1}"),
        weight      = cms.vstring("MC Event Weight", "0", "9999", ""),
        #d0_b        = cms.vstring("Global (Pat) Muon d0_{Beamspot}"  , "-20", "20", "cm"),
        #dz_b        = cms.vstring("Global (Pat) Muon dz_{Beamspot}"  , "-20", "20", "cm"),
        d0_v        = cms.vstring("Global (Pat) Muon d0_{Vertex}"  , "-2", "2", "cm"),
        dz_v        = cms.vstring("Global (Pat) Muon dz_{Vertex}"  , "-25", "25", "cm"),
        drjet       = cms.vstring("Global (Pat) Muon-Jet_{pt>40} #DeltaR", "0", "9999", ""),
        njet        = cms.vstring("Number of Jets_{pt>40}", "0", "9999", ""),
        ht          = cms.vstring("Event H_{T}", "0", "9999", "GeV/c"),
        met         = cms.vstring("Event #slash{E}_{T}", "0", "9999", "GeV"),
        st          = cms.vstring("Global (Pat) Muon S_{T}", "0", "9999", "GeV"),
        absdeltapt  = cms.vstring("Global (Pat) Muon |p_{T,Reco] - p_{T,PF}|", "0", "999999", "GeV/c"), # -9999: under 10 GeV, 999999: no global match
        reliso      = cms.vstring("Global (Pat) Muon Relative Isolation (PF)", "0", "9999", ""),
        event_passing_HLT_Mu40_HT200               = cms.vstring("Event Passing HLT_Mu40_HT200",              "0", "2", ""),
        event_passing_HLT_Mu40_FJHT200             = cms.vstring("Event Passing HLT_Mu40_FJHT200",            "0", "2", ""),
        event_passing_HLT_Mu40_PFHT350             = cms.vstring("Event Passing HLT_Mu40_PFHT350",            "0", "2", ""),
        event_passing_HLT_Mu60_PFHT350             = cms.vstring("Event Passing HLT_Mu60_PFHT350",            "0", "2", ""),
        event_passing_HLT_Mu40_PFNoPUHT350         = cms.vstring("Event Passing HLT_Mu40_PFNoPUHT350",        "0", "2", ""),
        event_passing_HLT_Mu60_PFNoPUHT350         = cms.vstring("Event Passing HLT_Mu60_PFNoPUHT350",        "0", "2", ""),
        event_passing_HLT_PFHT350_Mu15_PFMET45     = cms.vstring("Event Passing HLT_PFHT350_Mu15_PFMET45",    "0", "2", ""),
        event_passing_HLT_PFHT350_Mu15_PFMET50     = cms.vstring("Event Passing HLT_PFHT350_Mu15_PFMET50",    "0", "2", ""),
        event_passing_HLT_PFHT400_Mu5_PFMET45      = cms.vstring("Event Passing HLT_PFHT400_Mu5_PFMET45",     "0", "2", ""),
        event_passing_HLT_PFHT400_Mu5_PFMET50      = cms.vstring("Event Passing HLT_PFHT400_Mu5_PFMET50",     "0", "2", ""),
        event_passing_HLT_PFNoPUHT350_Mu15_PFMET45 = cms.vstring("Event Passing HLT_PFNoPUHT350_Mu15_PFMET45","0", "2", ""),
        event_passing_HLT_PFNoPUHT350_Mu15_PFMET50 = cms.vstring("Event Passing HLT_PFNoPUHT350_Mu15_PFMET50","0", "2", ""),
        event_passing_HLT_PFNoPUHT400_Mu5_PFMET45  = cms.vstring("Event Passing HLT_PFNoPUHT400_Mu5_PFMET45", "0", "2", ""),
        event_passing_HLT_PFNoPUHT400_Mu5_PFMET50  = cms.vstring("Event Passing HLT_PFNoPUHT400_Mu5_PFMET50", "0", "2", ""),
    ),
    
    Categories = cms.PSet(
        #mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        passing_HLT_Mu5                      = cms.vstring("Probe Passing HLT_Mu5_v*",                     "dummy[pass=1,fail=0]"),
        passing_HLT_Mu8                      = cms.vstring("Probe Passing HLT_Mu8_v*",                     "dummy[pass=1,fail=0]"),
        passing_HLT_Mu12                     = cms.vstring("Probe Passing HLT_Mu12_v*",                    "dummy[pass=1,fail=0]"),
        passing_HLT_Mu17                     = cms.vstring("Probe Passing HLT_Mu17_v*",                    "dummy[pass=1,fail=0]"),
        passing_HLT_Mu24                     = cms.vstring("Probe Passing HLT_Mu24_v*",                    "dummy[pass=1,fail=0]"),
        passing_HLT_Mu30                     = cms.vstring("Probe Passing HLT_Mu30_v*",                    "dummy[pass=1,fail=0]"),
        passing_HLT_Mu40                     = cms.vstring("Probe Passing HLT_Mu40_v*",                    "dummy[pass=1,fail=0]"),
        passing_HLT_Mu15_eta2p1              = cms.vstring("Probe Passing HLT_Mu15_eta2p1_v*",             "dummy[pass=1,fail=0]"),
        passing_HLT_Mu24_eta2p1              = cms.vstring("Probe Passing HLT_Mu24_eta2p1_v*",             "dummy[pass=1,fail=0]"),
        passing_HLT_Mu30_eta2p1              = cms.vstring("Probe Passing HLT_Mu30_eta2p1_v*",             "dummy[pass=1,fail=0]"),
        passing_HLT_Mu40_eta2p1              = cms.vstring("Probe Passing HLT_Mu40_eta2p1_v*",             "dummy[pass=1,fail=0]"),
        passing_HLT_Mu50_eta2p1              = cms.vstring("Probe Passing HLT_Mu50_eta2p1_v*",             "dummy[pass=1,fail=0]"),
        passing_HLT_IsoMu24                  = cms.vstring("Probe Passing HLT_IsoMu24_v*",                 "dummy[pass=1,fail=0]"),
        passing_HLT_IsoMu30                  = cms.vstring("Probe Passing HLT_IsoMu30_v*",                 "dummy[pass=1,fail=0]"),
        passing_HLT_IsoMu20_eta2p1           = cms.vstring("Probe Passing HLT_IsoMu20_eta2p1_v*",          "dummy[pass=1,fail=0]"),
        passing_HLT_IsoMu24_eta2p1           = cms.vstring("Probe Passing HLT_IsoMu24_eta2p1_v*",          "dummy[pass=1,fail=0]"),
        passing_HLT_IsoMu30_eta2p1           = cms.vstring("Probe Passing HLT_IsoMu30_eta2p1_v*",          "dummy[pass=1,fail=0]"),
        passing_HLT_IsoMu34_eta2p1           = cms.vstring("Probe Passing HLT_IsoMu34_eta2p1_v*",          "dummy[pass=1,fail=0]"),
        passing_HLT_IsoMu40_eta2p1           = cms.vstring("Probe Passing HLT_IsoMu40_eta2p1_v*",          "dummy[pass=1,fail=0]"),
        passing_HLT_Mu40_HT200               = cms.vstring("Probe Passing HLT_Mu40_HT200_v*",              "dummy[pass=1,fail=0]"),
        passing_HLT_Mu40_FJHT200             = cms.vstring("Probe Passing HLT_Mu40_FJHT200_v*",            "dummy[pass=1,fail=0]"),
        passing_HLT_Mu40_PFHT350             = cms.vstring("Probe Passing HLT_Mu40_PFHT350_v*",            "dummy[pass=1,fail=0]"),
        passing_HLT_Mu60_PFHT350             = cms.vstring("Probe Passing HLT_Mu60_PFHT350_v*",            "dummy[pass=1,fail=0]"),
        passing_HLT_Mu40_PFNoPUHT350         = cms.vstring("Probe Passing HLT_Mu40_PFNoPUHT350_v*",        "dummy[pass=1,fail=0]"),
        passing_HLT_Mu60_PFNoPUHT350         = cms.vstring("Probe Passing HLT_Mu60_PFNoPUHT350_v*",        "dummy[pass=1,fail=0]"),
        passing_HLT_PFHT350_Mu15_PFMET45     = cms.vstring("Probe Passing HLT_PFHT350_Mu15_PFMET45_v*",    "dummy[pass=1,fail=0]"),
        passing_HLT_PFHT350_Mu15_PFMET50     = cms.vstring("Probe Passing HLT_PFHT350_Mu15_PFMET50_v*",    "dummy[pass=1,fail=0]"),
        passing_HLT_PFHT400_Mu5_PFMET45      = cms.vstring("Probe Passing HLT_PFHT400_Mu5_PFMET45_v*",     "dummy[pass=1,fail=0]"),
        passing_HLT_PFHT400_Mu5_PFMET50      = cms.vstring("Probe Passing HLT_PFHT400_Mu5_PFMET50_v*",     "dummy[pass=1,fail=0]"),
        passing_HLT_PFNoPUHT350_Mu15_PFMET45 = cms.vstring("Probe Passing HLT_PFNoPUHT350_Mu15_PFMET45_v*","dummy[pass=1,fail=0]"), 
        passing_HLT_PFNoPUHT350_Mu15_PFMET50 = cms.vstring("Probe Passing HLT_PFNoPUHT350_Mu15_PFMET50_v*","dummy[pass=1,fail=0]"), 
        passing_HLT_PFNoPUHT400_Mu5_PFMET45  = cms.vstring("Probe Passing HLT_PFNoPUHT400_Mu5_PFMET45_v*", "dummy[pass=1,fail=0]"),
        passing_HLT_PFNoPUHT400_Mu5_PFMET50  = cms.vstring("Probe Passing HLT_PFNoPUHT400_Mu5_PFMET50_v*", "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Mu5                      = cms.vstring("Tag Passing HLT_Mu5_v*",                     "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Mu8                      = cms.vstring("Tag Passing HLT_Mu8_v*",                     "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Mu12                     = cms.vstring("Tag Passing HLT_Mu12_v*",                    "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Mu17                     = cms.vstring("Tag Passing HLT_Mu17_v*",                    "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Mu24                     = cms.vstring("Tag Passing HLT_Mu24_v*",                    "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Mu30                     = cms.vstring("Tag Passing HLT_Mu30_v*",                    "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Mu40                     = cms.vstring("Tag Passing HLT_Mu40_v*",                    "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Mu15_eta2p1              = cms.vstring("Tag Passing HLT_Mu15_eta2p1_v*",             "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Mu24_eta2p1              = cms.vstring("Tag Passing HLT_Mu24_eta2p1_v*",             "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Mu30_eta2p1              = cms.vstring("Tag Passing HLT_Mu30_eta2p1_v*",             "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Mu40_eta2p1              = cms.vstring("Tag Passing HLT_Mu40_eta2p1_v*",             "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Mu50_eta2p1              = cms.vstring("Tag Passing HLT_Mu50_eta2p1_v*",             "dummy[pass=1,fail=0]"),
        tag_passing_HLT_IsoMu24                  = cms.vstring("Tag Passing HLT_IsoMu24_v*",                 "dummy[pass=1,fail=0]"),
        tag_passing_HLT_IsoMu30                  = cms.vstring("Tag Passing HLT_IsoMu30_v*",                 "dummy[pass=1,fail=0]"),
        tag_passing_HLT_IsoMu20_eta2p1           = cms.vstring("Tag Passing HLT_IsoMu20_eta2p1_v*",          "dummy[pass=1,fail=0]"),
        tag_passing_HLT_IsoMu24_eta2p1           = cms.vstring("Tag Passing HLT_IsoMu24_eta2p1_v*",          "dummy[pass=1,fail=0]"),
        tag_passing_HLT_IsoMu30_eta2p1           = cms.vstring("Tag Passing HLT_IsoMu30_eta2p1_v*",          "dummy[pass=1,fail=0]"),
        tag_passing_HLT_IsoMu34_eta2p1           = cms.vstring("Tag Passing HLT_IsoMu34_eta2p1_v*",          "dummy[pass=1,fail=0]"),
        tag_passing_HLT_IsoMu40_eta2p1           = cms.vstring("Tag Passing HLT_IsoMu40_eta2p1_v*",          "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Mu40_HT200               = cms.vstring("Tag Passing HLT_Mu40_HT200_v*",              "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Mu40_FJHT200             = cms.vstring("Tag Passing HLT_Mu40_FJHT200_v*",            "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Mu40_PFHT350             = cms.vstring("Tag Passing HLT_Mu40_PFHT350_v*",            "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Mu60_PFHT350             = cms.vstring("Tag Passing HLT_Mu60_PFHT350_v*",            "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Mu40_PFNoPUHT350         = cms.vstring("Tag Passing HLT_Mu40_PFNoPUHT350_v*",        "dummy[pass=1,fail=0]"),
        tag_passing_HLT_Mu60_PFNoPUHT350         = cms.vstring("Tag Passing HLT_Mu60_PFNoPUHT350_v*",        "dummy[pass=1,fail=0]"),
        tag_passing_HLT_PFHT350_Mu15_PFMET45     = cms.vstring("Tag Passing HLT_PFHT350_Mu15_PFMET45_v*",    "dummy[pass=1,fail=0]"),
        tag_passing_HLT_PFHT350_Mu15_PFMET50     = cms.vstring("Tag Passing HLT_PFHT350_Mu15_PFMET50_v*",    "dummy[pass=1,fail=0]"),
        tag_passing_HLT_PFHT400_Mu5_PFMET45      = cms.vstring("Tag Passing HLT_PFHT400_Mu5_PFMET45_v*",     "dummy[pass=1,fail=0]"),
        tag_passing_HLT_PFHT400_Mu5_PFMET50      = cms.vstring("Tag Passing HLT_PFHT400_Mu5_PFMET50_v*",     "dummy[pass=1,fail=0]"),
        tag_passing_HLT_PFNoPUHT350_Mu15_PFMET45 = cms.vstring("Tag Passing HLT_PFNoPUHT350_Mu15_PFMET45_v*","dummy[pass=1,fail=0]"), 
        tag_passing_HLT_PFNoPUHT350_Mu15_PFMET50 = cms.vstring("Tag Passing HLT_PFNoPUHT350_Mu15_PFMET50_v*","dummy[pass=1,fail=0]"), 
        tag_passing_HLT_PFNoPUHT400_Mu5_PFMET45  = cms.vstring("Tag Passing HLT_PFNoPUHT400_Mu5_PFMET45_v*", "dummy[pass=1,fail=0]"),
        tag_passing_HLT_PFNoPUHT400_Mu5_PFMET50  = cms.vstring("Tag Passing HLT_PFNoPUHT400_Mu5_PFMET50_v*", "dummy[pass=1,fail=0]"),
    ),
)

#################################################################

#  ____   _                _               
# |  _ \ (_)              (_)              
# | |_) | _  _ __   _ __   _  _ __    __ _ 
# |  _ < | || '_ \ | '_ \ | || '_ \  / _` |
# | |_) || || | | || | | || || | | || (_| |
# |____/ |_||_| |_||_| |_||_||_| |_| \__, |
#                                     __/ |
#                                    |___/ 

# Common Bins
CommonBins = cms.PSet(
    pt       = cms.vdouble( 10, 20, 30, 40, 50, 65, 80, 200 ),
    nvtx     = cms.vdouble( 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50 ),
    pileup   = cms.vdouble( 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50 ),
    instlumi = cms.vdouble( 0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000 ),
    njet     = cms.vdouble( 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5 ),
    drjet    = cms.vdouble( 0.0, 0.3, 0.6, 0.9, 1.3, 1.7, 2.1, 2.5, 3, 3.5, 4, 4.5),
    reliso   = cms.vdouble( 0, 0.01, 0.02, 0.04, 0.06, 0.09, 0.12, 0.15, 0.2, 0.25),
    ht       = cms.vdouble( 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 1000),
    met      = cms.vdouble( 0, 10, 20, 30, 40, 60, 150),
    st       = cms.vdouble( 0,  25, 50, 75, 100, 150, 200, 250, 300, 400, 500 ),
)
pt5Bins  = cms.PSet( pt = cms.vdouble(  0,  2,  4,  6,  8, 12, 20, 30, 45, 70, 100) )
pt15Bins = cms.PSet( pt = cms.vdouble(     10, 12, 14, 16, 20, 25, 35, 50, 70, 100) )
pt24Bins = cms.PSet( pt = cms.vdouble(     15, 18, 22, 26, 30, 35, 40, 55, 75, 100) )
pt27Bins = cms.PSet( pt = cms.vdouble(         20, 22, 25, 29, 35, 40, 55, 75, 100) )
pt30Bins = cms.PSet( pt = cms.vdouble(         20, 24, 28, 32, 36, 40, 55, 75, 100) )
pt40Bins = cms.PSet( pt = cms.vdouble(         30, 34, 38, 42, 46, 50, 60, 80, 100) )
pt80Bins = cms.PSet( pt = cms.vdouble( 70, 74, 78, 82, 86, 90, 100, 125, 150) )

RelIso001Bins = cms.PSet( reliso = cms.vdouble(0, 0.01) )
RelIso002Bins = cms.PSet( reliso = cms.vdouble(0, 0.02) )
RelIso004Bins = cms.PSet( reliso = cms.vdouble(0, 0.04) )
RelIso006Bins = cms.PSet( reliso = cms.vdouble(0, 0.06) )
RelIso009Bins = cms.PSet( reliso = cms.vdouble(0, 0.09) )
RelIso012Bins = cms.PSet( reliso = cms.vdouble(0, 0.12) )
RelIso015Bins = cms.PSet( reliso = cms.vdouble(0, 0.15) )
RelIso020Bins = cms.PSet( reliso = cms.vdouble(0, 0.20) )
RelIso025Bins = cms.PSet( reliso = cms.vdouble(0, 0.25) )

# Electron Bins
EleBins = cms.PSet(
    eta = cms.vdouble( -2.5,  -2.1, -1.6,  -1.1,      -0.6,          0.0,         0.6,     1.1,  1.6, 2.1,  2.5 ),
)
EleGsfBins = cms.PSet( pt  = cms.vdouble( 20, 200 ), eta = cms.vdouble( -2.5, 2.5 ) )
EleIdBins  = cms.PSet( pt  = cms.vdouble( 20, 200 ), abs_eta = cms.vdouble( 0, 1.5, 2.5 ) ) # 1.479 is the EE-EB Gap
EleHLTBins = cms.PSet(
    eta = cms.vdouble(-2.5, 2.5),
    d0_v = cms.vdouble(-0.02, 0.02),
    dz_v = cms.vdouble(-0.1, 0.1),
    drjet = cms.vdouble(0.3, 10),
    absdeltapt = cms.vdouble(-10000, 10.0),
    reliso = cms.vdouble(0, 0.15),
)
EleHLTpt15Bins  = EleHLTBins.clone( pt = cms.vdouble( 15, 200) )
EleHLTpt20Bins  = EleHLTBins.clone( pt = cms.vdouble( 20, 200) )
EleHLTpt40Bins  = EleHLTBins.clone( pt = cms.vdouble( 40, 200) )
EleHLTpt50Bins  = EleHLTBins.clone( pt = cms.vdouble( 50, 200) )
EleHLTpt100Bins = EleHLTBins.clone( pt = cms.vdouble(100, 200) )

# Muon Bins
MuBins  = cms.PSet(
    eta = cms.vdouble(  -2.4, -2.1, -1.6, -1.2, -0.9, -0.6, -0.3, -0.2, 0.2, 0.3, 0.6, 0.9, 1.2, 1.6, 2.1, 2.4 ),
)
MuTrkBins  = cms.PSet( pt = cms.vdouble( 20, 200 ), eta = cms.vdouble( -2.4, 2.4 ) )
MuIdBins  = MuTrkBins
MuIsoBins  = MuIdBins.clone(
    passing_pat = cms.vstring("pass"),
    d0_v = cms.vdouble(-0.02, 0.02),
    dz_v = cms.vdouble(-0.5, 0.5),
    drjet = cms.vdouble(0.3, 10),
    absdeltapt = cms.vdouble(-10000, 5.0),
)
MuHLTBins = cms.PSet(
    eta = cms.vdouble(-2.4, 2.4),
    d0_v = cms.vdouble(-0.02, 0.02),
    dz_v = cms.vdouble(-0.5, 0.5),
    drjet = cms.vdouble(0.3, 10),
    absdeltapt = cms.vdouble(-10000, 5.0),
    reliso = cms.vdouble(0, 0.12),
)
MuHLTpt10Bins = MuHLTBins.clone( pt = cms.vdouble(10, 200) )
MuHLTpt20Bins = MuHLTBins.clone( pt = cms.vdouble(20, 200) )
MuHLTpt30Bins = MuHLTBins.clone( pt = cms.vdouble(30, 200) )
MuHLTpt45Bins = MuHLTBins.clone( pt = cms.vdouble(45, 200) )
MuHLTpt30eta2p1Bins = MuHLTBins.clone( pt = cms.vdouble(30, 200), eta = cms.vdouble(-2.1, 2.1) )

#################################################################

#  ______   __   __  _        _                     _            
# |  ____| / _| / _|(_)      (_)                   (_)           
# | |__   | |_ | |_  _   ___  _   ___  _ __    ___  _   ___  ___ 
# |  __|  |  _||  _|| | / __|| | / _ \| '_ \  / __|| | / _ \/ __|
# | |____ | |  | |  | || (__ | ||  __/| | | || (__ | ||  __/\__ \
# |______||_|  |_|  |_| \___||_| \___||_| |_| \___||_| \___||___/

EffTemplate = cms.PSet(
    UnbinnedVariables = theUnbinned,
    BinToPDFmap = cms.vstring("vpvPlusExpo"),
)

# Electron - Gsf Efficiency
EleGsf_avg = EffTemplate.clone(
    BinnedVariables = EleGsfBins,
    EfficiencyCategoryAndState = cms.vstring("passing_gsf","pass"),
)

# Electron - Id Efficiency
EleId_avg = EffTemplate.clone(
    BinnedVariables = EleIdBins,
    EfficiencyCategoryAndState = cms.vstring(
        "passing_pat","pass",
        "d0_v_m002","above","d0_v_p002","below",
        "dz_v_m01", "above","dz_v_p01", "below",
        "absdeltapt_m10000","above","absdeltapt_p10","below",
        "drjet_03","above",
        "reliso_015","below"
    ),
)

# Electron - Trigger Efficiency
HLT_Ele27_avg = EffTemplate.clone(
    EfficiencyCategoryAndState = cms.vstring("passing_HLT_Ele27_CaloIdL_CaloIsoVL", "pass"),
    BinnedVariables = EleHLTpt40Bins.clone( tag_passing_HLT_Ele27_CaloIdL_CaloIsoVL = cms.vstring("pass") ),
)
HLT_Ele27_WP80_avg = EffTemplate.clone(
    EfficiencyCategoryAndState = cms.vstring("passing_HLT_Ele27_WP80", "pass"),
    BinnedVariables = EleHLTpt40Bins.clone( tag_passing_HLT_Ele27_WP80 = cms.vstring("pass") ),
)
HLT_Ele30_avg = EffTemplate.clone(
    EfficiencyCategoryAndState = cms.vstring("passing_HLT_Ele30_CaloIdVT_TrkIdT", "pass"),
    BinnedVariables = EleHLTpt40Bins.clone( tag_passing_HLT_Ele30_CaloIdVT_TrkIdT = cms.vstring("pass") ),
)
HLT_Ele80_avg = EffTemplate.clone(
    EfficiencyCategoryAndState = cms.vstring("passing_HLT_Ele80_CaloIdVT_GsfTrkIdT", "pass"),
    BinnedVariables = EleHLTpt100Bins.clone( tag_passing_HLT_Ele90_CaloIdVT_GsfTrkIdT = cms.vstring("pass") ),
)
HLT_CleanPFHT350_Ele5_PFMET45_avg = EffTemplate.clone(
    EfficiencyCategoryAndState = cms.vstring("passing_HLT_CleanPFHT350_Ele5_PFMET45", "pass"),
    BinnedVariables = EleHLTpt15Bins.clone(
        tag_passing_HLT_CleanPFHT350_Ele5_PFMET45 = cms.vstring("pass"),   
        event_passing_HLT_CleanPFHT350_Ele5_PFMET45 = cms.vdouble(0.5, 1.5),                               
    ),
)
HLT_CleanPFHT300_Ele15_PFMET45_avg = EffTemplate.clone(
    EfficiencyCategoryAndState = cms.vstring("passing_HLT_CleanPFHT300_Ele15_PFMET45", "pass"),
    BinnedVariables = EleHLTpt20Bins.clone(
        tag_passing_HLT_CleanPFHT300_Ele15_PFMET45 = cms.vstring("pass"),
        event_passing_HLT_CleanPFHT300_Ele15_PFMET45 = cms.vdouble(0.5, 1.5),
    ),
)
HLT_CleanPFHT300_Ele40_avg = EffTemplate.clone(
    EfficiencyCategoryAndState = cms.vstring("passing_HLT_CleanPFHT300_Ele40", "pass"),
    BinnedVariables = EleHLTpt50Bins.clone(
        tag_passing_HLT_CleanPFHT300_Ele40 = cms.vstring("pass"),   
        event_passing_HLT_CleanPFHT300_Ele40 = cms.vdouble(0.5, 1.5),                               
    ),
)

# Muon - Tracking Efficiency
MuTrk_avg = EffTemplate.clone(
    BinnedVariables = MuTrkBins,
    EfficiencyCategoryAndState = cms.vstring("passing_trk","pass"),
)

# Muon - Id Efficiency (all cuts except isolation)
MuId_avg = EffTemplate.clone(
    BinnedVariables = MuIdBins,
    EfficiencyCategoryAndState = cms.vstring(
        "passing_glb","pass",
        "d0_v_m002","above","d0_v_p002","below",
        "dz_v_m05", "above","dz_v_p05", "below",
        "absdeltapt_m10000","above","absdeltapt_p5","below",
        "drjet_03","above",
        "reliso_012","below"
    ),
)

# Muon - Isolation efficiency
MuIso_avg = EffTemplate.clone(
    BinnedVariables = MuIsoBins,
    EfficiencyCategoryAndState = cms.vstring("reliso_012","below"),
)

# Muon - Trigger Efficiency
HLT_IsoMu24_avg = EffTemplate.clone(
    EfficiencyCategoryAndState = cms.vstring("passing_HLT_IsoMu24", "pass"),
    BinnedVariables = MuHLTpt30Bins.clone( tag_passing_HLT_IsoMu24 = cms.vstring("pass") ),
)
HLT_IsoMu24_eta2p1_avg = EffTemplate.clone(
    EfficiencyCategoryAndState = cms.vstring("passing_HLT_IsoMu24_eta2p1", "pass"),
    BinnedVariables = MuHLTpt30eta2p1Bins.clone( tag_passing_HLT_IsoMu24_eta2p1 = cms.vstring("pass") ),
)
HLT_Mu40_PFHT350_avg = EffTemplate.clone(
    EfficiencyCategoryAndState = cms.vstring("passing_HLT_Mu40_PFHT350", "pass"),
    BinnedVariables = MuHLTpt45Bins.clone(
        tag_passing_HLT_Mu40_PFHT350 = cms.vstring("pass"),   
        event_passing_HLT_Mu40_PFHT350 = cms.vdouble(0.5, 1.5),                               
    ),
)
HLT_PFHT400_Mu5_PFMET45_avg = EffTemplate.clone(
    EfficiencyCategoryAndState = cms.vstring("passing_HLT_PFHT400_Mu5_PFMET45", "pass"),
    BinnedVariables = MuHLTpt10Bins.clone(
        tag_passing_HLT_PFHT400_Mu5_PFMET45 = cms.vstring("pass"),   
        event_passing_HLT_PFHT400_Mu5_PFMET45 = cms.vdouble(0.5, 1.5),                               
    ),
)
HLT_PFHT350_Mu15_PFMET45_avg = EffTemplate.clone(
    EfficiencyCategoryAndState = cms.vstring("passing_HLT_PFHT350_Mu15_PFMET45", "pass"),
    BinnedVariables = MuHLTpt20Bins.clone(
        tag_passing_HLT_PFHT350_Mu15_PFMET45 = cms.vstring("pass"),   
        event_passing_HLT_PFHT350_Mu15_PFMET45 = cms.vdouble(0.5, 1.5),                               
    ),
)
HLT_Mu40_PFHT350_avg = EffTemplate.clone(
    EfficiencyCategoryAndState = cms.vstring("passing_HLT_Mu40_PFHT350", "pass"),
    BinnedVariables = MuHLTpt45Bins.clone(
        tag_passing_HLT_Mu40_PFHT350 = cms.vstring("pass"),   
        event_passing_HLT_Mu40_PFHT350 = cms.vdouble(0.5, 1.5),                               
    ),
)

#################################################################

#  _____   _         _        
# |  __ \ | |       | |       
# | |__) || |  ___  | |_  ___ 
# |  ___/ | | / _ \ | __|/ __|
# | |     | || (_) || |_ \__ \
# |_|     |_| \___/  \__||___/

#    ____   ____       __     ____      __ 
#   / ___| / ___|      \ \   / ___|___ / _|
#   \___ \| |      _____\ \ | |  _/ __| |_ 
#    ___) | |___  |_____/ / | |_| \__ \  _|
#   |____/ \____|      /_/   \____|___/_|  

process.FitEleGsfEff.Efficiencies = cms.PSet(
    avg      = EleGsf_avg,
    #eta      = EleGsf_avg.clone( BinnedVariables = EleGsfBins.clone( eta      = EleBins.eta ) ),
    #pt       = EleGsf_avg.clone( BinnedVariables = EleGsfBins.clone( pt       = CommonBins.pt ) ),
    #nvtx     = EleGsf_avg.clone( BinnedVariables = EleGsfBins.clone( nvtx     = CommonBins.nvtx ) ),
    #pileup   = EleGsf_avg.clone( BinnedVariables = EleGsfBins.clone( pileup   = CommonBins.pileup ) ),
    #instlumi = EleGsf_avg.clone( BinnedVariables = EleGsfBins.clone( instlumi = CommonBins.instlumi ) ),
)

#   ____      __       __     ___    _ 
#  / ___|___ / _|      \ \   |_ _|__| |
# | |  _/ __| |_   _____\ \   | |/ _` |
# | |_| \__ \  _| |_____/ /   | | (_| |
#  \____|___/_|        /_/   |___\__,_|
#

process.FitEleIdEff.Efficiencies = cms.PSet(
    avg          = EleId_avg,
    #eta          = EleId_avg.clone( BinnedVariables = cms.PSet( pt  = cms.vdouble( 20, 200 ), eta = EleBins.eta ) ),
    #pt           = EleId_avg.clone( BinnedVariables = EleIdBins.clone( pt       = CommonBins.pt ) ),
    #nvtx         = EleId_avg.clone( BinnedVariables = EleIdBins.clone( nvtx     = CommonBins.nvtx ) ),
    #pileup       = EleId_avg.clone( BinnedVariables = EleIdBins.clone( pileup   = CommonBins.pileup ) ),
    #instlumi     = EleId_avg.clone( BinnedVariables = EleIdBins.clone( instlumi = CommonBins.instlumi ) ),
    #njet         = EleId_avg.clone( BinnedVariables = EleIdBins.clone( njet     = CommonBins.njet ) ),
)

#    ___    _       __    _   _ _   _____ 
#   |_ _|__| |      \ \  | | | | | |_   _|
#    | |/ _` |  _____\ \ | |_| | |   | |  
#    | | (_| | |_____/ / |  _  | |___| |  
#   |___\__,_|      /_/  |_| |_|_____|_|

process.FitEleHLTEff.Efficiencies = cms.PSet(
    HLT_Ele27_avg      = HLT_Ele27_avg,
    #HLT_Ele27_pt       = HLT_Ele27_avg.clone( BinnedVariables = HLT_Ele27_avg.BinnedVariables.clone( pt       = pt27Bins.pt ) ),
    #HLT_Ele27_eta      = HLT_Ele27_avg.clone( BinnedVariables = HLT_Ele27_avg.BinnedVariables.clone( eta      = EleBins.eta ) ),
    #HLT_Ele27_reliso   = HLT_Ele27_avg.clone( BinnedVariables = HLT_Ele27_avg.BinnedVariables.clone( reliso   = CommonBins.reliso ) ),
    #HLT_Ele27_nvtx     = HLT_Ele27_avg.clone( BinnedVariables = HLT_Ele27_avg.BinnedVariables.clone( nvtx     = CommonBins.nvtx ) ),
    #HLT_Ele27_pileup   = HLT_Ele27_avg.clone( BinnedVariables = HLT_Ele27_avg.BinnedVariables.clone( pileup   = CommonBins.pileup, ) ),
    #HLT_Ele27_instlumi = HLT_Ele27_avg.clone( BinnedVariables = HLT_Ele27_avg.BinnedVariables.clone( instlumi = CommonBins.instlumi ) ),
    #HLT_Ele27_njet     = HLT_Ele27_avg.clone( BinnedVariables = HLT_Ele27_avg.BinnedVariables.clone( njet     = CommonBins.njet ) ),
    #HLT_Ele27_drjet    = HLT_Ele27_avg.clone( BinnedVariables = HLT_Ele27_avg.BinnedVariables.clone( drjet    = CommonBins.drjet ) ),
    #HLT_Ele27_ht       = HLT_Ele27_avg.clone( BinnedVariables = HLT_Ele27_avg.BinnedVariables.clone( ht       = CommonBins.ht ) ),
    #HLT_Ele27_met      = HLT_Ele27_avg.clone( BinnedVariables = HLT_Ele27_avg.BinnedVariables.clone( met      = CommonBins.met ) ),
    #HLT_Ele27_st       = HLT_Ele27_avg.clone( BinnedVariables = HLT_Ele27_avg.BinnedVariables.clone( st       = CommonBins.st ) ),
    #
    #HLT_Ele27_WP80_avg      = HLT_Ele27_WP80_avg,
    #HLT_Ele27_WP80_pt       = HLT_Ele27_WP80_avg.clone( BinnedVariables = HLT_Ele27_WP80_avg.BinnedVariables.clone( pt       = pt27Bins.pt ) ),
    #HLT_Ele27_WP80_eta      = HLT_Ele27_WP80_avg.clone( BinnedVariables = HLT_Ele27_WP80_avg.BinnedVariables.clone( eta      = EleBins.eta ) ),
    #HLT_Ele27_WP80_reliso   = HLT_Ele27_WP80_avg.clone( BinnedVariables = HLT_Ele27_WP80_avg.BinnedVariables.clone( reliso   = CommonBins.reliso ) ),
    #HLT_Ele27_WP80_nvtx     = HLT_Ele27_WP80_avg.clone( BinnedVariables = HLT_Ele27_WP80_avg.BinnedVariables.clone( nvtx     = CommonBins.nvtx ) ),
    #HLT_Ele27_WP80_pileup   = HLT_Ele27_WP80_avg.clone( BinnedVariables = HLT_Ele27_WP80_avg.BinnedVariables.clone( pileup   = CommonBins.pileup, ) ),
    #HLT_Ele27_WP80_instlumi = HLT_Ele27_WP80_avg.clone( BinnedVariables = HLT_Ele27_WP80_avg.BinnedVariables.clone( instlumi = CommonBins.instlumi ) ),
    #HLT_Ele27_WP80_njet     = HLT_Ele27_WP80_avg.clone( BinnedVariables = HLT_Ele27_WP80_avg.BinnedVariables.clone( njet     = CommonBins.njet ) ),
    #HLT_Ele27_WP80_drjet    = HLT_Ele27_WP80_avg.clone( BinnedVariables = HLT_Ele27_WP80_avg.BinnedVariables.clone( drjet    = CommonBins.drjet ) ),
    #HLT_Ele27_WP80_ht       = HLT_Ele27_WP80_avg.clone( BinnedVariables = HLT_Ele27_WP80_avg.BinnedVariables.clone( ht       = CommonBins.ht ) ),
    #HLT_Ele27_WP80_met      = HLT_Ele27_WP80_avg.clone( BinnedVariables = HLT_Ele27_WP80_avg.BinnedVariables.clone( met      = CommonBins.met ) ),
    #HLT_Ele27_WP80_st       = HLT_Ele27_WP80_avg.clone( BinnedVariables = HLT_Ele27_WP80_avg.BinnedVariables.clone( st       = CommonBins.st ) ),
    #
    HLT_Ele30_avg      = HLT_Ele30_avg,
    #HLT_Ele30_pt       = HLT_Ele30_avg.clone( BinnedVariables = HLT_Ele30_avg.BinnedVariables.clone( pt       = pt30Bins.pt ) ),
    #HLT_Ele30_eta      = HLT_Ele30_avg.clone( BinnedVariables = HLT_Ele30_avg.BinnedVariables.clone( eta      = EleBins.eta ) ),
    #HLT_Ele30_reliso   = HLT_Ele30_avg.clone( BinnedVariables = HLT_Ele30_avg.BinnedVariables.clone( reliso   = CommonBins.reliso ) ),
    #HLT_Ele30_nvtx     = HLT_Ele30_avg.clone( BinnedVariables = HLT_Ele30_avg.BinnedVariables.clone( nvtx     = CommonBins.nvtx ) ),
    #HLT_Ele30_pileup   = HLT_Ele30_avg.clone( BinnedVariables = HLT_Ele30_avg.BinnedVariables.clone( pileup   = CommonBins.pileup, ) ),
    #HLT_Ele30_instlumi = HLT_Ele30_avg.clone( BinnedVariables = HLT_Ele30_avg.BinnedVariables.clone( instlumi = CommonBins.instlumi ) ),
    #HLT_Ele30_njet     = HLT_Ele30_avg.clone( BinnedVariables = HLT_Ele30_avg.BinnedVariables.clone( njet     = CommonBins.njet ) ),
    #HLT_Ele30_drjet    = HLT_Ele30_avg.clone( BinnedVariables = HLT_Ele30_avg.BinnedVariables.clone( drjet    = CommonBins.drjet ) ),
    #HLT_Ele30_ht       = HLT_Ele30_avg.clone( BinnedVariables = HLT_Ele30_avg.BinnedVariables.clone( ht       = CommonBins.ht ) ),
    #HLT_Ele30_met      = HLT_Ele30_avg.clone( BinnedVariables = HLT_Ele30_avg.BinnedVariables.clone( met      = CommonBins.met ) ),
    #HLT_Ele30_st       = HLT_Ele30_avg.clone( BinnedVariables = HLT_Ele30_avg.BinnedVariables.clone( st       = CommonBins.st ) ),
    #
    HLT_Ele80_avg      = HLT_Ele80_avg,
    #HLT_Ele80_pt       = HLT_Ele80_avg.clone( BinnedVariables = HLT_Ele80_avg.BinnedVariables.clone( pt       = pt80Bins.pt ) ),
    #HLT_Ele80_eta      = HLT_Ele80_avg.clone( BinnedVariables = HLT_Ele80_avg.BinnedVariables.clone( eta      = EleBins.eta ) ),
    #HLT_Ele80_reliso   = HLT_Ele80_avg.clone( BinnedVariables = HLT_Ele80_avg.BinnedVariables.clone( reliso   = CommonBins.reliso ) ),
    #HLT_Ele80_nvtx     = HLT_Ele80_avg.clone( BinnedVariables = HLT_Ele80_avg.BinnedVariables.clone( nvtx     = CommonBins.nvtx ) ),
    #HLT_Ele80_pileup   = HLT_Ele80_avg.clone( BinnedVariables = HLT_Ele80_avg.BinnedVariables.clone( pileup   = CommonBins.pileup, ) ),
    #HLT_Ele80_instlumi = HLT_Ele80_avg.clone( BinnedVariables = HLT_Ele80_avg.BinnedVariables.clone( instlumi = CommonBins.instlumi ) ),
    #HLT_Ele80_njet     = HLT_Ele80_avg.clone( BinnedVariables = HLT_Ele80_avg.BinnedVariables.clone( njet     = CommonBins.njet ) ),
    #HLT_Ele80_drjet    = HLT_Ele80_avg.clone( BinnedVariables = HLT_Ele80_avg.BinnedVariables.clone( drjet    = CommonBins.drjet ) ),
    #HLT_Ele80_ht       = HLT_Ele80_avg.clone( BinnedVariables = HLT_Ele80_avg.BinnedVariables.clone( ht       = CommonBins.ht ) ),
    #HLT_Ele80_met      = HLT_Ele80_avg.clone( BinnedVariables = HLT_Ele80_avg.BinnedVariables.clone( met      = CommonBins.met ) ),
    #HLT_Ele80_st       = HLT_Ele80_avg.clone( BinnedVariables = HLT_Ele80_avg.BinnedVariables.clone( st       = CommonBins.st ) ),
    #
    HLT_CleanPFHT350_Ele5_PFMET45_avg      = HLT_CleanPFHT350_Ele5_PFMET45_avg,
    #HLT_CleanPFHT350_Ele5_PFMET45_pt       = HLT_CleanPFHT350_Ele5_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT350_Ele5_PFMET45_avg.BinnedVariables.clone( pt       = pt5Bins.pt ) ),
    #HLT_CleanPFHT350_Ele5_PFMET45_eta      = HLT_CleanPFHT350_Ele5_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT350_Ele5_PFMET45_avg.BinnedVariables.clone( eta      = EleBins.eta ) ),
    #HLT_CleanPFHT350_Ele5_PFMET45_reliso   = HLT_CleanPFHT350_Ele5_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT350_Ele5_PFMET45_avg.BinnedVariables.clone( reliso   = CommonBins.reliso ) ),
    #HLT_CleanPFHT350_Ele5_PFMET45_nvtx     = HLT_CleanPFHT350_Ele5_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT350_Ele5_PFMET45_avg.BinnedVariables.clone( nvtx     = CommonBins.nvtx ) ),
    #HLT_CleanPFHT350_Ele5_PFMET45_pileup   = HLT_CleanPFHT350_Ele5_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT350_Ele5_PFMET45_avg.BinnedVariables.clone( pileup   = CommonBins.pileup, ) ),
    #HLT_CleanPFHT350_Ele5_PFMET45_instlumi = HLT_CleanPFHT350_Ele5_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT350_Ele5_PFMET45_avg.BinnedVariables.clone( instlumi = CommonBins.instlumi ) ),
    #HLT_CleanPFHT350_Ele5_PFMET45_njet     = HLT_CleanPFHT350_Ele5_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT350_Ele5_PFMET45_avg.BinnedVariables.clone( njet     = CommonBins.njet ) ),
    #HLT_CleanPFHT350_Ele5_PFMET45_drjet    = HLT_CleanPFHT350_Ele5_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT350_Ele5_PFMET45_avg.BinnedVariables.clone( drjet    = CommonBins.drjet ) ),
    #HLT_CleanPFHT350_Ele5_PFMET45_ht       = HLT_CleanPFHT350_Ele5_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT350_Ele5_PFMET45_avg.BinnedVariables.clone( ht       = CommonBins.ht ) ),
    #HLT_CleanPFHT350_Ele5_PFMET45_met      = HLT_CleanPFHT350_Ele5_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT350_Ele5_PFMET45_avg.BinnedVariables.clone( met      = CommonBins.met ) ),
    #HLT_CleanPFHT350_Ele5_PFMET45_st       = HLT_CleanPFHT350_Ele5_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT350_Ele5_PFMET45_avg.BinnedVariables.clone( st       = CommonBins.st ) ),
    #
    HLT_CleanPFHT300_Ele15_PFMET45_avg      = HLT_CleanPFHT300_Ele15_PFMET45_avg,
    #HLT_CleanPFHT300_Ele15_PFMET45_pt       = HLT_CleanPFHT300_Ele15_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele15_PFMET45_avg.BinnedVariables.clone( pt       = pt15Bins.pt ) ),
    #HLT_CleanPFHT300_Ele15_PFMET45_eta      = HLT_CleanPFHT300_Ele15_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele15_PFMET45_avg.BinnedVariables.clone( eta      = EleBins.eta ) ),
    #HLT_CleanPFHT300_Ele15_PFMET45_reliso   = HLT_CleanPFHT300_Ele15_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele15_PFMET45_avg.BinnedVariables.clone( reliso   = CommonBins.reliso ) ),
    #HLT_CleanPFHT300_Ele15_PFMET45_nvtx     = HLT_CleanPFHT300_Ele15_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele15_PFMET45_avg.BinnedVariables.clone( nvtx     = CommonBins.nvtx ) ),
    #HLT_CleanPFHT300_Ele15_PFMET45_pileup   = HLT_CleanPFHT300_Ele15_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele15_PFMET45_avg.BinnedVariables.clone( pileup   = CommonBins.pileup, ) ),
    #HLT_CleanPFHT300_Ele15_PFMET45_instlumi = HLT_CleanPFHT300_Ele15_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele15_PFMET45_avg.BinnedVariables.clone( instlumi = CommonBins.instlumi ) ),
    #HLT_CleanPFHT300_Ele15_PFMET45_njet     = HLT_CleanPFHT300_Ele15_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele15_PFMET45_avg.BinnedVariables.clone( njet     = CommonBins.njet ) ),
    #HLT_CleanPFHT300_Ele15_PFMET45_drjet    = HLT_CleanPFHT300_Ele15_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele15_PFMET45_avg.BinnedVariables.clone( drjet    = CommonBins.drjet ) ),
    #HLT_CleanPFHT300_Ele15_PFMET45_ht       = HLT_CleanPFHT300_Ele15_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele15_PFMET45_avg.BinnedVariables.clone( ht       = CommonBins.ht ) ),
    #HLT_CleanPFHT300_Ele15_PFMET45_met      = HLT_CleanPFHT300_Ele15_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele15_PFMET45_avg.BinnedVariables.clone( met      = CommonBins.met ) ),
    #HLT_CleanPFHT300_Ele15_PFMET45_st       = HLT_CleanPFHT300_Ele15_PFMET45_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele15_PFMET45_avg.BinnedVariables.clone( st       = CommonBins.st ) ),
    #
    HLT_CleanPFHT300_Ele40_avg      = HLT_CleanPFHT300_Ele40_avg,
    #HLT_CleanPFHT300_Ele40_pt       = HLT_CleanPFHT300_Ele40_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele40_avg.BinnedVariables.clone( pt       = pt40Bins.pt ) ),
    #HLT_CleanPFHT300_Ele40_eta      = HLT_CleanPFHT300_Ele40_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele40_avg.BinnedVariables.clone( eta      = EleBins.eta ) ),
    #HLT_CleanPFHT300_Ele40_reliso   = HLT_CleanPFHT300_Ele40_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele40_avg.BinnedVariables.clone( reliso   = CommonBins.reliso ) ),
    #HLT_CleanPFHT300_Ele40_nvtx     = HLT_CleanPFHT300_Ele40_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele40_avg.BinnedVariables.clone( nvtx     = CommonBins.nvtx ) ),
    #HLT_CleanPFHT300_Ele40_pileup   = HLT_CleanPFHT300_Ele40_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele40_avg.BinnedVariables.clone( pileup   = CommonBins.pileup, ) ),
    #HLT_CleanPFHT300_Ele40_instlumi = HLT_CleanPFHT300_Ele40_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele40_avg.BinnedVariables.clone( instlumi = CommonBins.instlumi ) ),
    #HLT_CleanPFHT300_Ele40_njet     = HLT_CleanPFHT300_Ele40_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele40_avg.BinnedVariables.clone( njet     = CommonBins.njet ) ),
    #HLT_CleanPFHT300_Ele40_drjet    = HLT_CleanPFHT300_Ele40_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele40_avg.BinnedVariables.clone( drjet    = CommonBins.drjet ) ),
    #HLT_CleanPFHT300_Ele40_ht       = HLT_CleanPFHT300_Ele40_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele40_avg.BinnedVariables.clone( ht       = CommonBins.ht ) ),
    #HLT_CleanPFHT300_Ele40_met      = HLT_CleanPFHT300_Ele40_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele40_avg.BinnedVariables.clone( met      = CommonBins.met ) ),
    #HLT_CleanPFHT300_Ele40_st       = HLT_CleanPFHT300_Ele40_avg.clone( BinnedVariables = HLT_CleanPFHT300_Ele40_avg.BinnedVariables.clone( st       = CommonBins.st ) ),
)

################################################################

#   _____ _              __     _______   _    
#  / ____| |             \ \   |__   __| | |   
# | (___ | |_ __ _   _____\ \     | |_ __| | __
#  \___ \| __/ _` | |______> >    | | '__| |/ /
#  ____) | || (_| |       / /     | | |  |   < 
# |_____/ \__\__,_|      /_/      |_|_|  |_|\_\

process.FitMuTrkEff.Efficiencies = cms.PSet(
    avg      = MuTrk_avg,
    #eta      = MuTrk_avg.clone( BinnedVariables = MuTrkBins.clone( eta      = MuBins.eta ) ),
    #pt       = MuTrk_avg.clone( BinnedVariables = MuTrkBins.clone( pt       = CommonBins.pt ) ),
    #nvtx     = MuTrk_avg.clone( BinnedVariables = MuTrkBins.clone( nvtx     = CommonBins.nvtx ) ),
    #pileup   = MuTrk_avg.clone( BinnedVariables = MuTrkBins.clone( pileup   = CommonBins.pileup ) ),
    #instlumi = MuTrk_avg.clone( BinnedVariables = MuTrkBins.clone( instlumi = CommonBins.instlumi ) ),
    #njet     = MuTrk_avg.clone( BinnedVariables = MuTrkBins.clone( njet     = CommonBins.njet ) ),
    #drjet    = MuTrk_avg.clone( BinnedVariables = MuTrkBins.clone( drjet    = CommonBins.drjet ) ),
)

#  _______   _          __     _____            _____    _ 
# |__   __| | |         \ \   |_   _|          |_   _|  | |
#    | |_ __| | __  _____\ \    | |  ___  ___    | |  __| |
#    | | '__| |/ / |______> >   | | / __|/ _ \   | | / _` |
#    | | |  |   <        / /   _| |_\__ \ (_) | _| || (_| |
#    |_|_|  |_|\_\      /_/   |_____|___/\___( )_____\__,_|
#                                            |/

process.FitMuIdEff.Efficiencies = cms.PSet(
    avg      = MuId_avg,
    #eta      = MuId_avg.clone( BinnedVariables = MuIdBins.clone( eta      = MuBins.eta ) ),
    #pt       = MuId_avg.clone( BinnedVariables = MuIdBins.clone( pt       = CommonBins.pt ) ),
    #nvtx     = MuId_avg.clone( BinnedVariables = MuIdBins.clone( nvtx     = CommonBins.nvtx ) ),
    #pileup   = MuId_avg.clone( BinnedVariables = MuIdBins.clone( pileup   = CommonBins.pileup ) ),
    #instlumi = MuId_avg.clone( BinnedVariables = MuIdBins.clone( instlumi = CommonBins.instlumi ) ),
    #iso_avg      = MuIso_avg,
    #iso_eta      = MuIso_avg.clone( BinnedVariables = MuIsoBins.clone( eta      = MuBins.eta ) ),
    #iso_pt       = MuIso_avg.clone( BinnedVariables = MuIsoBins.clone( pt       = CommonBins.pt ) ),
    #iso_nvtx     = MuIso_avg.clone( BinnedVariables = MuIsoBins.clone( nvtx     = CommonBins.nvtx ) ),
    #iso_pileup   = MuIso_avg.clone( BinnedVariables = MuIsoBins.clone( pileup   = CommonBins.pileup ) ),
    #iso_instlumi = MuIso_avg.clone( BinnedVariables = MuIsoBins.clone( instlumi = CommonBins.instlumi ) ),
    #iso_njet     = MuIso_avg.clone( BinnedVariables = MuIsoBins.clone( njet     = CommonBins.njet ) ),
)

#  _____    _       __     _    _ _   _______ 
# |_   _|  | |      \ \   | |  | | | |__   __|
#   | |  __| |  _____\ \  | |__| | |    | |   
#   | | / _` | |______> > |  __  | |    | |   
#  _| || (_| |       / /  | |  | | |____| |   
# |_____\__,_|      /_/   |_|  |_|______|_|   

process.FitMuHLTEff.Efficiencies = cms.PSet(
    HLT_IsoMu24_eta2p1_avg      = HLT_IsoMu24_eta2p1_avg,
    #HLT_IsoMu24_eta2p1_pt       = HLT_IsoMu24_eta2p1_avg.clone( BinnedVariables = HLT_IsoMu24_eta2p1_avg.BinnedVariables.clone( pt       = pt24Bins.pt ) ),
    #HLT_IsoMu24_eta2p1_eta      = HLT_IsoMu24_eta2p1_avg.clone( BinnedVariables = HLT_IsoMu24_eta2p1_avg.BinnedVariables.clone( eta      = MuBins.eta ) ),
    #HLT_IsoMu24_eta2p1_reliso   = HLT_IsoMu24_eta2p1_avg.clone( BinnedVariables = HLT_IsoMu24_eta2p1_avg.BinnedVariables.clone( reliso   = CommonBins.reliso ) ),
    #HLT_IsoMu24_eta2p1_nvtx     = HLT_IsoMu24_eta2p1_avg.clone( BinnedVariables = HLT_IsoMu24_eta2p1_avg.BinnedVariables.clone( nvtx     = CommonBins.nvtx ) ),
    #HLT_IsoMu24_eta2p1_pileup   = HLT_IsoMu24_eta2p1_avg.clone( BinnedVariables = HLT_IsoMu24_eta2p1_avg.BinnedVariables.clone( pileup   = CommonBins.pileup, ) ),
    #HLT_IsoMu24_eta2p1_instlumi = HLT_IsoMu24_eta2p1_avg.clone( BinnedVariables = HLT_IsoMu24_eta2p1_avg.BinnedVariables.clone( instlumi = CommonBins.instlumi ) ),
    #HLT_IsoMu24_eta2p1_njet     = HLT_IsoMu24_eta2p1_avg.clone( BinnedVariables = HLT_IsoMu24_eta2p1_avg.BinnedVariables.clone( njet     = CommonBins.njet ) ),
    #HLT_IsoMu24_eta2p1_drjet    = HLT_IsoMu24_eta2p1_avg.clone( BinnedVariables = HLT_IsoMu24_eta2p1_avg.BinnedVariables.clone( drjet    = CommonBins.drjet ) ),
    #HLT_IsoMu24_eta2p1_ht       = HLT_IsoMu24_eta2p1_avg.clone( BinnedVariables = HLT_IsoMu24_eta2p1_avg.BinnedVariables.clone( ht       = CommonBins.ht ) ),
    #HLT_IsoMu24_eta2p1_met      = HLT_IsoMu24_eta2p1_avg.clone( BinnedVariables = HLT_IsoMu24_eta2p1_avg.BinnedVariables.clone( met      = CommonBins.met ) ),
    #HLT_IsoMu24_eta2p1_st       = HLT_IsoMu24_eta2p1_avg.clone( BinnedVariables = HLT_IsoMu24_eta2p1_avg.BinnedVariables.clone( st       = CommonBins.st ) ),
    #
    HLT_PFHT400_Mu5_PFMET45_avg      = HLT_PFHT400_Mu5_PFMET45_avg,
    #HLT_PFHT400_Mu5_PFMET45_pt       = HLT_PFHT400_Mu5_PFMET45_avg.clone( BinnedVariables = HLT_PFHT400_Mu5_PFMET45_avg.BinnedVariables.clone( pt       = pt5Bins.pt ) ),
    #HLT_PFHT400_Mu5_PFMET45_eta      = HLT_PFHT400_Mu5_PFMET45_avg.clone( BinnedVariables = HLT_PFHT400_Mu5_PFMET45_avg.BinnedVariables.clone( eta      = MuBins.eta ) ),
    #HLT_PFHT400_Mu5_PFMET45_reliso   = HLT_PFHT400_Mu5_PFMET45_avg.clone( BinnedVariables = HLT_PFHT400_Mu5_PFMET45_avg.BinnedVariables.clone( reliso   = CommonBins.reliso ) ),
    #HLT_PFHT400_Mu5_PFMET45_nvtx     = HLT_PFHT400_Mu5_PFMET45_avg.clone( BinnedVariables = HLT_PFHT400_Mu5_PFMET45_avg.BinnedVariables.clone( nvtx     = CommonBins.nvtx ) ),
    #HLT_PFHT400_Mu5_PFMET45_pileup   = HLT_PFHT400_Mu5_PFMET45_avg.clone( BinnedVariables = HLT_PFHT400_Mu5_PFMET45_avg.BinnedVariables.clone( pileup   = CommonBins.pileup, ) ),
    #HLT_PFHT400_Mu5_PFMET45_instlumi = HLT_PFHT400_Mu5_PFMET45_avg.clone( BinnedVariables = HLT_PFHT400_Mu5_PFMET45_avg.BinnedVariables.clone( instlumi = CommonBins.instlumi ) ),
    #HLT_PFHT400_Mu5_PFMET45_njet     = HLT_PFHT400_Mu5_PFMET45_avg.clone( BinnedVariables = HLT_PFHT400_Mu5_PFMET45_avg.BinnedVariables.clone( njet     = CommonBins.njet ) ),
    #HLT_PFHT400_Mu5_PFMET45_drjet    = HLT_PFHT400_Mu5_PFMET45_avg.clone( BinnedVariables = HLT_PFHT400_Mu5_PFMET45_avg.BinnedVariables.clone( drjet    = CommonBins.drjet ) ),
    #HLT_PFHT400_Mu5_PFMET45_ht       = HLT_PFHT400_Mu5_PFMET45_avg.clone( BinnedVariables = HLT_PFHT400_Mu5_PFMET45_avg.BinnedVariables.clone( ht       = CommonBins.ht ) ),
    #HLT_PFHT400_Mu5_PFMET45_met      = HLT_PFHT400_Mu5_PFMET45_avg.clone( BinnedVariables = HLT_PFHT400_Mu5_PFMET45_avg.BinnedVariables.clone( met      = CommonBins.met ) ),
    #HLT_PFHT400_Mu5_PFMET45_st       = HLT_PFHT400_Mu5_PFMET45_avg.clone( BinnedVariables = HLT_PFHT400_Mu5_PFMET45_avg.BinnedVariables.clone( st       = CommonBins.st ) ),
    #
    HLT_PFHT350_Mu15_PFMET45_avg      = HLT_PFHT350_Mu15_PFMET45_avg,
    #HLT_PFHT350_Mu15_PFMET45_pt       = HLT_PFHT350_Mu15_PFMET45_avg.clone( BinnedVariables = HLT_PFHT350_Mu15_PFMET45_avg.BinnedVariables.clone( pt       = pt15Bins.pt ) ),
    #HLT_PFHT350_Mu15_PFMET45_eta      = HLT_PFHT350_Mu15_PFMET45_avg.clone( BinnedVariables = HLT_PFHT350_Mu15_PFMET45_avg.BinnedVariables.clone( eta      = MuBins.eta ) ),
    #HLT_PFHT350_Mu15_PFMET45_reliso   = HLT_PFHT350_Mu15_PFMET45_avg.clone( BinnedVariables = HLT_PFHT350_Mu15_PFMET45_avg.BinnedVariables.clone( reliso   = CommonBins.reliso ) ),
    #HLT_PFHT350_Mu15_PFMET45_nvtx     = HLT_PFHT350_Mu15_PFMET45_avg.clone( BinnedVariables = HLT_PFHT350_Mu15_PFMET45_avg.BinnedVariables.clone( nvtx     = CommonBins.nvtx ) ),
    #HLT_PFHT350_Mu15_PFMET45_pileup   = HLT_PFHT350_Mu15_PFMET45_avg.clone( BinnedVariables = HLT_PFHT350_Mu15_PFMET45_avg.BinnedVariables.clone( pileup   = CommonBins.pileup, ) ),
    #HLT_PFHT350_Mu15_PFMET45_instlumi = HLT_PFHT350_Mu15_PFMET45_avg.clone( BinnedVariables = HLT_PFHT350_Mu15_PFMET45_avg.BinnedVariables.clone( instlumi = CommonBins.instlumi ) ),
    #HLT_PFHT350_Mu15_PFMET45_njet     = HLT_PFHT350_Mu15_PFMET45_avg.clone( BinnedVariables = HLT_PFHT350_Mu15_PFMET45_avg.BinnedVariables.clone( njet     = CommonBins.njet ) ),
    #HLT_PFHT350_Mu15_PFMET45_drjet    = HLT_PFHT350_Mu15_PFMET45_avg.clone( BinnedVariables = HLT_PFHT350_Mu15_PFMET45_avg.BinnedVariables.clone( drjet    = CommonBins.drjet ) ),
    #HLT_PFHT350_Mu15_PFMET45_ht       = HLT_PFHT350_Mu15_PFMET45_avg.clone( BinnedVariables = HLT_PFHT350_Mu15_PFMET45_avg.BinnedVariables.clone( ht       = CommonBins.ht ) ),
    #HLT_PFHT350_Mu15_PFMET45_met      = HLT_PFHT350_Mu15_PFMET45_avg.clone( BinnedVariables = HLT_PFHT350_Mu15_PFMET45_avg.BinnedVariables.clone( met      = CommonBins.met ) ),
    #HLT_PFHT350_Mu15_PFMET45_st       = HLT_PFHT350_Mu15_PFMET45_avg.clone( BinnedVariables = HLT_PFHT350_Mu15_PFMET45_avg.BinnedVariables.clone( st       = CommonBins.st ) ),
    #
    HLT_Mu40_PFHT350_avg      = HLT_Mu40_PFHT350_avg,
    #HLT_Mu40_PFHT350_pt       = HLT_Mu40_PFHT350_avg.clone( BinnedVariables = HLT_Mu40_PFHT350_avg.BinnedVariables.clone( pt       = pt40Bins.pt ) ),
    #HLT_Mu40_PFHT350_eta      = HLT_Mu40_PFHT350_avg.clone( BinnedVariables = HLT_Mu40_PFHT350_avg.BinnedVariables.clone( eta      = MuBins.eta ) ),
    #HLT_Mu40_PFHT350_reliso   = HLT_Mu40_PFHT350_avg.clone( BinnedVariables = HLT_Mu40_PFHT350_avg.BinnedVariables.clone( reliso   = CommonBins.reliso ) ),
    #HLT_Mu40_PFHT350_nvtx     = HLT_Mu40_PFHT350_avg.clone( BinnedVariables = HLT_Mu40_PFHT350_avg.BinnedVariables.clone( nvtx     = CommonBins.nvtx ) ),
    #HLT_Mu40_PFHT350_pileup   = HLT_Mu40_PFHT350_avg.clone( BinnedVariables = HLT_Mu40_PFHT350_avg.BinnedVariables.clone( pileup   = CommonBins.pileup, ) ),
    #HLT_Mu40_PFHT350_instlumi = HLT_Mu40_PFHT350_avg.clone( BinnedVariables = HLT_Mu40_PFHT350_avg.BinnedVariables.clone( instlumi = CommonBins.instlumi ) ),
    #HLT_Mu40_PFHT350_njet     = HLT_Mu40_PFHT350_avg.clone( BinnedVariables = HLT_Mu40_PFHT350_avg.BinnedVariables.clone( njet     = CommonBins.njet ) ),
    #HLT_Mu40_PFHT350_drjet    = HLT_Mu40_PFHT350_avg.clone( BinnedVariables = HLT_Mu40_PFHT350_avg.BinnedVariables.clone( drjet    = CommonBins.drjet ) ),
    #HLT_Mu40_PFHT350_ht       = HLT_Mu40_PFHT350_avg.clone( BinnedVariables = HLT_Mu40_PFHT350_avg.BinnedVariables.clone( ht       = CommonBins.ht ) ),
    #HLT_Mu40_PFHT350_met      = HLT_Mu40_PFHT350_avg.clone( BinnedVariables = HLT_Mu40_PFHT350_avg.BinnedVariables.clone( met      = CommonBins.met ) ),
    #HLT_Mu40_PFHT350_st       = HLT_Mu40_PFHT350_avg.clone( BinnedVariables = HLT_Mu40_PFHT350_avg.BinnedVariables.clone( st       = CommonBins.st ) ),
)

################################################################

#  _____        _    _     
# |  __ \      | |  | |    
# | |__) |__ _ | |_ | |__  
# |  ___// _` || __|| '_ \ 
# | |   | (_| || |_ | | | |
# |_|    \__,_| \__||_| |_|
                          
process.p = cms.Path()
if doEleGsf:
    process.p += process.FitEleGsfEff
if doEleId:
    process.p += process.FitEleIdEff
if doEleHLT:
    process.p += process.FitEleHLTEff
if doMuTrk:
    process.p += process.FitMuTrkEff
if doMuId:
    process.p += process.FitMuIdEff
if doMuHLT:
    process.p += process.FitMuHLTEff
