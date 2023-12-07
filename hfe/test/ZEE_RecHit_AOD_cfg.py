import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )


process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

options = VarParsing.VarParsing('standard')

options.register('inputFile',
        "~/",
        VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.string,
        "File containing a list of the EXACT location of the output file  (default = ~/)"
        )


options.parseArguments()
options.inputFile = 'root://eoscms//' + options.inputFile
print(options.inputFile)
process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
          'root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18RECO/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/260000/B29C79BD-90F3-784B-B80C-B2978529123F.root'
#          'root://cms-xrd-global.cern.ch///store/data/Run2022E/EGamma/AOD/27Jun2023-v1/25330005/d173e4c1-07a5-49d5-b2c4-605990493520.root'
#	  'root://cms-xrd-global.cern.ch///store/mc/RunIIWinter19PFCalibDR/DoubleElectron_FlatPt-1To300/AODSIM/2018ConditionsFlatPU0to70ECALGT_105X_upgrade2018_realistic_IdealEcalIC_v4-v1/40001/C6611D6B-CD06-C346-9828-BFA9EACE213E.root'
#                options.inputFile
                )
                            )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_v10', '')

'''
from CondCore.DBCommon.CondDBSetup_cfi import *

process.GlobalTag = cms.ESSource("PoolDBESSource",
                               CondDBSetup,
                               connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
                               globaltag = cms.string('106X_upgrade2018_realistic_v11_L1v1'),
                               toGet = cms.VPSet(

                                        cms.PSet(record = cms.string("GBRDWrapperRcd"),
                 tag = cms.string("pfscecal_EBCorrection_offline_v2_2018UL"),
                 label = cms.untracked.string("pfscecal_EBCorrection_offline_v2"),
                 connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
                 ),

        cms.PSet(record = cms.string("GBRDWrapperRcd"),
                 tag = cms.string("pfscecal_EBUncertainty_offline_v2_2018UL"),
                 label = cms.untracked.string("pfscecal_EBUncertainty_offline_v2"),
                 connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
                 ),


        cms.PSet(record = cms.string("GBRDWrapperRcd"),
                 tag = cms.string("pfscecal_EECorrection_offline_v2_2018UL"),
                 label = cms.untracked.string("pfscecal_EECorrection_offline_v2"),
                 connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
                 ),


        cms.PSet(record = cms.string("GBRDWrapperRcd"),
                 tag = cms.string("pfscecal_EEUncertainty_offline_v2_2018UL"),
                 label = cms.untracked.string("pfscecal_EEUncertainty_offline_v2"),
                 connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
                 ),

                                )
                              )
'''


process.nTuplelize = cms.EDAnalyzer('ZEE_RecHit_NTuplizer',
        electrons = cms.InputTag("gedGsfElectrons")
	)


process.TFileService = cms.Service("TFileService",
#     fileName = cms.string("nTuple_Data.root"),
#     fileName = cms.string("Tree_Gamma_ABCD.root"),
     fileName = cms.string("nTuple_MC.root"),
      closeFileFast = cms.untracked.bool(True)
  )


process.p = cms.Path(process.nTuplelize)

