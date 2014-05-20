from RegressionParametricJobLauncher import RegressionParametricJobLauncher
import glob

info = """
Based on 70pre11

SC regressions with weights
* PF-Mustache + Refined SC electron + EG
* Barrel + Endcap

RandomSeed=nEvent+123456*el_index
EventWeight=
min(1,exp(-(genPt-50)/50))

Correction + error

"""


inputFiles = []
inputFiles.extend(glob.glob("/grid_mnt/data1__data.polcms/cms/sauvan/Regression/Ntuples/SingleGammaPt300ExtRelVal610/CMSSW_6_1_0-PU_START61_V8-v1/vJan23_700pre11_FromLindsey23012014/runPhotonRegressionTrees_cfg-step4_RECO_EI*.root"))




commonVariables = [
    "nVtx",
    #"scRawEnergy",
    "scEta",
    "scPhi",
    "scEtaWidth",
    "scPhiWidth",
    "scSeedR9",
    "scSeedRawEnergy/scRawEnergy",
    "scSeedEmax/scRawEnergy",
    "scSeedE2nd/scRawEnergy",
    "scSeedLeftRightAsym",
    "scSeedTopBottomAsym",
    "scSeedSigmaIetaIeta",
    "scSeedSigmaIetaIphi",
    "scSeedSigmaIphiIphi",
    "N_ECALClusters",
    "clusterMaxDR",
    "clusterMaxDRDPhi",
    "clusterMaxDRDEta",
    "clusterMaxDRRawEnergy/scRawEnergy",
    "clusterRawEnergy[0]/scRawEnergy",
    "clusterRawEnergy[1]/scRawEnergy",
    "clusterRawEnergy[2]/scRawEnergy",
    "clusterDPhiToSeed[0]",
    "clusterDPhiToSeed[1]",
    "clusterDPhiToSeed[2]",
    "clusterDEtaToSeed[0]",
    "clusterDEtaToSeed[1]",
    "clusterDEtaToSeed[2]",
]


barrelVariables = [
    "scSeedCryEta",
    "scSeedCryPhi",
    "scSeedCryIeta",
    "scSeedCryIphi"
]

endcapVariables = [
    "scPreshowerEnergy/scRawEnergy"
]




commonVariablesEB = commonVariables + barrelVariables
commonVariablesEE = commonVariables + endcapVariables



batch = RegressionParametricJobLauncher()
batch.baseDir = "/home/llr/cms/sauvan/RegressionTraining/"
batch.baseName = "GBR_Clustering_70pre11_Photons"
batch.inputFiles = inputFiles
batch.outputDirectory = "/home/llr/cms/sauvan/DATA/sauvan/RegressionResults/StudyClustering/70pre11/Test/GBR_Clustering_Photons_MoreVariables_PtCut50_PtSlope50_Sig5/"
batch.commonOptions = ["MinEvents=200","Shrinkage=0.1","NTrees=1000","MinSignificance=5.0","RandomSeed=eventNumber+123456*scIndex","EventWeight=min(1,exp(-(genPt-50)/50))"]
batch.commonVariablesEB = commonVariablesEB
batch.commonVariablesEE = commonVariablesEE
batch.commonCutsEB      = ["scIsEB"]
batch.commonCutsEE      = ["!scIsEB"]
batch.commonCuts        = ["(eventNumber%2==0)&&(isMatched==1)"]
batch.commonCutsError   = ["(eventNumber%2!=0)&&(((eventNumber-1)/2)%4==3)&&(isMatched==1)"]
batch.doErrors = True
batch.doCombine = False

batch.target = "genEnergy/(scRawEnergy+scPreshowerEnergy)"
targets = {
    "EG":"genEnergy/(scRawEnergy+scPreshowerEnergy)",
    #"PFBox":"genEnergy/scCalibratedEnergy",
    "PFMustache":"genEnergy/scCalibratedEnergy",
    "GedPhoton":"genEnergy/scCalibratedEnergy",
}
specificVariables = {
    "EG":["scRawEnergy"],
    #"PFBox":["scCalibratedEnergy"],
    "PFMustache":["scCalibratedEnergy"],
    "GedPhoton":["scCalibratedEnergy"],
}



clusterings = [
    #("PFBox","regularBoxSCTree/SuperClusterTree"),
    ("PFMustache","mustacheSCTree/SuperClusterTree"),
    ("EG","egSCTree/SuperClusterTree"),
    ("GedPhoton","gedPhotonTree/RegressionTree"),
]

for tagClustering,treeClustering in clusterings:
    # define regressions
    batch.addRegression(tagClustering)
    batch.setInputTree(tagClustering, treeClustering)
    batch.setTarget(tagClustering, targets[tagClustering])
    batch.addVariablesEB(tagClustering, specificVariables[tagClustering])
    batch.addVariablesEE(tagClustering, specificVariables[tagClustering])

#batch.simulate = True

batch.info = info

