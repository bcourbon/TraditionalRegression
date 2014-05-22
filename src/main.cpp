/**
 *  @file  main.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    11/12/2012
 *
 *  @internal
 *     Created :  11/12/2012
 * Last update :  11/12/2012 10:56:02 AM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#include <vector>
#include <iostream>
#include <string>
#include "TFile.h"
#include "RegressionManager.h"
#include "RegressionTest.h"


int main(int argc, char** argv)
{
    if(argc!=2)
    {
        std::cout<<"Usage: regression.exe configurationFile\n";
        return 1;
    }

    std::string parameterFile(argv[1]);
    RegressionManager manager;
    bool status = true;
    status = manager.init(parameterFile);
    if(status)
    {
        manager.makeRegression();
    }

    if(!status)
        std::cout<<"FATAL: A fatal error occured - QUIT -\n";
    else
        std::cout<<"- Finish - All good -\n";

    RegressionTest test;
    test.init("GBR_Clustering_70pre11_Photons_EG_results.root","~/eos/cms/store/group/phys_egamma/lgray/PhotonRegressionTrees_23012014_v2-runPhotonRegressionTrees_cfg/runPhotonRegressionTrees_cfg-step4_RECO_EI-FA47A89F-906A-E211-9ABF-003048FFCBB0.root","egSCTree","SuperClusterTree");
    test.PlotResponse();
       
    return status;
}
