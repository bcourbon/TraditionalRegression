#ifndef REGRESSIONTEST_H
#define REGRESSIONTEST_H


#include <string>
#include <vector>
#include "GBRForest.h"
#include "TTree.h"

class RegressionTest
{
    public:
        RegressionTest();
        ~RegressionTest();

        void init(std::string RegFileName,std::string DataFileName,std::string DirName,std::string TreeName);
        void PlotResponse();
    const std::vector<std::string> *varlist;

    private:

    const GBRForest* forest;
    TTree *tree;
    int numvars;
    int numtrees;

};

#endif
