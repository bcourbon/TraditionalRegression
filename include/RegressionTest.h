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


    private:

    const std::vector<std::string> *varlistEB;
    const std::vector<std::string> *varlistEE;
    const GBRForest* forestEB;
    const GBRForest* forestEE;
    int numvarsEB;
    int numvarsEE;
    TTree *tree;



};

#endif
