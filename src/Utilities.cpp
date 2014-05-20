/**
 *  @file  Utilities.cxx
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    03/27/2010
 *
 *  @internal
 *     Created :  03/27/2010
 * Last update :  03/27/2010 01:00:52 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */


#include "Utilities.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"

//--- STL

using namespace std;




///*****************************************************************/
//void tokenize(const string& str,
//        vector<string>& tokens,
//        const string& delimiters)
///*****************************************************************/
//{
//    // Skip delimiters at beginning.
//    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
//    // Find first "non-delimiter".
//    string::size_type pos     = str.find_first_of(delimiters, lastPos);
//
//    while (string::npos != pos || string::npos != lastPos)
//    {
//        // Found a token, add it to the vector.
//        tokens.push_back(str.substr(lastPos, pos - lastPos));
//        // Skip delimiters.  Note the "not_of"
//        lastPos = str.find_first_not_of(delimiters, pos);
//        // Find next "non-delimiter"
//        pos = str.find_first_of(delimiters, lastPos);
//    }
//}

/*****************************************************************/
void tokenize(const string& str,
        vector<string>& tokens,
        const string& delimiter)
/*****************************************************************/
{
    string::size_type length = delimiter.size();
    string::size_type lastPos = 0;
    string::size_type pos     = str.find(delimiter, 0);


    while (string::npos != pos)
    {
        // Found a token, add it to the vector.
        if(str.substr(lastPos, pos - lastPos).size()>0)
            tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = pos + length;
        // Find next "non-delimiter"
        pos = str.find(delimiter, lastPos);
    }
    if(str.substr(lastPos).size()>0)
        tokens.push_back(str.substr(lastPos));
}



/*****************************************************************/
string intToString(int n)
/*****************************************************************/
{
    ostringstream oss;
    oss << n;
    return oss.str();
}



/*****************************************************************/
void findAndReplace(string& sInput, string sFind, string sReplace )
/*****************************************************************/
{
    size_t itPos = 0; 
    size_t itFindLen = sFind.length();
    size_t itReplaceLen = sReplace.length();

    if( itFindLen == 0 )
        return;

    while( (itPos = sInput.find( sFind, itPos )) != std::string::npos )
    {
        sInput.replace( itPos, itFindLen, sReplace );
        itPos += itReplaceLen;
    }

}


/*****************************************************************/
void strip(std::string& sInput)
/*****************************************************************/
{
	//-- removing blanks at the beginning and at the end
	while(*sInput.begin()==' ' || *sInput.begin()=='\t') sInput.erase(sInput.begin());
	while(*(sInput.end()-1)==' ' || *(sInput.end()-1)=='\t') sInput.erase(sInput.end()-1);
}

/*****************************************************************/
//effsigma function from Chris
double effSigma(TH1 * hist)
/*****************************************************************/
{

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
    return 0.;
  }
  
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
    return 0.;
  }
  //Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
//   if(total < 100.) {
//     cout << "effsigma: Too few entries " << total << endl;
//     return 0.;
//   }
  Int_t ierr=0;
  Int_t ismin=999;
  
  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }   
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;
  
  return widmin;
  
}


/*****************************************************************/
double DoubleCrystalBall(double *x, double *par)
/*****************************************************************/
{
double mean=par[0];
double sigma=par[1];
double alpha1=par[2];
double n1=par[3];
double alpha2=par[4];
double n2=par[5];

double t = (x[0]-mean)/sigma;
double val = -99.;
if(t>-alpha1 && t<alpha2){
     val = exp(-0.5*t*t);
}else if(t<=-alpha1){
     double alpha1invn1 = alpha1/n1;
     val = exp(-0.5*alpha1*alpha1)*pow(1. - alpha1invn1*(alpha1+t), -n1);

}else if(t>=alpha2){
     double alpha2invn2 = alpha2/n2;
     val = exp(-0.5*alpha2*alpha2)*pow(1. - alpha2invn2*(alpha2-t), -n2);     
}  

return val;

}

/*****************************************************************/
double GetMeanAfterFit(TH1 *histo)
/*****************************************************************/
{

TF1 *FitFunction = new TF1("FitFunction",DoubleCrystalBall,0,3,6);
FitFunction->SetParameter(0,1);
FitFunction->SetParLimits(0,0.5,1.5);
FitFunction->SetParameter(1,0.1);
FitFunction->SetParLimits(1,0.01,1);
FitFunction->SetParameter(2,2);
FitFunction->SetParLimits(2,2,2);
FitFunction->SetParameter(3,3);
FitFunction->SetParLimits(3,1.01,10);
FitFunction->SetParameter(4,1);
FitFunction->SetParLimits(4,1,1);
FitFunction->SetParameter(5,3);
FitFunction->SetParLimits(5,1.01,10);

histo->Fit("FitFunction");
return FitFunction->GetParameter(0);

}
