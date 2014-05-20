#include "RegressionTest.h"
#include "GBRMaker.h"
#include "GBRForest.h"
#include "GBRTree.h"
#include "Utilities.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TMath.h"
#include "TPaveLabel.h"
#include "TStyle.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooBreitWigner.h"
#include "RooArgList.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"

#include <TSystem.h>
#include <TFile.h>

#include <iostream>
#include <sstream>
#include <algorithm>
#include <vector>

using namespace std;
using namespace RooFit;

/*****************************************************************/
RegressionTest::RegressionTest()
/*****************************************************************/
{
}



/*****************************************************************/
RegressionTest::~RegressionTest()
/*****************************************************************/
{
}



/*****************************************************************/
void RegressionTest::init(string RegFileName, std::string DataFileName, std::string DirName, std::string TreeName)
/*****************************************************************/
{

  TFile *file = new TFile(RegFileName.c_str(),"READ");
  
  varlist = (std::vector<std::string>*)file->Get("varlistEB");
  forest = (GBRForest*)file->Get("EBCorrection");
  
  numvars = varlist->size();
  numtrees = forest->Trees().size();
  
  printf("liste des variables \n");
  for (int ivar = 0; ivar<numvars; ++ivar) {
    printf("%i: %s\n",ivar, varlist->at(ivar).c_str());
  }
  printf("nombre de trees : %i \n",numtrees); 
  
  TFile *datafile = TFile::Open(DataFileName.c_str());
  TDirectory *dir = (TDirectory*)datafile->FindObjectAny(DirName.c_str());
  tree = (TTree*)dir->Get(TreeName.c_str());  
  

}

/*****************************************************************/
void RegressionTest::PlotResponse()
/*****************************************************************/
{

  double Correction; // Ecor/Eraw
  double RawBias;    // Eraw/Etrue
  double CorBias;    // Ecor/Etrue

  //vector of variables of interest
  std::vector<std::string> InputVars;
  for (int ivar = 0; ivar<numvars; ++ivar) {
  	InputVars.push_back(varlist->at(ivar));    
  }
  InputVars.push_back("genEnergy");  
  InputVars.push_back("genPt");  
  unsigned int numvars2=InputVars.size();
 
  //to read the variables in the ttree
  std::vector<TTreeFormula*> inputforms;
  for (std::vector<std::string>::const_iterator it = InputVars.begin(); it != InputVars.end(); ++it) {
	inputforms.push_back(new TTreeFormula(it->c_str(),it->c_str(),tree));
  }
  
  //contains ttree events
  float* Event;
  Event = new float[numvars2];


  TH1D *hCor=new TH1D("Ecor/Eraw","Ecor/Eraw",500,0.8,1.2);
  TH1D *hBiasCor=new TH1D("Ecor/Etrue","Ecor/Etrue",500,0.8,1.2);
  TH1D *hBiasRaw=new TH1D("Eraw/Etrue","Eraw/Etrue",500,0.8,1.2);
  
  TH2D *hCor_Eta=new TH2D("Ecor/Eraw vs Eta","Ecor/Eraw vs Eta",500,-3,3,500,0,2);
  TH2D *hCor_Phi=new TH2D("Ecor/Eraw vs Phi","Ecor/Eraw vs Phi",500,-3.14,3.14,500,0,2);
  TH2D *hCor_EtaWidth=new TH2D("Ecor/Eraw vs Eta Width","Ecor/Eraw vs Eta Width",500,0,0.03,500,0,2);
  TH2D *hCor_PhiWidth=new TH2D("Ecor/Eraw vs Phi Width","Ecor/Eraw vs Phi Width",500,0,0.1,500,0,2);
  TH2D *hCor_R9=new TH2D("Ecor/Eraw vs R9","Ecor/Eraw vs R9",500,0,1.2,500,0,2);
  TH2D *hCor_Nvtx=new TH2D("Ecor/Eraw vs nVtx","Ecor/Eraw vs nVtx",500,0,50,500,0,2);

  std::vector<TH1*> HistPt;
  for (int i=0;i<10;i++) HistPt.push_back(new TH1D(Form("Ecor/Etrue for %i < Pt <%i ",i*10,(i+1)*10),Form("Ecor/Etrue for %i < Pt <%i ",i*10,(i+1)*10),500,0.5,1.5));

  
  for (Long64_t ievt=0; ievt<tree->GetEntries();ievt++){
        if(ievt%10000==0) printf("%lld / %lld \n",ievt,tree->GetEntries());
  	tree->LoadTree(ievt);
    	for (unsigned int ivar=0; ivar<numvars2; ++ivar) {
        	Event[ivar] = inputforms[ivar]->EvalInstance();
    	}
        if (fabs(Event[1])<1.48){    	
	Correction=forest->GetResponse(Event);
        double Eraw=Event[32];
	double Etrue=Event[33];
        RawBias=Eraw/Etrue;
	CorBias=RawBias*Correction;
    	hCor->Fill(Correction);
    	hBiasRaw->Fill(RawBias);
    	hBiasCor->Fill(CorBias);
	hCor_Eta->Fill(Event[1],Correction);
	hCor_Phi->Fill(Event[2],Correction);
	hCor_EtaWidth->Fill(Event[3],Correction);
	hCor_PhiWidth->Fill(Event[4],Correction);
	hCor_R9->Fill(Event[5],Correction);
	hCor_Nvtx->Fill(Event[0],Correction);
	double Pt=Event[34];
	for(int i=0;i<10;i++){
	if (Pt>i*10 && Pt<(i+1)*10) HistPt[i]->Fill(CorBias);
	}
  }
  }

  //Compute bias and resolution in function of pt bu fitting the distributions by a double Crystal Ball
  
  double pt[10];
  double mean[10];
  double sigeff[10];

  //create fit function and datahist
  RooRealVar *bias=new RooRealVar("Ecor/Etrue","Ecor/Etrue",0.5,1.5);
  const RooArgList *var=new RooArgList(*bias,"");
  RooRealVar *mu=new RooRealVar("mu","mu",1,0.5,1.5);
  RooRealVar *sig=new RooRealVar("sig","sig",0.05,0.001,1);
  RooRealVar *alpha1=new RooRealVar("alpha1","alpha1",0.2,0.,1.);
  RooRealVar *alpha2=new RooRealVar("alpha2","alpha2",-0.2,-1.,0.);
  RooRealVar *n1=new RooRealVar("n1","n2",2,1,10);
  RooRealVar *n2=new RooRealVar("n2","n2",2,1,10);
  RooRealVar *frac=new RooRealVar("frac","frac",0.5);
  RooCBShape *CB1=new RooCBShape("cb1","cb1",*bias,*mu,*sig,*alpha1,*n1);
  RooCBShape *CB2=new RooCBShape("cb2","cb2",*bias,*mu,*sig,*alpha2,*n2);
  RooArgList *cbs = new RooArgList();
  RooArgList *coeffs = new RooArgList();
  cbs->add(*CB1);
  cbs->add(*CB2);
  coeffs->add(*frac);
  RooAddPdf *DoubleCB=new RooAddPdf("Double CB","Double CB",*cbs,*coeffs);

  TCanvas *can7=new TCanvas("c7","c7",1200,600);
  can7->Divide(5,2); 

  for (int i=0;i<10;i++){
  pt[i]=5.+i*10.;
  RooDataHist *HistPtclone=new RooDataHist("dh","dh",*var,HistPt[i]);
  DoubleCB->fitTo(*HistPtclone); 
  RooArgSet* paramSet=DoubleCB->getParameters(HistPtclone) ;
  mean[i]=paramSet->getRealValue("mu"); 
  sigeff[i]=paramSet->getRealValue("sig"); 
  //sigeff[i]=effSigma(HistPt[i]);
  RooPlot* frame = bias->frame(Title(Form("Ecor/Etrue for %i<Pt<%i",i*10,(i+1)*10)));
  HistPtclone->plotOn(frame) ; 
  DoubleCB->plotOn(frame);
  can7->cd(i+1);
  frame->Draw();
  delete HistPtclone;
  }

  can7->SaveAs("plots/PtDistributions.png");
    
  //plot results
  TCanvas *can = new TCanvas;
  hCor->Draw();
  hCor->GetXaxis()->SetTitle("Ecor/Eraw");
  hCor->SetLineWidth(2);
  can->SaveAs("plots/Correction.png");

  TCanvas *can2 = new TCanvas;
  hBiasCor->Draw();
  hBiasRaw->Draw("same");
  hBiasCor->SetLineColor(kRed);
  hBiasRaw->SetLineColor(kBlue);
  hBiasCor->SetLineWidth(2);
  hBiasRaw->SetLineWidth(2);
  TLegend *leg = new TLegend(0.12,0.78,0.32,0.88);
  leg->AddEntry(hBiasRaw,"Eraw/Etrue","l");
  leg->AddEntry(hBiasCor,"Ecor/Etrue","l");
  leg->Draw();
  hBiasCor->SetTitle("E/Etrue");
  hBiasCor->GetXaxis()->SetTitle("E/Etrue");
  gStyle->SetOptStat(0);
  TPaveLabel *OldPerf = new TPaveLabel(0.63,0.88,0.88,0.78, Form("Peak = %.3f , RMS = %.3f ",hBiasRaw->GetXaxis()->GetBinCenter(hBiasRaw->GetMaximumBin()),hBiasRaw->GetRMS()),"NDC");
  TPaveLabel *NewPerf = new TPaveLabel(0.63,0.76,0.88,0.66, Form("Peak = %.3f , RMS = %.3f ",hBiasCor->GetXaxis()->GetBinCenter(hBiasCor->GetMaximumBin()),hBiasCor->GetRMS()),"NDC");
  OldPerf->SetBorderSize(1);
  OldPerf->SetTextSize(0.25);
  OldPerf->SetTextColor(kBlue);
  NewPerf->SetBorderSize(1);
  NewPerf->SetTextSize(0.25);
  NewPerf->SetTextColor(kRed);
  OldPerf->Draw();
  NewPerf->Draw();
  can2->SaveAs("plots/Performance.png");

  TCanvas *can3 = new TCanvas("","",1200,700);
  can3->Divide(3,2);
  can3->cd(1);
  TProfile *pCor_Eta= new TProfile();
  pCor_Eta=hCor_Eta->ProfileX("",1,-1,"s");
  pCor_Eta->SetMinimum(0.9);
  pCor_Eta->SetMaximum(1.1);
  pCor_Eta->Draw();
  pCor_Eta->GetXaxis()->SetTitle("Eta");
  pCor_Eta->GetYaxis()->SetTitle("Ecor/Eraw");
  pCor_Eta->SetTitle("Profile of Ecor/Eraw vs Eta");
  can3->cd(2);
  TProfile *pCor_Phi= new TProfile();
  pCor_Phi=hCor_Phi->ProfileX("",1,-1,"s");
  pCor_Phi->SetMinimum(0.9);
  pCor_Phi->SetMaximum(1.1);
  pCor_Phi->Draw();
  pCor_Phi->GetXaxis()->SetTitle("Phi");
  pCor_Phi->GetYaxis()->SetTitle("Ecor/Eraw");
  pCor_Phi->SetTitle("Profile of Ecor/Eraw vs Phi");
  can3->cd(3);
  TProfile *pCor_EtaWidth= new TProfile();
  pCor_EtaWidth=hCor_EtaWidth->ProfileX("",1,-1,"s");
  pCor_EtaWidth->SetMinimum(0.9);
  pCor_EtaWidth->SetMaximum(1.2);
  pCor_EtaWidth->Draw();
  pCor_EtaWidth->GetXaxis()->SetTitle("EtaWidth");
  pCor_EtaWidth->GetYaxis()->SetTitle("Ecor/Eraw");
  pCor_EtaWidth->SetTitle("Profile of Ecor/Eraw vs EtaWidth");
  can3->cd(4);
  TProfile *pCor_PhiWidth= new TProfile();
  pCor_PhiWidth=hCor_PhiWidth->ProfileX("",1,-1,"s");
  pCor_PhiWidth->SetMinimum(0.9);
  pCor_PhiWidth->SetMaximum(1.2);
  pCor_PhiWidth->Draw();
  pCor_PhiWidth->GetXaxis()->SetTitle("PhiWidth");
  pCor_PhiWidth->GetYaxis()->SetTitle("Ecor/Eraw");
  pCor_PhiWidth->SetTitle("Profile of Ecor/Eraw vs PhiWidth");
  can3->cd(5);
  TProfile *pCor_R9= new TProfile();
  pCor_R9=hCor_R9->ProfileX("",1,-1,"s");
  pCor_R9->SetMinimum(0.9);
  pCor_R9->SetMaximum(1.2);
  pCor_R9->Draw();
  pCor_R9->GetXaxis()->SetTitle("R9");
  pCor_R9->GetYaxis()->SetTitle("Ecor/Eraw");
  pCor_R9->SetTitle("Profile of Ecor/Eraw vs R9");
  can3->cd(6);
  TProfile *pCor_Nvtx= new TProfile();
  pCor_Nvtx=hCor_Nvtx->ProfileX("",1,-1,"s");
  pCor_Nvtx->SetMinimum(0.9);
  pCor_Nvtx->SetMaximum(1.1);
  pCor_Nvtx->Draw();
  pCor_Nvtx->GetXaxis()->SetTitle("Nvtx");
  pCor_Nvtx->GetYaxis()->SetTitle("Ecor/Eraw");
  pCor_Nvtx->SetTitle("Profile of Ecor/Eraw vs Nvtx");
  can3->SaveAs("plots/Correction_Variables.png");

  TCanvas *can4=new TCanvas("c4","c4",1100,500);
  can4->Divide(2,0);
  can4->cd(1);
  TGraph *graphCor_Pt = new TGraph(10,pt,mean);
  graphCor_Pt->Draw("AL*");
  graphCor_Pt->SetMarkerStyle(8);
  graphCor_Pt->GetXaxis()->SetTitle("pt");
  graphCor_Pt->GetYaxis()->SetTitle("Ecor/Etrue mean");
  graphCor_Pt->SetTitle("Ecor/Etrue peak value vs Pt");  
  can4->cd(2);
  TGraph *graphRes_Pt = new TGraph(10,pt,sigeff);
  graphRes_Pt->Draw("AL*");
  graphRes_Pt->SetMarkerStyle(8);
  graphRes_Pt->GetXaxis()->SetTitle("pt");
  graphRes_Pt->GetYaxis()->SetTitle("effective sigma");
  graphRes_Pt->SetTitle("Ecor/Etrue effective sigma vs Pt");  
  can4->SaveAs("plots/Correction_Pt.png");

  TCanvas *can5=new TCanvas("c5","c5",1200,600);
  can5->Divide(5,2);
  for (int i=0;i<10;i++){
  can5->cd(i+1);
  HistPt[i]->Draw();
  }
  gStyle->SetOptStat(1);
  can5->SaveAs("plots/CorPtRaw.png");

}
