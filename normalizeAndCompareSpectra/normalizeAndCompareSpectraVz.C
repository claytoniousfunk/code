#include "/home/clayton/Analysis/code/myProcesses/hiforest/plugin/eventMap_hiForest.h"
#include <iostream>
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include <TF1.h>
#include "assert.h"
#include <fstream>
#include "TMath.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TLatex.h>
#include <TCut.h>
#include "TDatime.h"
#include <vector>
#include "TCanvas.h"
#include <array>

TH1D *pythia_vz, *data_vz, *r;

double func_temp(double *x, double *par){
    double vz = x[0];
    //int bin = r->FindBin(pt);
    //double fitval = par[0] + par[1]*r->GetBinContent(bin) + par[2]*TMath::Exp(-1*par[3]*r->GetBinContent(bin));
    //double fitval = par[0] + par[1]*pt + par[2]*TMath::Exp(par[3]*pt);
    double fitval = par[0] + par[1]*vz + par[2]*vz*vz + par[3]*vz*vz*vz ;
    return fitval;
}



void normalizeAndCompareSpectraVz(){

	TFile* f_pythia =TFile::Open("/home/clayton/Analysis/data/ppskim/pythia_skim_newWeights_muptcut_10.root"); 
	TFile* f_data =TFile::Open("/home/clayton/Analysis/code/makeData/V2/rootFiles/makeDataV2_pp_data_muptcut_10.root"); 


	
	
  // pt
  
	f_pythia->GetObject("h_vz",pythia_vz);
	f_data->GetObject("h_vz",data_vz);
    
    TF1 *func = new TF1("func",func_temp,-15.0,15.0,4);
  

  // normalize

  pythia_vz->Scale(1.0/pythia_vz->Integral());
  data_vz->Scale(1.0/data_vz->Integral());

  double vz_axis[31] = {-15.0,-14.0,-13.0,-12.0,-11.0,-10.0,-9.0,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0};
  
  
  
//////////////////////////////////////////////////////////////////////////////////    spectrum plot    /////////////////////////////////////////////////////////////////////
    

// Define the Canvas
    
      TCanvas *c2 = new TCanvas("c2", "canvas", 500, 500);
      
      c2->cd();
      
      //c2->SetLogy();
      TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.3, 1., 1.);
      //pad1->SetGridx();         // Vertical grid
      pad2->SetLeftMargin(0.15);
      pad2->SetLogy();
      pad2->Draw();             // Draw the upper pad: pad1
      pad2->cd();
      

      // jet pT 
    
    
           

      // normalize pythia_pt by bin width
  
    int Nbin1 = pythia_vz->GetSize();
    TH1D *pythia_vz_xnorm = (TH1D*) pythia_vz->Clone("pythia_vz_xnorm");

    for (int i=0;i<Nbin1;i++){

      double x1 = pythia_vz->GetBinWidth(i);
      double y1 = pythia_vz->GetBinContent(i);
      double z1 = pythia_vz->GetBinError(i);
      double r1 = 0.0;
      double e1 = 0.0;
      
      if(x1!=0){

        r1 = y1/x1;    
        e1 = z1/x1;
        pythia_vz_xnorm->SetBinContent(i,r1); 
        pythia_vz_xnorm->SetBinError(i,e1);
 
               
      }
      
    }
      
  // normalize data_pt by bin width
  
    int Nbin2 = data_vz->GetSize();
    TH1D* data_vz_xnorm = (TH1D*) data_vz->Clone("data_vz_xnorm");

    for (int j=0;j<Nbin2;j++){

      double x2 = data_vz->GetBinWidth(j);
      double y2 = data_vz->GetBinContent(j);
      double z2= data_vz->GetBinError(j);
      double r2= 0.0;
      double e2= 0.0;
      
      if(x2!=0){

        r2 = y2/x2;    
        e2 = z2/x2;
        data_vz_xnorm->SetBinContent(j,r2); 
        data_vz_xnorm->SetBinError(j,e2);
 
               
      }
      
    }

      double msize = 0.8;

      pythia_vz_xnorm->SetMarkerStyle(8);
      pythia_vz_xnorm->SetMarkerColor(kBlue);
      pythia_vz_xnorm->SetMarkerSize(msize);
      pythia_vz_xnorm->SetFillColorAlpha(kBlue,0.7);
      pythia_vz_xnorm->SetStats(0);
      

      data_vz_xnorm->SetMarkerStyle(8);
      data_vz_xnorm->SetMarkerColor(kRed);
      data_vz_xnorm->SetMarkerSize(msize);
      data_vz_xnorm->SetFillColorAlpha(kRed,0.7);
      data_vz_xnorm->SetStats(0);
      data_vz_xnorm->SetTitle("");   

   

      pythia_vz_xnorm->Draw("e2");
      data_vz_xnorm->Draw("e2 same");
      
      

      auto legend = new TLegend(0.91,0.4,0.99,0.6);
      legend->AddEntry(data_vz_xnorm,"Data","p");
      legend->AddEntry(pythia_vz_xnorm,"Pythia","p");
      legend->Draw();
      

     
      pythia_vz_xnorm->GetYaxis()->SetTitle("1/N^{jet} dN^{jet}/dv_{z}");
      pythia_vz_xnorm->GetYaxis()->SetTitleSize(0.04);


      pythia_vz_xnorm->SetTitle("");
      
	


///////////////////////////////////////////////////////////////////////////  ratio plot  ////////////////////////////////////////////////////////////////////////////


    c2->cd();
    TPad *pad1 = new TPad("pad1","pad1",0.,0.0,1.,0.3);
    //pad1->SetGridx();
    pad1->SetBottomMargin(0.3);
    pad1->SetLeftMargin(0.15);
    pad1->Draw();
    pad1->cd();
    r = (TH1D*)pythia_vz_xnorm->Clone("r");
    //TH1D *r2 = (TH1D*)hpr->Clone("r2");
    r->Divide(data_vz_xnorm);
    //r2->Divide(h1);
    r->SetMarkerStyle(8);
    r->SetMarkerColor(kBlack);
    r->SetMarkerSize(1.2);
    r->SetFillColorAlpha(kBlack,0.7);
    r->SetStats(0);
    r->GetXaxis()->SetTitle("jet v_{z} [mm]");
    r->GetXaxis()->SetTitleSize(0.09);
    r->GetYaxis()->SetTitle("Pythia / data");
    r->GetYaxis()->SetTitleSize(0.08);
    r->SetTitle("");
    r->Draw("e2");
    func->SetParNames("a_0","a_1","a_2","a_3","a_4");
    r->Fit("func","MR","N",-15.0,15.0);
    r->SetMinimum(0.8);
    r->SetMaximum(1.2);
    //gStyle->SetOptFit(1);
    //r->Fit("pol7");
    //r->Fit("expo");
    c2->SaveAs("/home/clayton/Analysis/code/normalizeAndCompareSpectra/figures/pythiaDataVzCompare.pdf");





    double a0 = func->GetParameter(0);
    double e_a0 = func->GetParError(0);
    double a1 = func->GetParameter(1);
    double e_a1 = func->GetParError(1);
    double a2 = func->GetParameter(2);
    double e_a2 = func->GetParError(2);
    double a3 = func->GetParameter(3);
    double e_a3 = func->GetParError(3);
    //double a4 = func->GetParameter(4);
    //double e_a4 = func->GetParError(4);
    TF1 *fit_fxn = new TF1("fit_fxn","[0]+[1]*x + [2]*x*x + [3]*x*x*x ",-15.0,15.0);
    fit_fxn->SetParameters(a0,a1,a2,a3);
    
    TCanvas *c3 = new TCanvas("c3","c3",500,500);
    c3->cd();
    TPad *pad3 = new TPad("pad3","pad3",0.0,0.3,1.0,1.0);
    TPad *pad4 = new TPad("pad4","pad4",0.0,0.0,1.0,0.3);
    pad4->SetLeftMargin(0.15);
    pad4->SetBottomMargin(0.3);
    pad4->Draw();
    pad3->SetLeftMargin(0.15);
    pad3->SetLogy();
    pad3->Draw();
    pad3->cd();

    
    //fit_fxn->Draw();

    TH1D* pythia_vz_xnorm_fitscaled = (TH1D*)pythia_vz_xnorm->Clone("pythia_vz_xnorm_fitscaled");
    int Nbin3 = pythia_vz_xnorm_fitscaled->GetSize();
    for(int k=0;k<Nbin3;k++){
        double pteval = pythia_vz_xnorm_fitscaled->GetXaxis()->GetBinCenter(k);
        double fitscale = fit_fxn->Eval(pteval);
        double oldval = pythia_vz_xnorm_fitscaled->GetBinContent(k);
        double olderror = pythia_vz_xnorm_fitscaled->GetBinError(k);
        double newval = oldval/fitscale;
        double newerror = newval*TMath::Sqrt(TMath::Power(olderror/oldval,2) + TMath::Power(e_a0/a0,2) + TMath::Power(e_a1/a1,2) + TMath::Power(e_a2/a2,2) + TMath::Power(e_a3/a3,2));
        pythia_vz_xnorm_fitscaled->SetBinContent(k,newval);
        pythia_vz_xnorm_fitscaled->SetBinError(k,newerror);
    }    
   
   pythia_vz_xnorm_fitscaled->SetMarkerColor(kGreen+2);
   pythia_vz_xnorm_fitscaled->SetFillColorAlpha(kGreen+2,0.7);
   pythia_vz_xnorm_fitscaled->SetMarkerStyle(24);
   data_vz_xnorm->GetYaxis()->SetTitle("1/N^{jet} dN^{jet}/dv_{z}");
   data_vz_xnorm->GetYaxis()->SetTitleSize(0.04);
   
   data_vz_xnorm->SetMarkerStyle(24);

   data_vz_xnorm->Draw("ep");
   pythia_vz_xnorm_fitscaled->Draw("ep same");

   pad4->cd();
   TH1D *r2 = (TH1D*) pythia_vz_xnorm_fitscaled->Clone("r2");
   r2->Divide(data_vz_xnorm);
   r2->SetMarkerColor(kBlack);
   
   r2->Draw("ep");
   r2->SetMinimum(0.8);
   r2->SetMaximum(1.2);
   r2->GetXaxis()->SetTitle("jet v_{z} [mm]");
   r2->GetYaxis()->SetTitle("Pythia / data");
   r2->GetXaxis()->SetTitleSize(0.09);
   r2->GetYaxis()->SetTitleSize(0.08);
   c3->SaveAs("/home/clayton/Analysis/code/normalizeAndCompareSpectra/figures/pythiaDataVzRescaleCompare.pdf");


    

    









}
