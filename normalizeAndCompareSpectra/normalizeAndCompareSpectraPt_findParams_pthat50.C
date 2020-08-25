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
#include "TLegend.h"

TH1D *pythia_pt, *data_pt, *r;

double func_temp(double *x, double *par){
    double pt = x[0];
    //int bin = r->FindBin(pt);
    //double fitval = par[0] + par[1]*r->GetBinContent(bin) + par[2]*TMath::Exp(-1*par[3]*r->GetBinContent(bin));
    //double fitval = par[0] + par[1]*pt + par[2]*TMath::Exp(par[3]*pt);
    double fitval = par[0]*TMath::Exp(par[1]*(pt-100)) + par[2]*TMath::Exp(par[3]*(pt-100));
    //double fitval = par[0] + par[1]*pt + par[2]*pt*pt + par[3]*pt*pt*pt;
    return fitval;
}



void normalizeAndCompareSpectraPt_findParams_pthat50(){

	//TFile* f_pythia =TFile::Open("/home/clayton/Analysis/data/ppskim/clayton_pp_mc_skim_muptcut_10.root"); 
	TFile* f_data =TFile::Open("/home/clayton/Analysis/code/makeData/V2/rootFiles/makeDataV2_pp_data_corrpt_muptcut_10.root"); 
	//TFile *f_pythia = TFile::Open("/home/clayton/Analysis/data/ppskim/pythia_skim_newWeights_pthatcut_50_muptcut_10.root");
	TFile *f_pythia = TFile::Open("/home/clayton/Analysis/data/ppskim/pythia_skim_newWeights_pthatcut_50_muptcut_10_16Aug2020.root");

	
	
  // pt
  
	f_pythia->GetObject("h_jetpt_raw",pythia_pt);
	f_data->GetObject("h_jetpt",data_pt);
    
    TF1 *func = new TF1("func",func_temp,100.0,500.0,4);
  

  // normalize

  pythia_pt->Scale(1.0/pythia_pt->Integral());
  data_pt->Scale(1.0/data_pt->Integral());

  double pt_axis[16] = {50.0,60.0,70.0,80.0,90.0,100.0,120.0,140.0,160.0,180.0,200.0,250.0,300.0,350.0,400.0,500.0};
  cout << pt_axis[1] << endl;
  

//////////////////////////////////////////////////////////////////////////////////    spectrum plot    /////////////////////////////////////////////////////////////////////
    

// Define the Canvas
    
      TCanvas *c2 = new TCanvas("c2", "Unweighted pt vs data", 500, 500);
      
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
  
    int Nbin1 = pythia_pt->GetSize();
    TH1D *pythia_pt_xnorm = (TH1D*) pythia_pt->Clone("pythia_pt_xnorm");

    for (int i=0;i<Nbin1;i++){

      double x1 = pythia_pt_xnorm->GetBinWidth(i);
      double y1 = pythia_pt_xnorm->GetBinContent(i);
      double z1 = pythia_pt_xnorm->GetBinError(i);
      double r1 = 0.0;
      double e1 = 0.0;
      
      if(x1!=0){

        r1 = y1/x1;    
        e1 = z1/x1;
        pythia_pt_xnorm->SetBinContent(i,r1); 
        pythia_pt_xnorm->SetBinError(i,e1);
 
               
      }
      
    }
      
  // normalize data_pt by bin width
  
    int Nbin2 = data_pt->GetSize();
    TH1D* data_pt_xnorm = (TH1D*) data_pt->Clone("data_pt_xnorm");

    for (int j=0;j<Nbin2;j++){

      double x2 = data_pt_xnorm->GetBinWidth(j);
      double y2 = data_pt_xnorm->GetBinContent(j);
      double z2= data_pt_xnorm->GetBinError(j);
      double r2= 0.0;
      double e2= 0.0;
      
      if(x2!=0){

        r2 = y2/x2;    
        e2 = z2/x2;
        data_pt_xnorm->SetBinContent(j,r2); 
        data_pt_xnorm->SetBinError(j,e2);
 
               
      }
      
    }
	TH1D *pythia_pt_xnorm_rebin = (TH1D*) pythia_pt_xnorm->Rebin(15,"pythia_pt_xnorm_rebin",pt_axis);
	TH1D *data_pt_xnorm_rebin = (TH1D*) data_pt_xnorm->Rebin(15,"data_pt_xnorm_rebin",pt_axis);
    
    
    // normalize pythia_pt_xnorm_rebin by bin width
    
    int Nbin3 = pythia_pt_xnorm_rebin->GetSize();
    TH1D *pythia_pt_xnorm_rebin_xnorm = (TH1D*) pythia_pt_xnorm_rebin->Clone("pythia_pt_xnorm_rebin_xnorm");

    for (int k=0;k<Nbin3;k++){

      double x3 = pythia_pt_xnorm_rebin_xnorm->GetBinWidth(k);
      double y3 = pythia_pt_xnorm_rebin_xnorm->GetBinContent(k);
      double z3 = pythia_pt_xnorm_rebin_xnorm->GetBinError(k);
      double r3 = 0.0;
      double e3 = 0.0;
      
      if(x3!=0){

        r3 = y3/x3;    
        e3 = z3/x3;
        pythia_pt_xnorm_rebin_xnorm->SetBinContent(k,r3); 
        pythia_pt_xnorm_rebin_xnorm->SetBinError(k,e3);
 
               
      }
      
    }
      
// normalize pythia_pt_xnorm_rebin by bin width
    
    int Nbin4 = data_pt_xnorm_rebin->GetSize();
    TH1D *data_pt_xnorm_rebin_xnorm = (TH1D*) data_pt_xnorm_rebin->Clone("data_pt_xnorm_rebin_xnorm");

    for (int l=0;l<Nbin4;l++){

      double x4 = data_pt_xnorm_rebin_xnorm->GetBinWidth(l);
      double y4 = data_pt_xnorm_rebin_xnorm->GetBinContent(l);
      double z4 = data_pt_xnorm_rebin_xnorm->GetBinError(l);
      double r4 = 0.0;
      double e4 = 0.0;
      
      if(x4!=0){

        r4 = y4/x4;    
        e4 = z4/x4;
        data_pt_xnorm_rebin_xnorm->SetBinContent(l,r4); 
        data_pt_xnorm_rebin_xnorm->SetBinError(l,e4);
 
               
      }
      
    }

      double msize = 1.0;

      pythia_pt_xnorm_rebin_xnorm->SetMarkerStyle(24);
      pythia_pt_xnorm_rebin_xnorm->SetMarkerColor(kBlue);
      pythia_pt_xnorm_rebin_xnorm->SetMarkerSize(msize);
      pythia_pt_xnorm_rebin_xnorm->SetFillColorAlpha(kBlue,0.7);
      pythia_pt_xnorm_rebin_xnorm->SetStats(0);
      

      data_pt_xnorm_rebin_xnorm->SetMarkerStyle(24);
      data_pt_xnorm_rebin_xnorm->SetMarkerColor(kRed);
      data_pt_xnorm_rebin_xnorm->SetMarkerSize(msize);
      data_pt_xnorm_rebin_xnorm->SetFillColorAlpha(kRed,0.7);
      data_pt_xnorm_rebin_xnorm->SetStats(0);
      data_pt_xnorm_rebin_xnorm->SetTitle("");   

   

      pythia_pt_xnorm_rebin_xnorm->Draw("e2");
      data_pt_xnorm_rebin_xnorm->Draw("e2 same");
      
      

      auto legend = new TLegend(0.5,0.65,0.89,0.82);
      legend->AddEntry(data_pt_xnorm_rebin_xnorm,"Data","p");
      legend->AddEntry(pythia_pt_xnorm_rebin_xnorm,"Pythia","p");
      legend->SetBorderSize(0);
      legend->Draw();
      

     
      pythia_pt_xnorm_rebin_xnorm->GetYaxis()->SetTitle("1/N^{jet} dN^{jet}/dp_{T}");
      pythia_pt_xnorm_rebin_xnorm->GetYaxis()->SetTitleSize(0.06);


      pythia_pt_xnorm_rebin_xnorm->SetTitle("");
      
	


///////////////////////////////////////////////////////////////////////////  ratio plot  ////////////////////////////////////////////////////////////////////////////


    c2->cd();
    TPad *pad1 = new TPad("pad1","pad1",0.,0.0,1.,0.3);
    //pad1->SetGridx();
    pad1->SetBottomMargin(0.3);
    pad1->SetLeftMargin(0.15);
    pad1->Draw();
    pad1->cd();
    r = (TH1D*)pythia_pt_xnorm_rebin_xnorm->Clone("r");
    r->Divide(data_pt_xnorm_rebin_xnorm);
   
    r->SetMarkerStyle(24);
    r->SetMarkerColor(kBlack);
    r->SetMarkerSize(msize);
    r->SetFillColorAlpha(kBlack,0.7);
    r->SetStats(0);
    r->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
    r->GetXaxis()->SetTitleSize(0.10);
    r->GetYaxis()->SetTitle("Pythia / data");
    r->GetYaxis()->SetTitleSize(0.10);
    r->GetXaxis()->SetLabelSize(0.08);
    r->GetYaxis()->SetLabelSize(0.08);
    r->SetTitle("");
    r->Draw("ep");

    c2->SaveAs("/home/clayton/Analysis/code/normalizeAndCompareSpectra/figures/pythiaDataPtCompare_findParams.pdf");
    
    func->SetParNames("a","#alpha","b","#beta");
    TH1D* r_raw = (TH1D*) pythia_pt->Clone("r_raw");
    r_raw->Divide(data_pt);
    r_raw->Fit("func","N","MR",100.0,500.0);
    

   





    double a = func->GetParameter(0);
    double e_a = func->GetParError(0);
    double alpha = func->GetParameter(1);
    double e_alpha = func->GetParError(1);
    double b = func->GetParameter(2);
    double e_b = func->GetParError(2);
    double beta = func->GetParameter(3);
    double e_beta = func->GetParError(3);
    TF1 *fit_fxn = new TF1("fit_fxn","[0]*exp([1]*(x-100)) + [2]*exp([3]*(x-100))",50.0,500.0);
    fit_fxn->SetParameters(a,alpha,b,beta);
    
    TCanvas *c3 = new TCanvas("c3","Reweighted pt vs data",500,500);
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

    TH1D* pythia_pt_xnorm_fitscaled = (TH1D*)pythia_pt_xnorm->Clone("pythia_pt_xnorm_fitscaled");
    int Nbin6 = pythia_pt_xnorm_fitscaled->GetSize();
    for(int k=0;k<Nbin6;k++){
        double pteval = pythia_pt_xnorm_fitscaled->GetXaxis()->GetBinCenter(k);
        double fitscale = fit_fxn->Eval(pteval);
        double oldval = pythia_pt_xnorm_fitscaled->GetBinContent(k);
        double olderror = pythia_pt_xnorm_fitscaled->GetBinError(k);
        double newval = oldval/fitscale;
        double newerror = newval*TMath::Sqrt(TMath::Power(olderror/oldval,2) + TMath::Power(e_a/a,2) + TMath::Power(e_alpha/alpha,2) + TMath::Power(e_b/b,2) + TMath::Power(e_beta/beta,2));
        pythia_pt_xnorm_fitscaled->SetBinContent(k,newval);
        pythia_pt_xnorm_fitscaled->SetBinError(k,newerror);
    }    
   
   pythia_pt_xnorm_fitscaled->SetMarkerColor(kGreen+2);
   pythia_pt_xnorm_fitscaled->SetFillColorAlpha(kGreen+2,0.7);
   pythia_pt_xnorm_fitscaled->SetMarkerStyle(24);
   data_pt_xnorm->GetYaxis()->SetTitle("1/N^{jet} dN^{jet}/dp_{T}");
   data_pt_xnorm->GetYaxis()->SetTitleSize(0.06);
   data_pt_xnorm->SetTitle("#hat{p}_{T} > 50 GeV/c, reweighted");
   data_pt_xnorm->SetStats(0);
   
   data_pt_xnorm->SetMarkerStyle(24);

   TH1D *pythia_pt_xnorm_fitscaled_xnorm = (TH1D*) pythia_pt_xnorm_fitscaled->Rebin(15,"pythia_pt_xnorm_fitscaled_xnorm",pt_axis);
   TH1D *data_pt_xnorm_xnorm = (TH1D*) data_pt_xnorm->Rebin(15,"data_pt_xnorm_xnorm",pt_axis);
   int N_pythia_pt_xnorm_fitscaled_xnorm = pythia_pt_xnorm_fitscaled_xnorm->GetSize();
   for(int l=0;l<N_pythia_pt_xnorm_fitscaled_xnorm;l++){
       double xr2 = pythia_pt_xnorm_fitscaled_xnorm->GetBinWidth(l);
       double yr2 = pythia_pt_xnorm_fitscaled_xnorm->GetBinContent(l);
       double zr2 = pythia_pt_xnorm_fitscaled_xnorm->GetBinError(l);
       double rr2 = 0.0;
       double er2 = 0.0;
       if(xr2!=0){
           rr2 = yr2/xr2;
           er2 = zr2/xr2;
           pythia_pt_xnorm_fitscaled_xnorm->SetBinContent(l,rr2);
           pythia_pt_xnorm_fitscaled_xnorm->SetBinError(l,er2);
       }
   }

   int N_data_pt_xnorm_xnorm = data_pt_xnorm_xnorm->GetSize();
   for(int m=0;m<N_data_pt_xnorm_xnorm;m++){
       double xxr2 = data_pt_xnorm_xnorm->GetBinWidth(m);
       double yyr2 = data_pt_xnorm_xnorm->GetBinContent(m);
       double zzr2 = data_pt_xnorm_xnorm->GetBinError(m);
       double rrr2 = 0.0;
       double eer2 = 0.0;
       if(xxr2!=0){
           rrr2 = yyr2/xxr2;
           eer2 = zzr2/xxr2;
           data_pt_xnorm_xnorm->SetBinContent(m,rrr2);
           data_pt_xnorm_xnorm->SetBinError(m,eer2);
       }
   }
   


   
   data_pt_xnorm_xnorm->Draw("e2");
   pythia_pt_xnorm_fitscaled_xnorm->Draw("e2 same");


   auto legend2 = new TLegend(0.5,0.65,0.89,0.82);
   legend2->AddEntry(data_pt_xnorm_xnorm,"Data","p");
   legend2->AddEntry(pythia_pt_xnorm_fitscaled_xnorm,"Reweighted pythia","p");
   legend2->SetBorderSize(0);
   legend2->Draw();




   pad4->cd();
   TH1D *r2 = (TH1D*) pythia_pt_xnorm_fitscaled_xnorm->Clone("r2");
   r2->Divide(data_pt_xnorm_xnorm);
   r2->SetMarkerColor(kBlack);
   
   r2->Draw("ep");
   r2->SetMinimum(0.8);
   r2->SetMaximum(1.2);
   r2->SetTitle("");
   r2->SetStats(0);
   r2->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
   r2->GetYaxis()->SetTitle("Pythia / data");
   r2->GetXaxis()->SetTitleSize(0.10);
   r2->GetYaxis()->SetTitleSize(0.10);
   r2->GetXaxis()->SetLabelSize(0.08);
   r2->GetYaxis()->SetLabelSize(0.08);
   r2->GetYaxis()->SetTitleOffset(0.6);
   c3->SaveAs("/home/clayton/Analysis/code/normalizeAndCompareSpectra/figures/pythiaDataPtRescaleCompare_findParams.pdf");

   


    

    









}
