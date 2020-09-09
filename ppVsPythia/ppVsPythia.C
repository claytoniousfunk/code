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


void ppVsPythia(){

    //TFile *f_pp = TFile::Open("/home/clayton/Analysis/code/makeData/V2/rootFiles/makeDataV2_pp_data_corrpt_muptcut_10.root");
    //TFile *f_PbPb = TFile::Open("/home/clayton/Analysis/code/makeData/V2/rootFiles/makeDataV2_PbPb_data_corrpt_muptcut_10.root");

    //TFile *f_pp = TFile::Open("/home/clayton/Analysis/data/ppDataSkim_27Aug20/merge.root");
    //TFile *f_PbPb = TFile::Open("/home/clayton/Analysis/data/PbPbDataSkim_27Aug20/merge.root");

    //TFile *f_pp = TFile::Open("/home/clayton/Analysis/data/ppDataSkim_31Aug20/merge.root");
    //TFile *f_PbPb = TFile::Open("/home/clayton/Analysis/data/PbPbDataSkim_31Aug20/merge.root");

    //TFile *f_pp = TFile::Open("/home/clayton/Analysis/code/skimming/pp_mc_skim/pp_mc_skim_noWeights_pthat30_muptcut10_2Sep20.root");
    //TFile *f_PbPb = TFile::Open("/home/clayton/Analysis/data/PbPbDataSkim_4Sep20/PbPb_mc_skim_slim_CsJets_noWeights_pthat30_muptcut10_4Sep20.root");

    TFile *f_pythia = TFile::Open("/home/clayton/Analysis/code/skimming/pp_mc_skim/pp_mc_skim_pthatWeight_pthatcut30_muptcut10_8Sep20.root");
    TFile *f_pp = TFile::Open("/home/clayton/Analysis/code/skimming/pp_data_skim/merge.root");

    TH1D *pp_jetpt, *pp_jetphi, *pp_jeteta, *pp_mujetpt, *pp_mujetphi, *pp_mujeteta, *pp_muPt, *pp_muPhi, *pp_muEta;
    TH1D *pythia_jetpt, *pythia_jetphi, *pythia_jeteta, *pythia_mujetpt, *pythia_mujetphi, *pythia_mujeteta, 
        *pythia_muPt, *pythia_muPhi, *pythia_muEta;

    f_pp->GetObject("h_jetpt",pp_jetpt);
    f_pp->GetObject("h_jetphi",pp_jetphi);
    f_pp->GetObject("h_jeteta",pp_jeteta);
    f_pp->GetObject("h_mujetpt",pp_mujetpt);
    f_pp->GetObject("h_mujetphi",pp_mujetphi);
    f_pp->GetObject("h_mujeteta",pp_mujeteta);
    f_pp->GetObject("h_muPt",pp_muPt);
    f_pp->GetObject("h_muPhi",pp_muPhi);
    f_pp->GetObject("h_muEta",pp_muEta);
    
    f_pythia->GetObject("h_jetpt",pythia_jetpt);
    f_pythia->GetObject("h_jetphi",pythia_jetphi);
    f_pythia->GetObject("h_jeteta",pythia_jeteta);
    f_pythia->GetObject("h_mujetpt",pythia_mujetpt);
    f_pythia->GetObject("h_mujetphi",pythia_mujetphi);
    f_pythia->GetObject("h_mujeteta",pythia_mujeteta);
    f_pythia->GetObject("h_muPt",pythia_muPt);
    f_pythia->GetObject("h_muPhi",pythia_muPhi);
    f_pythia->GetObject("h_muEta",pythia_muEta);
    
    // normalize everything by integral

    TH1D *pp_jetpt_ynorm = (TH1D*) pp_jetpt->Clone("pp_jetpt_ynorm");
    TH1D *pp_jetphi_ynorm = (TH1D*) pp_jetphi->Clone("pp_jetphi_ynorm");
    TH1D *pp_jeteta_ynorm = (TH1D*) pp_jeteta->Clone("pp_jeteta_ynorm");
    TH1D *pp_mujetpt_ynorm = (TH1D*) pp_mujetpt->Clone("pp_mujetpt_ynorm");
    TH1D *pp_mujetphi_ynorm = (TH1D*) pp_mujetphi->Clone("pp_mujetphi_ynorm");
    TH1D *pp_mujeteta_ynorm = (TH1D*) pp_mujeteta->Clone("pp_mujeteta_ynorm");
    TH1D *pp_muPt_ynorm = (TH1D*) pp_muPt->Clone("pp_muPt_ynorm");
    TH1D *pp_muPhi_ynorm = (TH1D*) pp_muPhi->Clone("pp_muPhi_ynorm");
    TH1D *pp_muEta_ynorm = (TH1D*) pp_muEta->Clone("pp_muEta_ynorm");
    
    TH1D *pythia_jetpt_ynorm = (TH1D*) pythia_jetpt->Clone("pythia_jetpt_ynorm");
    TH1D *pythia_jetphi_ynorm = (TH1D*) pythia_jetphi->Clone("pythia_jetphi_ynorm");
    TH1D *pythia_jeteta_ynorm = (TH1D*) pythia_jeteta->Clone("pythia_jeteta_ynorm");
    TH1D *pythia_mujetpt_ynorm = (TH1D*) pythia_mujetpt->Clone("pythia_mujetpt_ynorm");
    TH1D *pythia_mujetphi_ynorm = (TH1D*) pythia_mujetphi->Clone("pythia_mujetphi_ynorm");
    TH1D *pythia_mujeteta_ynorm = (TH1D*) pythia_mujeteta->Clone("pythia_mujeteta_ynorm");
    TH1D *pythia_muPt_ynorm = (TH1D*) pythia_muPt->Clone("pythia_muPt_ynorm");
    TH1D *pythia_muPhi_ynorm = (TH1D*) pythia_muPhi->Clone("pythia_muPhi_ynorm");
    TH1D *pythia_muEta_ynorm = (TH1D*) pythia_muEta->Clone("pythia_muEta_ynorm");
        
    
    pp_jetpt_ynorm->Scale(1.0/pp_jetpt_ynorm->Integral());
    pp_jetphi_ynorm->Scale(1.0/pp_jetphi_ynorm->Integral());
    pp_jeteta_ynorm->Scale(1.0/pp_jeteta_ynorm->Integral());
    pp_mujetphi_ynorm->Scale(1.0/pp_mujetphi_ynorm->Integral());
    pp_mujeteta_ynorm->Scale(1.0/pp_mujeteta_ynorm->Integral());
    pp_muPhi_ynorm->Scale(1.0/pp_muPhi_ynorm->Integral());
    pp_muEta_ynorm->Scale(1.0/pp_muEta_ynorm->Integral());
    
    pythia_jetpt_ynorm->Scale(1.0/pythia_jetpt_ynorm->Integral());
    pythia_jetphi_ynorm->Scale(1.0/pythia_jetphi_ynorm->Integral());
    pythia_jeteta_ynorm->Scale(1.0/pythia_jeteta_ynorm->Integral());
    pythia_mujetphi_ynorm->Scale(1.0/pythia_mujetphi_ynorm->Integral());
    pythia_mujeteta_ynorm->Scale(1.0/pythia_mujeteta_ynorm->Integral());
    pythia_muPhi_ynorm->Scale(1.0/pythia_muPhi_ynorm->Integral());
    pythia_muEta_ynorm->Scale(1.0/pythia_muEta_ynorm->Integral());
    

    // scale by integral
        pp_mujetpt_ynorm->Scale(1.0/pp_mujetpt_ynorm->Integral());
        pythia_mujetpt_ynorm->Scale(1.0/pythia_mujetpt_ynorm->Integral());
        
        pp_muPt_ynorm->Scale(1.0/pp_muPt_ynorm->Integral());
        pythia_muPt_ynorm->Scale(1.0/pythia_muPt_ynorm->Integral());
        

    // define rebin axes
     double pt_axis[16] = {50.0,60.0,70.0,80.0,90.0,100.0,120.0,140.0,160.0,180.0,200.0,250.0,300.0,350.0,400.0,500.0};
     double mupt_axis[21] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,120.0,140.0,160.0,180.0,200.0,250.0,300.0,350.0,400.0,500.0};
     double eta_axis[27] = {-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3};
     double phi_axis[13] = {-150*TMath::Pi()/150,-125*TMath::Pi()/150,-100*TMath::Pi()/150,-75*TMath::Pi()/150,-50*TMath::Pi()/150,-25*TMath::Pi()/150,0.0,25*TMath::Pi()/150,
     50*TMath::Pi()/150,75*TMath::Pi()/150,100*TMath::Pi()/150,125*TMath::Pi()/150,150*TMath::Pi()/150};

    // create rebinned histograms
    TH1D *pp_jetpt_ynorm_rebin = (TH1D*) pp_jetpt_ynorm->Rebin(15,"pp_jetpt_ynorm_rebin",pt_axis);
    TH1D *pythia_jetpt_ynorm_rebin = (TH1D*) pythia_jetpt_ynorm->Rebin(15,"pythia_jetpt_ynorm_rebin",pt_axis);
    
    TH1D *pp_mujetpt_ynorm_rebin = (TH1D*) pp_mujetpt_ynorm->Rebin(15,"pp_mujetpt_ynorm_rebin",pt_axis);
    TH1D *pythia_mujetpt_ynorm_rebin = (TH1D*) pythia_mujetpt_ynorm->Rebin(15,"pythia_mujetpt_ynorm_rebin",pt_axis);
    
    TH1D *pp_muPt_ynorm_rebin = (TH1D*) pp_muPt_ynorm->Rebin(20,"pp_muPt_ynorm_rebin",mupt_axis);
    TH1D *pythia_muPt_ynorm_rebin = (TH1D*) pythia_muPt_ynorm->Rebin(20,"pythia_muPt_ynorm_rebin",mupt_axis);


    TH1D *pp_jeteta_ynorm_rebin = (TH1D*) pp_jeteta_ynorm->Rebin(26,"pp_jeteta_ynorm_rebin",eta_axis);
    TH1D *pp_mujeteta_ynorm_rebin = (TH1D*) pp_mujeteta_ynorm->Rebin(26,"pp_mujeteta_ynorm_rebin",eta_axis);
    TH1D *pp_muEta_ynorm_rebin = (TH1D*) pp_muEta_ynorm->Rebin(26,"pp_muEta_ynorm_rebin",eta_axis);

    TH1D *pythia_jeteta_ynorm_rebin = (TH1D*) pythia_jeteta_ynorm->Rebin(26,"pythia_jeteta_ynorm_rebin",eta_axis);
    TH1D *pythia_mujeteta_ynorm_rebin = (TH1D*) pythia_mujeteta_ynorm->Rebin(26,"pythia_mujeteta_ynorm_rebin",eta_axis);
    TH1D *pythia_muEta_ynorm_rebin = (TH1D*) pythia_muEta_ynorm->Rebin(26,"pythia_muEta_ynorm_rebin",eta_axis);
  
    TH1D *pp_jetphi_ynorm_rebin = (TH1D*) pp_jetphi_ynorm->Rebin(12,"pp_jetphi_ynorm_rebin",phi_axis);
    TH1D *pp_mujetphi_ynorm_rebin = (TH1D*) pp_mujetphi_ynorm->Rebin(12,"pp_mujetphi_ynorm_rebin",phi_axis);
    TH1D *pp_muPhi_ynorm_rebin = (TH1D*) pp_muPhi_ynorm->Rebin(12,"pp_muPhi_ynorm_rebin",phi_axis);

    TH1D *pythia_jetphi_ynorm_rebin = (TH1D*) pythia_jetphi_ynorm->Rebin(12,"pythia_jetphi_ynorm_rebin",phi_axis);
    TH1D *pythia_mujetphi_ynorm_rebin = (TH1D*) pythia_mujetphi_ynorm->Rebin(12,"pythia_mujetphi_ynorm_rebin",phi_axis);
    TH1D *pythia_muPhi_ynorm_rebin = (TH1D*) pythia_muPhi_ynorm->Rebin(12,"pythia_muPhi_ynorm_rebin",phi_axis);
   

    // normalize by bin width
     // pt loop

    TH1D *pp_jetpt_ynorm_rebin_xnorm = (TH1D*) pp_jetpt_ynorm_rebin->Clone("pp_jetpt_ynorm_rebin_xnorm");
    TH1D *pp_mujetpt_ynorm_rebin_xnorm = (TH1D*) pp_mujetpt_ynorm_rebin->Clone("pp_mujetpt_ynorm_rebin_xnorm");
    TH1D *pp_muPt_ynorm_rebin_xnorm = (TH1D*) pp_muPt_ynorm_rebin->Clone("pp_muPt_ynorm_rebin_xnorm");

    TH1D *pythia_jetpt_ynorm_rebin_xnorm = (TH1D*) pythia_jetpt_ynorm_rebin->Clone("pythia_jetpt_ynorm_rebin_xnorm");
    TH1D *pythia_mujetpt_ynorm_rebin_xnorm = (TH1D*) pythia_mujetpt_ynorm_rebin->Clone("pythia_mujetpt_ynorm_rebin_xnorm");
    TH1D *pythia_muPt_ynorm_rebin_xnorm = (TH1D*) pythia_muPt_ynorm_rebin->Clone("pythia_muPt_ynorm_rebin_xnorm");
   

    int N_pt = pp_jetpt_ynorm_rebin_xnorm->GetSize();
    for(int i=0;i<N_pt;i++){
        double w1 = pp_jetpt_ynorm_rebin_xnorm->GetBinWidth(i);
        double y1 = pp_jetpt_ynorm_rebin_xnorm->GetBinContent(i);
        double e1 = pp_jetpt_ynorm_rebin_xnorm->GetBinError(i);
        if(w1!=0){
            pp_jetpt_ynorm_rebin_xnorm->SetBinContent(i,y1/w1);
            pp_jetpt_ynorm_rebin_xnorm->SetBinError(i,e1/w1);
        }
        double w2 = pp_mujetpt_ynorm_rebin_xnorm->GetBinWidth(i);
        double y2 = pp_mujetpt_ynorm_rebin_xnorm->GetBinContent(i);
        double e2 = pp_mujetpt_ynorm_rebin_xnorm->GetBinError(i);
        if(w2!=0){
            pp_mujetpt_ynorm_rebin_xnorm->SetBinContent(i,y2/w2);
            pp_mujetpt_ynorm_rebin_xnorm->SetBinError(i,e2/w2);
        }
        double w3 = pythia_jetpt_ynorm_rebin_xnorm->GetBinWidth(i);
        double y3 = pythia_jetpt_ynorm_rebin_xnorm->GetBinContent(i);
        double e3 = pythia_jetpt_ynorm_rebin_xnorm->GetBinError(i);
        if(w3!=0){
            pythia_jetpt_ynorm_rebin_xnorm->SetBinContent(i,y3/w3);
            pythia_jetpt_ynorm_rebin_xnorm->SetBinError(i,e3/w3);
        }
        double w4 = pythia_mujetpt_ynorm_rebin_xnorm->GetBinWidth(i);
        double y4 = pythia_mujetpt_ynorm_rebin_xnorm->GetBinContent(i);
        double e4 = pythia_mujetpt_ynorm_rebin_xnorm->GetBinError(i);
        if(w4!=0){
            pythia_mujetpt_ynorm_rebin_xnorm->SetBinContent(i,y4/w4);
            pythia_mujetpt_ynorm_rebin_xnorm->SetBinError(i,e4/w4);
        }
        
    }

    // mupt loop

    int N_muPt = pp_muPt_ynorm_rebin_xnorm->GetSize();
    for(int i=0;i<N_muPt;i++){
        double w1 = pp_muPt_ynorm_rebin_xnorm->GetBinWidth(i);
        double y1 = pp_muPt_ynorm_rebin_xnorm->GetBinContent(i);
        double e1 = pp_muPt_ynorm_rebin_xnorm->GetBinError(i);
        if(w1!=0){
            pp_muPt_ynorm_rebin_xnorm->SetBinContent(i,y1/w1);
            pp_muPt_ynorm_rebin_xnorm->SetBinError(i,e1/w1);
        }
        double w2 = pythia_muPt_ynorm_rebin_xnorm->GetBinWidth(i);
        double y2 = pythia_muPt_ynorm_rebin_xnorm->GetBinContent(i);
        double e2 = pythia_muPt_ynorm_rebin_xnorm->GetBinError(i);
        if(w2!=0){
            pythia_muPt_ynorm_rebin_xnorm->SetBinContent(i,y2/w2);
            pythia_muPt_ynorm_rebin_xnorm->SetBinError(i,e2/w2);
        }
            
    }

    // eta loop

    TH1D *pp_jeteta_ynorm_rebin_xnorm = (TH1D*) pp_jeteta_ynorm_rebin->Clone("pp_jeteta_ynorm_rebin_xnorm");
    TH1D *pp_mujeteta_ynorm_rebin_xnorm = (TH1D*) pp_mujeteta_ynorm_rebin->Clone("pp_mujeteta_ynorm_rebin_xnorm");
    TH1D *pp_muEta_ynorm_rebin_xnorm = (TH1D*) pp_muEta_ynorm_rebin->Clone("pp_muEta_ynorm_rebin_xnorm");

    TH1D *pythia_jeteta_ynorm_rebin_xnorm = (TH1D*) pythia_jeteta_ynorm_rebin->Clone("pp_jeteta_ynorm_rebin_xnorm");
    TH1D *pythia_mujeteta_ynorm_rebin_xnorm = (TH1D*) pythia_mujeteta_ynorm_rebin->Clone("pp_mujeteta_ynorm_rebin_xnorm");
    TH1D *pythia_muEta_ynorm_rebin_xnorm = (TH1D*) pythia_muEta_ynorm_rebin->Clone("pp_muEta_ynorm_rebin_xnorm");
    

    int N_eta = pp_jeteta_ynorm_rebin_xnorm->GetSize();
    for(int i=0;i<N_eta;i++){
        double w1 = pp_jeteta_ynorm_rebin_xnorm->GetBinWidth(i);
        double y1 = pp_jeteta_ynorm_rebin_xnorm->GetBinContent(i);
        double e1 = pp_jeteta_ynorm_rebin_xnorm->GetBinError(i);
        if(w1!=0){
            pp_jeteta_ynorm_rebin_xnorm->SetBinContent(i,y1/w1);
            pp_jeteta_ynorm_rebin_xnorm->SetBinError(i,e1/w1);
        }
        double w2 = pp_mujeteta_ynorm_rebin_xnorm->GetBinWidth(i);
        double y2 = pp_mujeteta_ynorm_rebin_xnorm->GetBinContent(i);
        double e2 = pp_mujeteta_ynorm_rebin_xnorm->GetBinError(i);
        if(w2!=0){
            pp_mujeteta_ynorm_rebin_xnorm->SetBinContent(i,y2/w2);
            pp_mujeteta_ynorm_rebin_xnorm->SetBinError(i,e2/w2);
        }
        double w3 = pp_muEta_ynorm_rebin_xnorm->GetBinWidth(i);
        double y3 = pp_muEta_ynorm_rebin_xnorm->GetBinContent(i);
        double e3 = pp_muEta_ynorm_rebin_xnorm->GetBinError(i);
        if(w3!=0){
            pp_muEta_ynorm_rebin_xnorm->SetBinContent(i,y3/w3);
            pp_muEta_ynorm_rebin_xnorm->SetBinError(i,e3/w3);
        }
        double w4 = pythia_jeteta_ynorm_rebin_xnorm->GetBinWidth(i);
        double y4 = pythia_jeteta_ynorm_rebin_xnorm->GetBinContent(i);
        double e4 = pythia_jeteta_ynorm_rebin_xnorm->GetBinError(i);
        if(w4!=0){
            pythia_jeteta_ynorm_rebin_xnorm->SetBinContent(i,y4/w4);
            pythia_jeteta_ynorm_rebin_xnorm->SetBinError(i,e4/w4);
        }
        double w5 = pythia_mujeteta_ynorm_rebin_xnorm->GetBinWidth(i);
        double y5 = pythia_mujeteta_ynorm_rebin_xnorm->GetBinContent(i);
        double e5 = pythia_mujeteta_ynorm_rebin_xnorm->GetBinError(i);
        if(w5!=0){
            pythia_mujeteta_ynorm_rebin_xnorm->SetBinContent(i,y5/w5);
            pythia_mujeteta_ynorm_rebin_xnorm->SetBinError(i,e5/w5);
        }
        double w6 = pythia_muEta_ynorm_rebin_xnorm->GetBinWidth(i);
        double y6 = pythia_muEta_ynorm_rebin_xnorm->GetBinContent(i);
        double e6 = pythia_muEta_ynorm_rebin_xnorm->GetBinError(i);
        if(w6!=0){
            pythia_muEta_ynorm_rebin_xnorm->SetBinContent(i,y6/w6);
            pythia_muEta_ynorm_rebin_xnorm->SetBinError(i,e6/w6);
        }
        
    }

     // phi loop
    TH1D *pp_jetphi_ynorm_rebin_xnorm = (TH1D*) pp_jetphi_ynorm_rebin->Clone("pp_jetphi_ynorm_rebin_xnorm");
    TH1D *pp_mujetphi_ynorm_rebin_xnorm = (TH1D*) pp_mujetphi_ynorm_rebin->Clone("pp_mujetphi_ynorm_rebin_xnorm");
    TH1D *pp_muPhi_ynorm_rebin_xnorm = (TH1D*) pp_muPhi_ynorm_rebin->Clone("pp_muPhi_ynorm_rebin_xnorm");

    TH1D *pythia_jetphi_ynorm_rebin_xnorm = (TH1D*) pythia_jetphi_ynorm_rebin->Clone("pythia_jetphi_ynorm_rebin_xnorm");
    TH1D *pythia_mujetphi_ynorm_rebin_xnorm = (TH1D*) pythia_mujetphi_ynorm_rebin->Clone("pythia_mujetphi_ynorm_rebin_xnorm");
    TH1D *pythia_muPhi_ynorm_rebin_xnorm = (TH1D*) pythia_muPhi_ynorm_rebin->Clone("pythia_muPhi_ynorm_rebin_xnorm");
    

    int N_phi = pp_jetphi_ynorm_rebin_xnorm->GetSize();
    for(int i=0;i<N_phi;i++){
        double w1 = pp_jetphi_ynorm_rebin_xnorm->GetBinWidth(i);
        double y1 = pp_jetphi_ynorm_rebin_xnorm->GetBinContent(i);
        double e1 = pp_jetphi_ynorm_rebin_xnorm->GetBinError(i);
        if(w1!=0){
            pp_jetphi_ynorm_rebin_xnorm->SetBinContent(i,y1/w1);
            pp_jetphi_ynorm_rebin_xnorm->SetBinError(i,e1/w1);
        }
        double w2 = pp_mujetphi_ynorm_rebin_xnorm->GetBinWidth(i);
        double y2 = pp_mujetphi_ynorm_rebin_xnorm->GetBinContent(i);
        double e2 = pp_mujetphi_ynorm_rebin_xnorm->GetBinError(i);
        if(w2!=0){
            pp_mujetphi_ynorm_rebin_xnorm->SetBinContent(i,y2/w2);
            pp_mujetphi_ynorm_rebin_xnorm->SetBinError(i,e2/w2);
        }
        double w3 = pp_muPhi_ynorm_rebin_xnorm->GetBinWidth(i);
        double y3 = pp_muPhi_ynorm_rebin_xnorm->GetBinContent(i);
        double e3 = pp_muPhi_ynorm_rebin_xnorm->GetBinError(i);
        if(w3!=0){
            pp_muPhi_ynorm_rebin_xnorm->SetBinContent(i,y3/w3);
            pp_muPhi_ynorm_rebin_xnorm->SetBinError(i,e3/w3);
        }
        double w4 = pythia_jetphi_ynorm_rebin_xnorm->GetBinWidth(i);
        double y4 = pythia_jetphi_ynorm_rebin_xnorm->GetBinContent(i);
        double e4 = pythia_jetphi_ynorm_rebin_xnorm->GetBinError(i);
        if(w4!=0){
            pythia_jetphi_ynorm_rebin_xnorm->SetBinContent(i,y4/w4);
            pythia_jetphi_ynorm_rebin_xnorm->SetBinError(i,e4/w4);
        }
        double w5 = pythia_mujetphi_ynorm_rebin_xnorm->GetBinWidth(i);
        double y5 = pythia_mujetphi_ynorm_rebin_xnorm->GetBinContent(i);
        double e5 = pythia_mujetphi_ynorm_rebin_xnorm->GetBinError(i);
        if(w5!=0){
            pythia_mujetphi_ynorm_rebin_xnorm->SetBinContent(i,y5/w5);
            pythia_mujetphi_ynorm_rebin_xnorm->SetBinError(i,e5/w5);
        }
        double w6 = pythia_muPhi_ynorm_rebin_xnorm->GetBinWidth(i);
        double y6 = pythia_muPhi_ynorm_rebin_xnorm->GetBinContent(i);
        double e6 = pythia_muPhi_ynorm_rebin_xnorm->GetBinError(i);
        if(w6!=0){
            pythia_muPhi_ynorm_rebin_xnorm->SetBinContent(i,y6/w6);
            pythia_muPhi_ynorm_rebin_xnorm->SetBinError(i,e6/w6);
        }
       
    }
 double msize = 0.9;



      TCanvas *c1 = new TCanvas("c1", "jet pt", 500, 500); 
      c1->cd();
      TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.3, 1., 1.);
      pad2->SetLeftMargin(0.15);
      pad2->SetLogy();
      pad2->Draw();             // Draw the upper pad: pad1
      pad2->cd();
      // jet pT 
      pp_jetpt_ynorm_rebin_xnorm->SetMarkerStyle(8);
      pp_jetpt_ynorm_rebin_xnorm->SetMarkerColor(kRed);
      pp_jetpt_ynorm_rebin_xnorm->SetMarkerSize(msize);
      pp_jetpt_ynorm_rebin_xnorm->SetFillColorAlpha(kRed,0.7);
      pp_jetpt_ynorm_rebin_xnorm->SetStats(0);
      pp_jetpt_ynorm_rebin_xnorm->SetTitle("");

      pythia_jetpt_ynorm_rebin_xnorm->SetMarkerStyle(8);
      pythia_jetpt_ynorm_rebin_xnorm->SetMarkerColor(kBlue);
      pythia_jetpt_ynorm_rebin_xnorm->SetMarkerSize(msize);
      pythia_jetpt_ynorm_rebin_xnorm->SetFillColorAlpha(kBlue,0.7);
      pythia_jetpt_ynorm_rebin_xnorm->SetStats(0);
      pythia_jetpt_ynorm_rebin_xnorm->SetTitle("");

     


      pp_jetpt_ynorm_rebin_xnorm->Draw("e2");
      pythia_jetpt_ynorm_rebin_xnorm->Draw("e2 same");
      

      auto legend = new TLegend(0.7,0.6,0.85,0.85);
      legend->AddEntry(pp_jetpt_ynorm_rebin_xnorm,"pp data","p");
      legend->AddEntry(pythia_jetpt_ynorm_rebin_xnorm,"pp MC","p");
      legend->SetBorderSize(0);
      legend->Draw();


      pp_jetpt_ynorm_rebin_xnorm->GetYaxis()->SetTitle("1/N^{jet}_{tot} dN^{jet}/dp^{jet}_{T}");
      pp_jetpt_ynorm_rebin_xnorm->GetYaxis()->SetTitleSize(0.06);
      pp_jetpt_ynorm_rebin_xnorm->SetMinimum(1e-7);
      pp_jetpt_ynorm_rebin_xnorm->SetMaximum(1e-1);
      
	


///////////////////////////////////////////////////////////////////////////  ratio plot  ////////////////////////////////////////////////////////////////////////////

    c1->cd();
    TPad *pad1 = new TPad("pad1","pad1",0.,0.0,1.,0.3);
    pad1->SetBottomMargin(0.3);
    pad1->SetLeftMargin(0.15);
    pad1->Draw();
    pad1->cd();
    TH1D *r1 = (TH1D*)pythia_jetpt_ynorm_rebin_xnorm->Clone("r1"); 
    r1->Divide(pp_jetpt_ynorm_rebin_xnorm);
    r1->SetMarkerStyle(24);
    r1->SetMarkerColor(kBlack);
    r1->SetFillColorAlpha(kBlack,0.7);
    r1->SetStats(0);
    r1->SetMinimum(0.0);
    r1->SetMaximum(2.0);
    r1->GetXaxis()->SetTitle("p^{jet}_{T} [GeV/c]");
    r1->GetXaxis()->SetTitleSize(0.10);
    r1->GetXaxis()->SetLabelSize(0.08);
    r1->GetYaxis()->SetTitle("MC / data");
    r1->GetYaxis()->SetTitleSize(0.10);
    r1->GetYaxis()->SetLabelSize(0.08);
    r1->SetTitle("");
    TLine *line = new TLine(50,1,500,1);
    line->SetLineStyle(7);
    r1->Draw("e2");
    line->Draw("same");
    

    c1->SaveAs("/home/clayton/Analysis/code/ppVsPythia/figures/ppVsPythia_alljets_pt.pdf");

/*
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// muon jet pt
    TCanvas *c2 = new TCanvas("c2","muon jet pt",500,500);
    c2->cd();
    TPad *pad3 = new TPad("pad3", "pad3", 0.0, 0.3, 1., 1.);
    pad3->SetLeftMargin(0.15);
    pad3->SetLogy();
    pad3->Draw();             // Draw the upper pad: pad1
    pad3->cd();

    pp_mujetpt_ynorm_rebin_xnorm->SetMarkerStyle(8);
    pp_mujetpt_ynorm_rebin_xnorm->SetMarkerColor(kBlack);
    pp_mujetpt_ynorm_rebin_xnorm->SetMarkerSize(msize);
    pp_mujetpt_ynorm_rebin_xnorm->SetFillColorAlpha(kBlack,0.7);
    pp_mujetpt_ynorm_rebin_xnorm->SetStats(0);
    pp_mujetpt_ynorm_rebin_xnorm->SetTitle("");

    PbPb_mujetpt_cent50to90_ynorm_rebin_xnorm->SetMarkerStyle(8);
    PbPb_mujetpt_cent50to90_ynorm_rebin_xnorm->SetMarkerColor(kBlue);
    PbPb_mujetpt_cent50to90_ynorm_rebin_xnorm->SetMarkerSize(msize);
    PbPb_mujetpt_cent50to90_ynorm_rebin_xnorm->SetFillColorAlpha(kBlue,0.7);
    PbPb_mujetpt_cent50to90_ynorm_rebin_xnorm->SetStats(0);
    PbPb_mujetpt_cent50to90_ynorm_rebin_xnorm->SetTitle("");

    PbPb_mujetpt_cent30to50_ynorm_rebin_xnorm->SetMarkerStyle(8);
    PbPb_mujetpt_cent30to50_ynorm_rebin_xnorm->SetMarkerColor(kGreen+2);
    PbPb_mujetpt_cent30to50_ynorm_rebin_xnorm->SetMarkerSize(msize);
    PbPb_mujetpt_cent30to50_ynorm_rebin_xnorm->SetFillColorAlpha(kGreen+2,0.7);
    PbPb_mujetpt_cent30to50_ynorm_rebin_xnorm->SetStats(0);
    PbPb_mujetpt_cent30to50_ynorm_rebin_xnorm->SetTitle("");

    PbPb_mujetpt_cent10to30_ynorm_rebin_xnorm->SetMarkerStyle(8);
    PbPb_mujetpt_cent10to30_ynorm_rebin_xnorm->SetMarkerColor(kRed);
    PbPb_mujetpt_cent10to30_ynorm_rebin_xnorm->SetMarkerSize(msize);
    PbPb_mujetpt_cent10to30_ynorm_rebin_xnorm->SetFillColorAlpha(kRed,0.7);
    PbPb_mujetpt_cent10to30_ynorm_rebin_xnorm->SetStats(0);
    PbPb_mujetpt_cent10to30_ynorm_rebin_xnorm->SetTitle(""); 

    PbPb_mujetpt_cent0to10_ynorm_rebin_xnorm->SetMarkerStyle(8);
    PbPb_mujetpt_cent0to10_ynorm_rebin_xnorm->SetMarkerColor(kBlue-10);
    PbPb_mujetpt_cent0to10_ynorm_rebin_xnorm->SetMarkerSize(msize);
    PbPb_mujetpt_cent0to10_ynorm_rebin_xnorm->SetFillColorAlpha(kBlue-10,0.7);
    PbPb_mujetpt_cent0to10_ynorm_rebin_xnorm->SetStats(0);
    PbPb_mujetpt_cent0to10_ynorm_rebin_xnorm->SetTitle("");

    pp_mujetpt_ynorm_rebin_xnorm->Draw("e2");
    PbPb_mujetpt_cent50to90_ynorm_rebin_xnorm->Draw("e2 same");
    PbPb_mujetpt_cent30to50_ynorm_rebin_xnorm->Draw("e2 same");
    PbPb_mujetpt_cent10to30_ynorm_rebin_xnorm->Draw("e2 same");
    PbPb_mujetpt_cent0to10_ynorm_rebin_xnorm->Draw("e2 same");
    
  
    auto legend2 = new TLegend(0.7,0.6,0.85,0.85);
    legend2->AddEntry(pp_mujetpt_ynorm_rebin_xnorm,"pp","p");
    legend2->AddEntry(PbPb_mujetpt_cent50to90_ynorm_rebin_xnorm,"50-90 %","p");
    legend2->AddEntry(PbPb_mujetpt_cent30to50_ynorm_rebin_xnorm,"30-50 %","p");
    legend2->AddEntry(PbPb_mujetpt_cent10to30_ynorm_rebin_xnorm,"10-30 %","p");
    legend2->AddEntry(PbPb_mujetpt_cent0to10_ynorm_rebin_xnorm,"0-10 %","p");
    legend2->SetBorderSize(0);
    legend2->Draw();


    pp_mujetpt_ynorm_rebin_xnorm->GetYaxis()->SetTitle("1/N^{#mu-jet}_{tot} dN^{#mu-jet}/dp^{jet}_{T}");
    pp_mujetpt_ynorm_rebin_xnorm->GetYaxis()->SetTitleSize(0.06);
   // pp_mujetpt_ynorm_rebin_xnorm->SetMinimum(1e-8);
   // pp_mujetpt_ynorm_rebin_xnorm->SetMaximum(1e4);

    c2->cd();
    TPad *pad4 = new TPad("pad4","pad4",0.,0.0,1.,0.3);
    pad4->SetBottomMargin(0.3);
    pad4->SetLeftMargin(0.15);
    pad4->Draw();
    pad4->cd();
    TH1D *r5 = (TH1D*)PbPb_mujetpt_cent0to10_ynorm_rebin_xnorm->Clone("r5"); 
    TH1D *r6 = (TH1D*)PbPb_mujetpt_cent50to90_ynorm_rebin_xnorm->Clone("r6");
    TH1D *r7 = (TH1D*)PbPb_mujetpt_cent30to50_ynorm_rebin_xnorm->Clone("r7");
    TH1D *r8 = (TH1D*)PbPb_mujetpt_cent10to30_ynorm_rebin_xnorm->Clone("r8");
    r5->Divide(pp_mujetpt_ynorm_rebin_xnorm);
    r6->Divide(pp_mujetpt_ynorm_rebin_xnorm);
    r7->Divide(pp_mujetpt_ynorm_rebin_xnorm);
    r8->Divide(pp_mujetpt_ynorm_rebin_xnorm);
    r5->SetMarkerStyle(24);
    r6->SetMarkerStyle(24);
    r7->SetMarkerStyle(24);
    r8->SetMarkerStyle(24);
    r5->SetStats(0);
    r5->SetMinimum(0.0);
    r5->SetMaximum(2.0);
    r5->GetXaxis()->SetTitle("p^{jet}_{T} [GeV/c]");
    r5->GetXaxis()->SetTitleSize(0.10);
    r5->GetXaxis()->SetLabelSize(0.08);
    r5->GetYaxis()->SetTitle("Centrality region / pp");
    r5->GetYaxis()->SetTitleSize(0.10);
    r5->GetYaxis()->SetLabelSize(0.08);
    r5->SetTitle("");
    TLine *line2 = new TLine(50.0,1.0,500.0,1.0);
    line2->SetLineStyle(7);
    r5->Draw("e2");
    r6->Draw("e2 same");
    r7->Draw("e2 same");
    r8->Draw("e2 same");
    line2->Draw("same");
    

    c2->SaveAs("/home/clayton/Analysis/code/normalizeAndCompareSpectra/figures/MCPPvsPbPbCentralityRegions_mujets_pt.pdf");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TCanvas *c3 = new TCanvas("c3","muon jet eta",500,500);
    c3->cd();
    TPad *pad5 = new TPad("pad5", "pad5", 0.0, 0.3, 1., 1.);
    pad5->SetLeftMargin(0.15);
    pad5->SetLogy();
    pad5->Draw();             // Draw the upper pad: pad1
    pad5->cd();

    pp_mujeteta_ynorm_rebin_xnorm->SetMarkerStyle(8);
    pp_mujeteta_ynorm_rebin_xnorm->SetMarkerColor(kBlack);
    pp_mujeteta_ynorm_rebin_xnorm->SetMarkerSize(msize);
    pp_mujeteta_ynorm_rebin_xnorm->SetFillColorAlpha(kBlack,0.7);
    pp_mujeteta_ynorm_rebin_xnorm->SetStats(0);
    pp_mujeteta_ynorm_rebin_xnorm->SetTitle("");

    PbPb_mujeteta_cent50to90_ynorm_rebin_xnorm->SetMarkerStyle(8);
    PbPb_mujeteta_cent50to90_ynorm_rebin_xnorm->SetMarkerColor(kBlue);
    PbPb_mujeteta_cent50to90_ynorm_rebin_xnorm->SetMarkerSize(msize);
    PbPb_mujeteta_cent50to90_ynorm_rebin_xnorm->SetFillColorAlpha(kBlue,0.7);
    PbPb_mujeteta_cent50to90_ynorm_rebin_xnorm->SetStats(0);
    PbPb_mujeteta_cent50to90_ynorm_rebin_xnorm->SetTitle("");

    PbPb_mujeteta_cent30to50_ynorm_rebin_xnorm->SetMarkerStyle(8);
    PbPb_mujeteta_cent30to50_ynorm_rebin_xnorm->SetMarkerColor(kGreen+2);
    PbPb_mujeteta_cent30to50_ynorm_rebin_xnorm->SetMarkerSize(msize);
    PbPb_mujeteta_cent30to50_ynorm_rebin_xnorm->SetFillColorAlpha(kGreen+2,0.7);
    PbPb_mujeteta_cent30to50_ynorm_rebin_xnorm->SetStats(0);
    PbPb_mujeteta_cent30to50_ynorm_rebin_xnorm->SetTitle("");

    PbPb_mujeteta_cent10to30_ynorm_rebin_xnorm->SetMarkerStyle(8);
    PbPb_mujeteta_cent10to30_ynorm_rebin_xnorm->SetMarkerColor(kRed);
    PbPb_mujeteta_cent10to30_ynorm_rebin_xnorm->SetMarkerSize(msize);
    PbPb_mujeteta_cent10to30_ynorm_rebin_xnorm->SetFillColorAlpha(kRed,0.7);
    PbPb_mujeteta_cent10to30_ynorm_rebin_xnorm->SetStats(0);
    PbPb_mujeteta_cent10to30_ynorm_rebin_xnorm->SetTitle(""); 

    PbPb_mujeteta_cent0to10_ynorm_rebin_xnorm->SetMarkerStyle(8);
    PbPb_mujeteta_cent0to10_ynorm_rebin_xnorm->SetMarkerColor(kBlue-10);
    PbPb_mujeteta_cent0to10_ynorm_rebin_xnorm->SetMarkerSize(msize);
    PbPb_mujeteta_cent0to10_ynorm_rebin_xnorm->SetFillColorAlpha(kBlue-10,0.7);
    PbPb_mujeteta_cent0to10_ynorm_rebin_xnorm->SetStats(0);
    PbPb_mujeteta_cent0to10_ynorm_rebin_xnorm->SetTitle("");

    pp_mujeteta_ynorm_rebin_xnorm->Draw("e2");
    PbPb_mujeteta_cent50to90_ynorm_rebin_xnorm->Draw("e2 same");
    PbPb_mujeteta_cent30to50_ynorm_rebin_xnorm->Draw("e2 same");
    PbPb_mujeteta_cent10to30_ynorm_rebin_xnorm->Draw("e2 same");
    PbPb_mujeteta_cent0to10_ynorm_rebin_xnorm->Draw("e2 same");
    
  
    auto legend3 = new TLegend(0.75,0.65,0.89,0.89);
    legend3->AddEntry(pp_mujeteta_ynorm_rebin_xnorm,"pp","p");
    legend3->AddEntry(PbPb_mujeteta_cent50to90_ynorm_rebin_xnorm,"50-90 %","p");
    legend3->AddEntry(PbPb_mujeteta_cent30to50_ynorm_rebin_xnorm,"30-50 %","p");
    legend3->AddEntry(PbPb_mujeteta_cent10to30_ynorm_rebin_xnorm,"10-30 %","p");
    legend3->AddEntry(PbPb_mujeteta_cent0to10_ynorm_rebin_xnorm,"0-10 %","p");
    legend3->SetBorderSize(0);
    legend3->Draw();


    pp_mujeteta_ynorm_rebin_xnorm->GetYaxis()->SetTitle("1/N^{#mu-jet}_{tot} dN^{#mu-jet}/d#eta^{jet}");
    pp_mujeteta_ynorm_rebin_xnorm->GetYaxis()->SetTitleSize(0.06);
    pp_mujeteta_ynorm_rebin_xnorm->SetMinimum(1e-1);
    pp_mujeteta_ynorm_rebin_xnorm->SetMaximum(1e0);

    c3->cd();
    TPad *pad6 = new TPad("pad6","pad6",0.,0.0,1.,0.3);
    pad6->SetBottomMargin(0.3);
    pad6->SetLeftMargin(0.15);
    pad6->Draw();
    pad6->cd();
    TH1D *r9 = (TH1D*)PbPb_mujeteta_cent0to10_ynorm_rebin_xnorm->Clone("r9"); 
    TH1D *r10 = (TH1D*)PbPb_mujeteta_cent50to90_ynorm_rebin_xnorm->Clone("r10");
    TH1D *r11 = (TH1D*)PbPb_mujeteta_cent30to50_ynorm_rebin_xnorm->Clone("r11");
    TH1D *r12 = (TH1D*)PbPb_mujeteta_cent10to30_ynorm_rebin_xnorm->Clone("r12");
    r9->Divide(pp_mujeteta_ynorm_rebin_xnorm);
    r10->Divide(pp_mujeteta_ynorm_rebin_xnorm);
    r11->Divide(pp_mujeteta_ynorm_rebin_xnorm);
    r12->Divide(pp_mujeteta_ynorm_rebin_xnorm);
    r9->SetMarkerStyle(24);
    r10->SetMarkerStyle(24);
    r11->SetMarkerStyle(24);
    r12->SetMarkerStyle(24);
    r9->SetStats(0);
    r9->SetMinimum(0.0);
    r9->SetMaximum(2.0);
    r9->GetXaxis()->SetTitle("#eta^{jet}");
    r9->GetXaxis()->SetTitleSize(0.10);
    r9->GetXaxis()->SetLabelSize(0.08);
    r9->GetYaxis()->SetTitle("Centrality region / pp");
    r9->GetYaxis()->SetTitleSize(0.10);
    r9->GetYaxis()->SetLabelSize(0.08);
    r9->SetTitle("");
    TLine *line3 = new TLine(-1.5,1.0,1.5,1.0);
    line3->SetLineStyle(7);
    r9->Draw("e2");
    r10->Draw("e2 same");
    r11->Draw("e2 same");
    r12->Draw("e2 same");
    line3->Draw("same");
    

    c3->SaveAs("/home/clayton/Analysis/code/normalizeAndCompareSpectra/figures/MCPPvsPbPbCentralityRegions_mujets_eta.pdf");
    


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TCanvas *c4 = new TCanvas("c4","jet eta",500,500);
    c4->cd();
    TPad *pad7 = new TPad("pad7", "pad7", 0.0, 0.3, 1., 1.);
    pad7->SetLeftMargin(0.15);
    pad7->SetLogy();
    pad7->Draw();             // Draw the upper pad: pad1
    pad7->cd();

    pp_jeteta_ynorm_rebin_xnorm->SetMarkerStyle(8);
    pp_jeteta_ynorm_rebin_xnorm->SetMarkerColor(kBlack);
    pp_jeteta_ynorm_rebin_xnorm->SetMarkerSize(msize);
    pp_jeteta_ynorm_rebin_xnorm->SetFillColorAlpha(kBlack,0.7);
    pp_jeteta_ynorm_rebin_xnorm->SetStats(0);
    pp_jeteta_ynorm_rebin_xnorm->SetTitle("");

    PbPb_jeteta_cent50to90_ynorm_rebin_xnorm->SetMarkerStyle(8);
    PbPb_jeteta_cent50to90_ynorm_rebin_xnorm->SetMarkerColor(kBlue);
    PbPb_jeteta_cent50to90_ynorm_rebin_xnorm->SetMarkerSize(msize);
    PbPb_jeteta_cent50to90_ynorm_rebin_xnorm->SetFillColorAlpha(kBlue,0.7);
    PbPb_jeteta_cent50to90_ynorm_rebin_xnorm->SetStats(0);
    PbPb_jeteta_cent50to90_ynorm_rebin_xnorm->SetTitle("");

    PbPb_jeteta_cent30to50_ynorm_rebin_xnorm->SetMarkerStyle(8);
    PbPb_jeteta_cent30to50_ynorm_rebin_xnorm->SetMarkerColor(kGreen+2);
    PbPb_jeteta_cent30to50_ynorm_rebin_xnorm->SetMarkerSize(msize);
    PbPb_jeteta_cent30to50_ynorm_rebin_xnorm->SetFillColorAlpha(kGreen+2,0.7);
    PbPb_jeteta_cent30to50_ynorm_rebin_xnorm->SetStats(0);
    PbPb_jeteta_cent30to50_ynorm_rebin_xnorm->SetTitle("");

    PbPb_jeteta_cent10to30_ynorm_rebin_xnorm->SetMarkerStyle(8);
    PbPb_jeteta_cent10to30_ynorm_rebin_xnorm->SetMarkerColor(kRed);
    PbPb_jeteta_cent10to30_ynorm_rebin_xnorm->SetMarkerSize(msize);
    PbPb_jeteta_cent10to30_ynorm_rebin_xnorm->SetFillColorAlpha(kRed,0.7);
    PbPb_jeteta_cent10to30_ynorm_rebin_xnorm->SetStats(0);
    PbPb_jeteta_cent10to30_ynorm_rebin_xnorm->SetTitle(""); 

    PbPb_jeteta_cent0to10_ynorm_rebin_xnorm->SetMarkerStyle(8);
    PbPb_jeteta_cent0to10_ynorm_rebin_xnorm->SetMarkerColor(kBlue-10);
    PbPb_jeteta_cent0to10_ynorm_rebin_xnorm->SetMarkerSize(msize);
    PbPb_jeteta_cent0to10_ynorm_rebin_xnorm->SetFillColorAlpha(kBlue-10,0.7);
    PbPb_jeteta_cent0to10_ynorm_rebin_xnorm->SetStats(0);
    PbPb_jeteta_cent0to10_ynorm_rebin_xnorm->SetTitle("");

    pp_jeteta_ynorm_rebin_xnorm->Draw("e2");
    PbPb_jeteta_cent50to90_ynorm_rebin_xnorm->Draw("e2 same");
    PbPb_jeteta_cent30to50_ynorm_rebin_xnorm->Draw("e2 same");
    PbPb_jeteta_cent10to30_ynorm_rebin_xnorm->Draw("e2 same");
    PbPb_jeteta_cent0to10_ynorm_rebin_xnorm->Draw("e2 same");
    
  
    auto legend4 = new TLegend(0.75,0.65,0.89,0.89);
    legend4->AddEntry(pp_jeteta_ynorm_rebin_xnorm,"pp","p");
    legend4->AddEntry(PbPb_jeteta_cent50to90_ynorm_rebin_xnorm,"50-90 %","p");
    legend4->AddEntry(PbPb_jeteta_cent30to50_ynorm_rebin_xnorm,"30-50 %","p");
    legend4->AddEntry(PbPb_jeteta_cent10to30_ynorm_rebin_xnorm,"10-30 %","p");
    legend4->AddEntry(PbPb_jeteta_cent0to10_ynorm_rebin_xnorm,"0-10 %","p");
    legend4->SetBorderSize(0);
    legend4->Draw();


    pp_jeteta_ynorm_rebin_xnorm->GetYaxis()->SetTitle("1/N^{jet}_{tot} dN^{jet}/d#eta^{jet}");
    pp_jeteta_ynorm_rebin_xnorm->GetYaxis()->SetTitleSize(0.06);
    pp_jeteta_ynorm_rebin_xnorm->SetMinimum(1e-1);
    pp_jeteta_ynorm_rebin_xnorm->SetMaximum(1e0);

    c4->cd();
    TPad *pad8 = new TPad("pad8","pad8",0.,0.0,1.,0.3);
    pad8->SetBottomMargin(0.3);
    pad8->SetLeftMargin(0.15);
    pad8->Draw();
    pad8->cd();
    TH1D *r13 = (TH1D*)PbPb_jeteta_cent0to10_ynorm_rebin_xnorm->Clone("r13"); 
    TH1D *r14 = (TH1D*)PbPb_jeteta_cent50to90_ynorm_rebin_xnorm->Clone("r14");
    TH1D *r15 = (TH1D*)PbPb_jeteta_cent30to50_ynorm_rebin_xnorm->Clone("r15");
    TH1D *r16 = (TH1D*)PbPb_jeteta_cent10to30_ynorm_rebin_xnorm->Clone("r16");
    r13->Divide(pp_jeteta_ynorm_rebin_xnorm);
    r14->Divide(pp_jeteta_ynorm_rebin_xnorm);
    r15->Divide(pp_jeteta_ynorm_rebin_xnorm);
    r16->Divide(pp_jeteta_ynorm_rebin_xnorm);
    r13->SetMarkerStyle(24);
    r14->SetMarkerStyle(24);
    r15->SetMarkerStyle(24);
    r16->SetMarkerStyle(24);
    r13->SetStats(0);
    r13->SetMinimum(0.0);
    r13->SetMaximum(2.0);
    r13->GetXaxis()->SetTitle("#eta^{jet}");
    r13->GetXaxis()->SetTitleSize(0.10);
    r13->GetXaxis()->SetLabelSize(0.08);
    r13->GetYaxis()->SetTitle("Centrality region / pp");
    r13->GetYaxis()->SetTitleSize(0.10);
    r13->GetYaxis()->SetLabelSize(0.08);
    r13->SetTitle("");
    TLine *line4 = new TLine(-1.5,1.0,1.5,1.0);
    line4->SetLineStyle(7);
    r13->Draw("e2");
    r14->Draw("e2 same");
    r15->Draw("e2 same");
    r16->Draw("e2 same");
    line4->Draw("same");
    

    c4->SaveAs("/home/clayton/Analysis/code/normalizeAndCompareSpectra/figures/MCPPvsPbPbCentralityRegions_alljets_eta.pdf");
    

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TCanvas *c5 = new TCanvas("c5","muon pt",500,500);
    c5->cd();
    TPad *pad9 = new TPad("pad9", "pad9", 0.0, 0.3, 1., 1.);
    pad9->SetLeftMargin(0.15);
    pad9->SetLogy();
    pad9->Draw();             // Draw the upper pad: pad1
    pad9->cd();

    pp_muPt_ynorm_rebin_xnorm->SetMarkerStyle(8);
    pp_muPt_ynorm_rebin_xnorm->SetMarkerColor(kBlack);
    pp_muPt_ynorm_rebin_xnorm->SetMarkerSize(msize);
    pp_muPt_ynorm_rebin_xnorm->SetFillColorAlpha(kBlack,0.7);
    pp_muPt_ynorm_rebin_xnorm->SetStats(0);
    pp_muPt_ynorm_rebin_xnorm->SetTitle("");

    PbPb_muPt_cent50to90_ynorm_rebin_xnorm->SetMarkerStyle(8);
    PbPb_muPt_cent50to90_ynorm_rebin_xnorm->SetMarkerColor(kBlue);
    PbPb_muPt_cent50to90_ynorm_rebin_xnorm->SetMarkerSize(msize);
    PbPb_muPt_cent50to90_ynorm_rebin_xnorm->SetFillColorAlpha(kBlue,0.7);
    PbPb_muPt_cent50to90_ynorm_rebin_xnorm->SetStats(0);
    PbPb_muPt_cent50to90_ynorm_rebin_xnorm->SetTitle("");

    PbPb_muPt_cent30to50_ynorm_rebin_xnorm->SetMarkerStyle(8);
    PbPb_muPt_cent30to50_ynorm_rebin_xnorm->SetMarkerColor(kGreen+2);
    PbPb_muPt_cent30to50_ynorm_rebin_xnorm->SetMarkerSize(msize);
    PbPb_muPt_cent30to50_ynorm_rebin_xnorm->SetFillColorAlpha(kGreen+2,0.7);
    PbPb_muPt_cent30to50_ynorm_rebin_xnorm->SetStats(0);
    PbPb_muPt_cent30to50_ynorm_rebin_xnorm->SetTitle("");

    PbPb_muPt_cent10to30_ynorm_rebin_xnorm->SetMarkerStyle(8);
    PbPb_muPt_cent10to30_ynorm_rebin_xnorm->SetMarkerColor(kRed);
    PbPb_muPt_cent10to30_ynorm_rebin_xnorm->SetMarkerSize(msize);
    PbPb_muPt_cent10to30_ynorm_rebin_xnorm->SetFillColorAlpha(kRed,0.7);
    PbPb_muPt_cent10to30_ynorm_rebin_xnorm->SetStats(0);
    PbPb_muPt_cent10to30_ynorm_rebin_xnorm->SetTitle(""); 

    PbPb_muPt_cent0to10_ynorm_rebin_xnorm->SetMarkerStyle(8);
    PbPb_muPt_cent0to10_ynorm_rebin_xnorm->SetMarkerColor(kBlue-10);
    PbPb_muPt_cent0to10_ynorm_rebin_xnorm->SetMarkerSize(msize);
    PbPb_muPt_cent0to10_ynorm_rebin_xnorm->SetFillColorAlpha(kBlue-10,0.7);
    PbPb_muPt_cent0to10_ynorm_rebin_xnorm->SetStats(0);
    PbPb_muPt_cent0to10_ynorm_rebin_xnorm->SetTitle("");

    pp_muPt_ynorm_rebin_xnorm->Draw("e2");
    PbPb_muPt_cent50to90_ynorm_rebin_xnorm->Draw("e2 same");
    PbPb_muPt_cent30to50_ynorm_rebin_xnorm->Draw("e2 same");
    PbPb_muPt_cent10to30_ynorm_rebin_xnorm->Draw("e2 same");
    PbPb_muPt_cent0to10_ynorm_rebin_xnorm->Draw("e2 same");
    
  
    auto legend5 = new TLegend(0.75,0.65,0.89,0.89);
    legend5->AddEntry(pp_muPt_ynorm_rebin_xnorm,"pp","p");
    legend5->AddEntry(PbPb_muPt_cent50to90_ynorm_rebin_xnorm,"50-90 %","p");
    legend5->AddEntry(PbPb_muPt_cent30to50_ynorm_rebin_xnorm,"30-50 %","p");
    legend5->AddEntry(PbPb_muPt_cent10to30_ynorm_rebin_xnorm,"10-30 %","p");
    legend5->AddEntry(PbPb_muPt_cent0to10_ynorm_rebin_xnorm,"0-10 %","p");
    legend5->SetBorderSize(0);
    legend5->Draw();


    pp_muPt_ynorm_rebin_xnorm->GetYaxis()->SetTitle("1/N^{#mu}_{tot} dN^{#mu}/dp^{#mu}_{T}");
    pp_muPt_ynorm_rebin_xnorm->GetYaxis()->SetTitleSize(0.06);
    pp_muPt_ynorm_rebin_xnorm->SetMinimum(1e-8);
    pp_muPt_ynorm_rebin_xnorm->SetMaximum(1e-1);

    c5->cd();
    TPad *pad10 = new TPad("pad10","pad10",0.,0.0,1.,0.3);
    pad10->SetBottomMargin(0.3);
    pad10->SetLeftMargin(0.15);
    pad10->Draw();
    pad10->cd();
    TH1D *r17 = (TH1D*)PbPb_muPt_cent0to10_ynorm_rebin_xnorm->Clone("r17"); 
    TH1D *r18 = (TH1D*)PbPb_muPt_cent50to90_ynorm_rebin_xnorm->Clone("r18");
    TH1D *r19 = (TH1D*)PbPb_muPt_cent30to50_ynorm_rebin_xnorm->Clone("r19");
    TH1D *r20 = (TH1D*)PbPb_muPt_cent10to30_ynorm_rebin_xnorm->Clone("r20");
    r17->Divide(pp_muPt_ynorm_rebin_xnorm);
    r18->Divide(pp_muPt_ynorm_rebin_xnorm);
    r19->Divide(pp_muPt_ynorm_rebin_xnorm);
    r20->Divide(pp_muPt_ynorm_rebin_xnorm);
    r17->SetMarkerStyle(24);
    r18->SetMarkerStyle(24);
    r19->SetMarkerStyle(24);
    r20->SetMarkerStyle(24);
    r17->SetStats(0);
    r17->SetMinimum(0.0);
    r17->SetMaximum(6.0);
    r17->GetXaxis()->SetTitle("p_{T}^{#mu} [GeV/c]");
    r17->GetXaxis()->SetTitleSize(0.10);
    r17->GetXaxis()->SetLabelSize(0.08);
    r17->GetYaxis()->SetTitle("Centrality region / pp");
    r17->GetYaxis()->SetTitleSize(0.10);
    r17->GetYaxis()->SetLabelSize(0.08);
    r17->SetTitle("");
    TLine *line5 = new TLine(0.0,1.0,500.0,1.0);
    line5->SetLineStyle(7);
    r17->Draw("e2");
    r18->Draw("e2 same");
    r19->Draw("e2 same");
    r20->Draw("e2 same");
    line5->Draw("same");
    

    c5->SaveAs("/home/clayton/Analysis/code/normalizeAndCompareSpectra/figures/MCPPvsPbPbCentralityRegions_muons_pt.pdf");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


*/

    

} // end program