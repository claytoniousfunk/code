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


void DataPPvsPbPbCentralityRegionsPt(){

    //TFile *f_pp = TFile::Open("/home/clayton/Analysis/code/makeData/V2/rootFiles/makeDataV2_pp_data_corrpt_muptcut_10.root");
    //TFile *f_PbPb = TFile::Open("/home/clayton/Analysis/code/makeData/V2/rootFiles/makeDataV2_PbPb_data_corrpt_muptcut_10.root");

    TFile *f_pp = TFile::Open("/home/clayton/Analysis/data/ppDataSkim_27Aug20/merge.root");
    TFile *f_PbPb = TFile::Open("/home/clayton/Analysis/data/PbPbDataSkim_27Aug20/merge.root");

    TH1D *pp_jetpt, *pp_jetphi, *pp_jeteta, *pp_mupt, *pp_muphi, *pp_mueta;
    TH1D *PbPb_jetpt_cent0to10, *PbPb_jetphi_cent0to10, *PbPb_jeteta_cent0to10, *PbPb_mupt_cent0to10, *PbPb_muphi_cent0to10, *PbPb_mueta_cent0to10;
    TH1D *PbPb_jetpt_cent10to30, *PbPb_jetphi_cent10to30, *PbPb_jeteta_cent10to30, *PbPb_mupt_cent10to30, *PbPb_muphi_cent10to30, *PbPb_mueta_cent10to30;
    TH1D *PbPb_jetpt_cent30to50, *PbPb_jetphi_cent30to50, *PbPb_jeteta_cent30to50, *PbPb_mupt_cent30to50, *PbPb_muphi_cent30to50, *PbPb_mueta_cent30to50;
    TH1D *PbPb_jetpt_cent50to90, *PbPb_jetphi_cent50to90, *PbPb_jeteta_cent50to90, *PbPb_mupt_cent50to90, *PbPb_muphi_cent50to90, *PbPb_mueta_cent50to90;

    f_pp->GetObject("h_jetpt",pp_jetpt);
    f_pp->GetObject("h_jetphi",pp_jetphi);
    f_pp->GetObject("h_jeteta",pp_jeteta);
    f_pp->GetObject("h_muPt",pp_mupt);
    f_pp->GetObject("h_muPhi",pp_muphi);
    f_pp->GetObject("h_muEta",pp_mueta);

    f_PbPb->GetObject("h_jetpt_cent0to10",PbPb_jetpt_cent0to10);
    f_PbPb->GetObject("h_jetphi_cent0to10",PbPb_jetphi_cent0to10);
    f_PbPb->GetObject("h_jeteta_cent0to10",PbPb_jeteta_cent0to10);
    f_PbPb->GetObject("h_muPt_cent0to10",PbPb_mupt_cent0to10);
    f_PbPb->GetObject("h_muPhi_cent0to10",PbPb_muphi_cent0to10);
    f_PbPb->GetObject("h_muEta_cent0to10",PbPb_mueta_cent0to10);

    f_PbPb->GetObject("h_jetpt_cent10to30",PbPb_jetpt_cent10to30);
    f_PbPb->GetObject("h_jetphi_cent10to30",PbPb_jetphi_cent10to30);
    f_PbPb->GetObject("h_jeteta_cent10to30",PbPb_jeteta_cent10to30);
    f_PbPb->GetObject("h_muPt_cent10to30",PbPb_mupt_cent10to30);
    f_PbPb->GetObject("h_muPhi_cent10to30",PbPb_muphi_cent10to30);
    f_PbPb->GetObject("h_muEta_cent10to30",PbPb_mueta_cent10to30);

    f_PbPb->GetObject("h_jetpt_cent30to50",PbPb_jetpt_cent30to50);
    f_PbPb->GetObject("h_jetphi_cent30to50",PbPb_jetphi_cent30to50);
    f_PbPb->GetObject("h_jeteta_cent30to50",PbPb_jeteta_cent30to50);
    f_PbPb->GetObject("h_muPt_cent30to50",PbPb_mupt_cent30to50);
    f_PbPb->GetObject("h_muPhi_cent30to50",PbPb_muphi_cent30to50);
    f_PbPb->GetObject("h_muEta_cent30to50",PbPb_mueta_cent30to50);

    f_PbPb->GetObject("h_jetpt_cent50to90",PbPb_jetpt_cent50to90);
    f_PbPb->GetObject("h_jetphi_cent50to90",PbPb_jetphi_cent50to90);
    f_PbPb->GetObject("h_jeteta_cent50to90",PbPb_jeteta_cent50to90);
    f_PbPb->GetObject("h_muPt_cent50to90",PbPb_mupt_cent50to90);
    f_PbPb->GetObject("h_muPhi_cent50to90",PbPb_muphi_cent50to90);
    f_PbPb->GetObject("h_muEta_cent50to90",PbPb_mueta_cent50to90);


    // normalize everything by integral

    TH1D *pp_jetpt_ynorm = (TH1D*) pp_jetpt->Clone("pp_jetpt_ynorm");
    TH1D *pp_jetphi_ynorm = (TH1D*) pp_jetphi->Clone("pp_jetphi_ynorm");
    TH1D *pp_jeteta_ynorm = (TH1D*) pp_jeteta->Clone("pp_jeteta_ynorm");
    TH1D *pp_mupt_ynorm = (TH1D*) pp_mupt->Clone("pp_mupt_ynorm");
    TH1D *pp_muphi_ynorm = (TH1D*) pp_muphi->Clone("pp_muphi_ynorm");
    TH1D *pp_mueta_ynorm = (TH1D*) pp_mueta->Clone("pp_mueta_ynorm");

    TH1D *PbPb_jetpt_cent0to10_ynorm = (TH1D*) PbPb_jetpt_cent0to10->Clone("PbPb_jetpt_cent0to10_ynorm");
    TH1D *PbPb_jetphi_cent0to10_ynorm = (TH1D*) PbPb_jetphi_cent0to10->Clone("PbPb_jetphi_cent0to10_ynorm");
    TH1D *PbPb_jeteta_cent0to10_ynorm = (TH1D*) PbPb_jeteta_cent0to10->Clone("PbPb_jeteta_cent0to10_ynorm");
    TH1D *PbPb_mupt_cent0to10_ynorm = (TH1D*) PbPb_mupt_cent0to10->Clone("PbPb_mupt_cent0to10_ynorm");
    TH1D *PbPb_muphi_cent0to10_ynorm = (TH1D*) PbPb_muphi_cent0to10->Clone("PbPb_muphi_cent0to10_ynorm");
    TH1D *PbPb_mueta_cent0to10_ynorm = (TH1D*) PbPb_mueta_cent0to10->Clone("PbPb_mueta_cent0to10_ynorm");

    TH1D *PbPb_jetpt_cent10to30_ynorm = (TH1D*) PbPb_jetpt_cent10to30->Clone("PbPb_jetpt_cent10to30_ynorm");
    TH1D *PbPb_jetphi_cent10to30_ynorm = (TH1D*) PbPb_jetphi_cent10to30->Clone("PbPb_jetphi_cent10to30_ynorm");
    TH1D *PbPb_jeteta_cent10to30_ynorm = (TH1D*) PbPb_jeteta_cent10to30->Clone("PbPb_jeteta_cent10to30_ynorm");
    TH1D *PbPb_mupt_cent10to30_ynorm = (TH1D*) PbPb_mupt_cent10to30->Clone("PbPb_mupt_cent10to30_ynorm");
    TH1D *PbPb_muphi_cent10to30_ynorm = (TH1D*) PbPb_muphi_cent10to30->Clone("PbPb_muphi_cent10to30_ynorm");
    TH1D *PbPb_mueta_cent10to30_ynorm = (TH1D*) PbPb_mueta_cent10to30->Clone("PbPb_mueta_cent10to30_ynorm");

    TH1D *PbPb_jetpt_cent30to50_ynorm = (TH1D*) PbPb_jetpt_cent30to50->Clone("PbPb_jetpt_cent30to50_ynorm");
    TH1D *PbPb_jetphi_cent30to50_ynorm = (TH1D*) PbPb_jetphi_cent30to50->Clone("PbPb_jetphi_cent30to50_ynorm");
    TH1D *PbPb_jeteta_cent30to50_ynorm = (TH1D*) PbPb_jeteta_cent30to50->Clone("PbPb_jeteta_cent30to50_ynorm");
    TH1D *PbPb_mupt_cent30to50_ynorm = (TH1D*) PbPb_mupt_cent30to50->Clone("PbPb_mupt_cent30to50_ynorm");
    TH1D *PbPb_muphi_cent30to50_ynorm = (TH1D*) PbPb_muphi_cent30to50->Clone("PbPb_muphi_cent30to50_ynorm");
    TH1D *PbPb_mueta_cent30to50_ynorm = (TH1D*) PbPb_mueta_cent30to50->Clone("PbPb_mueta_cent30to50_ynorm");

    TH1D *PbPb_jetpt_cent50to90_ynorm = (TH1D*) PbPb_jetpt_cent50to90->Clone("PbPb_jetpt_cent50to90_ynorm");
    TH1D *PbPb_jetphi_cent50to90_ynorm = (TH1D*) PbPb_jetphi_cent50to90->Clone("PbPb_jetphi_cent50to90_ynorm");
    TH1D *PbPb_jeteta_cent50to90_ynorm = (TH1D*) PbPb_jeteta_cent50to90->Clone("PbPb_jeteta_cent50to90_ynorm");
    TH1D *PbPb_mupt_cent50to90_ynorm = (TH1D*) PbPb_mupt_cent50to90->Clone("PbPb_mupt_cent50to90_ynorm");
    TH1D *PbPb_muphi_cent50to90_ynorm = (TH1D*) PbPb_muphi_cent50to90->Clone("PbPb_muphi_cent50to90_ynorm");
    TH1D *PbPb_mueta_cent50to90_ynorm = (TH1D*) PbPb_mueta_cent50to90->Clone("PbPb_mueta_cent50to90_ynorm");

    pp_jetpt_ynorm->Scale(1.0/pp_jetpt_ynorm->Integral());
    pp_jetphi_ynorm->Scale(1.0/pp_jetphi_ynorm->Integral());
    pp_jeteta_ynorm->Scale(1.0/pp_jeteta_ynorm->Integral());
    pp_mupt_ynorm->Scale(1.0/pp_mupt_ynorm->Integral());
    pp_muphi_ynorm->Scale(1.0/pp_muphi_ynorm->Integral());
    pp_mueta_ynorm->Scale(1.0/pp_mueta_ynorm->Integral());

    PbPb_jetpt_cent0to10_ynorm->Scale(1.0/PbPb_jetpt_cent0to10_ynorm->Integral());
    PbPb_jetphi_cent0to10_ynorm->Scale(1.0/PbPb_jetphi_cent0to10_ynorm->Integral());
    PbPb_jeteta_cent0to10_ynorm->Scale(1.0/PbPb_jeteta_cent0to10_ynorm->Integral());
    PbPb_mupt_cent0to10_ynorm->Scale(1.0/PbPb_mupt_cent0to10_ynorm->Integral());
    PbPb_muphi_cent0to10_ynorm->Scale(1.0/PbPb_muphi_cent0to10_ynorm->Integral());
    PbPb_mueta_cent0to10_ynorm->Scale(1.0/PbPb_mueta_cent0to10_ynorm->Integral());

    PbPb_jetpt_cent10to30_ynorm->Scale(1.0/PbPb_jetpt_cent10to30_ynorm->Integral());
    PbPb_jetphi_cent10to30_ynorm->Scale(1.0/PbPb_jetphi_cent10to30_ynorm->Integral());
    PbPb_jeteta_cent10to30_ynorm->Scale(1.0/PbPb_jeteta_cent10to30_ynorm->Integral());
    PbPb_mupt_cent10to30_ynorm->Scale(1.0/PbPb_mupt_cent10to30_ynorm->Integral());
    PbPb_muphi_cent10to30_ynorm->Scale(1.0/PbPb_muphi_cent10to30_ynorm->Integral());
    PbPb_mueta_cent10to30_ynorm->Scale(1.0/PbPb_mueta_cent10to30_ynorm->Integral());

    PbPb_jetpt_cent30to50_ynorm->Scale(1.0/PbPb_jetpt_cent30to50_ynorm->Integral());
    PbPb_jetphi_cent30to50_ynorm->Scale(1.0/PbPb_jetphi_cent30to50_ynorm->Integral());
    PbPb_jeteta_cent30to50_ynorm->Scale(1.0/PbPb_jeteta_cent30to50_ynorm->Integral());
    PbPb_mupt_cent30to50_ynorm->Scale(1.0/PbPb_mupt_cent30to50_ynorm->Integral());
    PbPb_muphi_cent30to50_ynorm->Scale(1.0/PbPb_muphi_cent30to50_ynorm->Integral());
    PbPb_mueta_cent30to50_ynorm->Scale(1.0/PbPb_mueta_cent30to50_ynorm->Integral());

    PbPb_jetpt_cent50to90_ynorm->Scale(1.0/PbPb_jetpt_cent50to90_ynorm->Integral());
    PbPb_jetphi_cent50to90_ynorm->Scale(1.0/PbPb_jetphi_cent50to90_ynorm->Integral());
    PbPb_jeteta_cent50to90_ynorm->Scale(1.0/PbPb_jeteta_cent50to90_ynorm->Integral());
    PbPb_mupt_cent50to90_ynorm->Scale(1.0/PbPb_mupt_cent50to90_ynorm->Integral());
    PbPb_muphi_cent50to90_ynorm->Scale(1.0/PbPb_muphi_cent50to90_ynorm->Integral());
    PbPb_mueta_cent50to90_ynorm->Scale(1.0/PbPb_mueta_cent50to90_ynorm->Integral());

    // define rebin axes
     double pt_axis[16] = {50.0,60.0,70.0,80.0,90.0,100.0,120.0,140.0,160.0,180.0,200.0,250.0,300.0,350.0,400.0,500.0};
     double eta_axis[31] = {-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};
     double phi_axis[13] = {-150*TMath::Pi()/150,-125*TMath::Pi()/150,-100*TMath::Pi()/150,-75*TMath::Pi()/150,-50*TMath::Pi()/150,-25*TMath::Pi()/150,0.0,25*TMath::Pi()/150,
     50*TMath::Pi()/150,75*TMath::Pi()/150,100*TMath::Pi()/150,125*TMath::Pi()/150,150*TMath::Pi()/150};

    // create rebinned histograms
    TH1D *pp_jetpt_ynorm_rebin = (TH1D*) pp_jetpt_ynorm->Rebin(15,"pp_jetpt_ynorm_rebin",pt_axis);
    TH1D *pp_mupt_ynorm_rebin = (TH1D*) pp_mupt_ynorm->Rebin(15,"pp_mupt_ynorm_rebin",pt_axis);
    TH1D *PbPb_jetpt_cent0to10_ynorm_rebin = (TH1D*) PbPb_jetpt_cent0to10_ynorm->Rebin(15,"PbPb_jetpt_cent0to10_ynorm_rebin",pt_axis);
    TH1D *PbPb_mupt_cent0to10_ynorm_rebin = (TH1D*) PbPb_mupt_cent0to10_ynorm->Rebin(15,"PbPb_mupt_cent0to10_ynorm_rebin",pt_axis);
    TH1D *PbPb_jetpt_cent10to30_ynorm_rebin = (TH1D*) PbPb_jetpt_cent10to30_ynorm->Rebin(15,"PbPb_jetpt_cent10to30_ynorm_rebin",pt_axis);
    TH1D *PbPb_mupt_cent10to30_ynorm_rebin = (TH1D*) PbPb_mupt_cent10to30_ynorm->Rebin(15,"PbPb_mupt_cent10to30_ynorm_rebin",pt_axis);
    TH1D *PbPb_jetpt_cent30to50_ynorm_rebin = (TH1D*) PbPb_jetpt_cent30to50_ynorm->Rebin(15,"PbPb_jetpt_cent30to50_ynorm_rebin",pt_axis);
    TH1D *PbPb_mupt_cent30to50_ynorm_rebin = (TH1D*) PbPb_mupt_cent30to50_ynorm->Rebin(15,"PbPb_mupt_cent30to50_ynorm_rebin",pt_axis);
    TH1D *PbPb_jetpt_cent50to90_ynorm_rebin = (TH1D*) PbPb_jetpt_cent50to90_ynorm->Rebin(15,"PbPb_jetpt_cent50to90_ynorm_rebin",pt_axis);
    TH1D *PbPb_mupt_cent50to90_ynorm_rebin = (TH1D*) PbPb_mupt_cent50to90_ynorm->Rebin(15,"PbPb_mupt_cent50to90_ynorm_rebin",pt_axis);

    TH1D *pp_jeteta_ynorm_rebin = (TH1D*) pp_jeteta_ynorm->Rebin(30,"pp_jeteta_ynorm_rebin",eta_axis);
    TH1D *pp_mueta_ynorm_rebin = (TH1D*) pp_mueta_ynorm->Rebin(30,"pp_mueta_ynorm_rebin",eta_axis);
    TH1D *PbPb_jeteta_cent0to10_ynorm_rebin = (TH1D*) PbPb_jeteta_cent0to10_ynorm->Rebin(30,"PbPb_jeteta_cent0to10_ynorm_rebin",eta_axis);
    TH1D *PbPb_mueta_cent0to10_ynorm_rebin = (TH1D*) PbPb_mueta_cent0to10_ynorm->Rebin(30,"PbPb_mueta_cent0to10_ynorm_rebin",eta_axis);
    TH1D *PbPb_jeteta_cent10to30_ynorm_rebin = (TH1D*) PbPb_jeteta_cent10to30_ynorm->Rebin(30,"PbPb_jeteta_cent10to30_ynorm_rebin",eta_axis);
    TH1D *PbPb_mueta_cent10to30_ynorm_rebin = (TH1D*) PbPb_mueta_cent10to30_ynorm->Rebin(30,"PbPb_mueta_cent10to30_ynorm_rebin",eta_axis);
    TH1D *PbPb_jeteta_cent30to50_ynorm_rebin = (TH1D*) PbPb_jeteta_cent30to50_ynorm->Rebin(30,"PbPb_jeteta_cent30to50_ynorm_rebin",eta_axis);
    TH1D *PbPb_mueta_cent30to50_ynorm_rebin = (TH1D*) PbPb_mueta_cent30to50_ynorm->Rebin(30,"PbPb_mueta_cent30to50_ynorm_rebin",eta_axis);
    TH1D *PbPb_jeteta_cent50to90_ynorm_rebin = (TH1D*) PbPb_jeteta_cent50to90_ynorm->Rebin(30,"PbPb_jeteta_cent50to90_ynorm_rebin",eta_axis);
    TH1D *PbPb_mueta_cent50to90_ynorm_rebin = (TH1D*) PbPb_mueta_cent50to90_ynorm->Rebin(30,"PbPb_mueta_cent50to90_ynorm_rebin",eta_axis);

    TH1D *pp_jetphi_ynorm_rebin = (TH1D*) pp_jetphi_ynorm->Rebin(12,"pp_jetphi_ynorm_rebin",phi_axis);
    TH1D *pp_muphi_ynorm_rebin = (TH1D*) pp_muphi_ynorm->Rebin(12,"pp_muphi_ynorm_rebin",phi_axis);
    TH1D *PbPb_jetphi_cent0to10_ynorm_rebin = (TH1D*) PbPb_jetphi_cent0to10_ynorm->Rebin(12,"PbPb_jetphi_cent0to10_ynorm_rebin",phi_axis);
    TH1D *PbPb_muphi_cent0to10_ynorm_rebin = (TH1D*) PbPb_muphi_cent0to10_ynorm->Rebin(12,"PbPb_muphi_cent0to10_ynorm_rebin",phi_axis);
    TH1D *PbPb_jetphi_cent10to30_ynorm_rebin = (TH1D*) PbPb_jetphi_cent10to30_ynorm->Rebin(12,"PbPb_jetphi_cent10to30_ynorm_rebin",phi_axis);
    TH1D *PbPb_muphi_cent10to30_ynorm_rebin = (TH1D*) PbPb_muphi_cent10to30_ynorm->Rebin(12,"PbPb_muphi_cent10to30_ynorm_rebin",phi_axis);
    TH1D *PbPb_jetphi_cent30to50_ynorm_rebin = (TH1D*) PbPb_jetphi_cent30to50_ynorm->Rebin(12,"PbPb_jetphi_cent30to50_ynorm_rebin",phi_axis);
    TH1D *PbPb_muphi_cent30to50_ynorm_rebin = (TH1D*) PbPb_muphi_cent30to50_ynorm->Rebin(12,"PbPb_muphi_cent30to50_ynorm_rebin",phi_axis);
    TH1D *PbPb_jetphi_cent50to90_ynorm_rebin = (TH1D*) PbPb_jetphi_cent50to90_ynorm->Rebin(12,"PbPb_jetphi_cent50to90_ynorm_rebin",phi_axis);
    TH1D *PbPb_muphi_cent50to90_ynorm_rebin = (TH1D*) PbPb_muphi_cent50to90_ynorm->Rebin(12,"PbPb_muphi_cent50to90_ynorm_rebin",phi_axis);

    // normalize by bin width
     // pt loop

    TH1D *pp_jetpt_ynorm_rebin_xnorm = (TH1D*) pp_jetpt_ynorm_rebin->Clone("pp_jetpt_ynorm_rebin_xnorm");
    TH1D *pp_mupt_ynorm_rebin_xnorm = (TH1D*) pp_mupt_ynorm_rebin->Clone("pp_mupt_ynorm_rebin_xnorm");
    TH1D *PbPb_jetpt_cent0to10_ynorm_rebin_xnorm = (TH1D*) PbPb_jetpt_cent0to10_ynorm_rebin->Clone("PbPb_jetpt_cent0to10_ynorm_rebin_xnorm");
    TH1D *PbPb_mupt_cent0to10_ynorm_rebin_xnorm = (TH1D*) PbPb_mupt_cent0to10_ynorm_rebin->Clone("PbPb_mupt_cent0to10_ynorm_rebin_xnorm");
    TH1D *PbPb_jetpt_cent10to30_ynorm_rebin_xnorm = (TH1D*) PbPb_jetpt_cent10to30_ynorm_rebin->Clone("PbPb_jetpt_cent10to30_ynorm_rebin_xnorm");
    TH1D *PbPb_mupt_cent10to30_ynorm_rebin_xnorm = (TH1D*) PbPb_mupt_cent10to30_ynorm_rebin->Clone("PbPb_mupt_cent10to30_ynorm_rebin_xnorm");
    TH1D *PbPb_jetpt_cent30to50_ynorm_rebin_xnorm = (TH1D*) PbPb_jetpt_cent30to50_ynorm_rebin->Clone("PbPb_jetpt_cent30to50_ynorm_rebin_xnorm");
    TH1D *PbPb_mupt_cent30to50_ynorm_rebin_xnorm = (TH1D*) PbPb_mupt_cent30to50_ynorm_rebin->Clone("PbPb_mupt_cent30to50_ynorm_rebin_xnorm");
    TH1D *PbPb_jetpt_cent50to90_ynorm_rebin_xnorm = (TH1D*) PbPb_jetpt_cent50to90_ynorm_rebin->Clone("PbPb_jetpt_cent50to90_ynorm_rebin_xnorm");
    TH1D *PbPb_mupt_cent50to90_ynorm_rebin_xnorm = (TH1D*) PbPb_mupt_cent50to90_ynorm_rebin->Clone("PbPb_mupt_cent50to90_ynorm_rebin_xnorm");

    int N_pt = pp_jetpt_ynorm_rebin_xnorm->GetSize();
    for(int i=0;i<N_pt;i++){
        double w1 = pp_jetpt_ynorm_rebin_xnorm->GetBinWidth(i);
        double y1 = pp_jetpt_ynorm_rebin_xnorm->GetBinContent(i);
        double e1 = pp_jetpt_ynorm_rebin_xnorm->GetBinError(i);
        if(w1!=0){
            pp_jetpt_ynorm_rebin_xnorm->SetBinContent(i,y1/w1);
            pp_jetpt_ynorm_rebin_xnorm->SetBinError(i,e1/w1);
        }
        double w2 = pp_mupt_ynorm_rebin_xnorm->GetBinWidth(i);
        double y2 = pp_mupt_ynorm_rebin_xnorm->GetBinContent(i);
        double e2 = pp_mupt_ynorm_rebin_xnorm->GetBinError(i);
        if(w2!=0){
            pp_mupt_ynorm_rebin_xnorm->SetBinContent(i,y2/w2);
            pp_mupt_ynorm_rebin_xnorm->SetBinError(i,e2/w2);
        }
        double w3 = PbPb_jetpt_cent0to10_ynorm_rebin_xnorm->GetBinWidth(i);
        double y3 = PbPb_jetpt_cent0to10_ynorm_rebin_xnorm->GetBinContent(i);
        double e3 = PbPb_jetpt_cent0to10_ynorm_rebin_xnorm->GetBinError(i);
        if(w3!=0){
            PbPb_jetpt_cent0to10_ynorm_rebin_xnorm->SetBinContent(i,y3/w3);
            PbPb_jetpt_cent0to10_ynorm_rebin_xnorm->SetBinError(i,e3/w3);
        }
        double w4 = PbPb_mupt_cent0to10_ynorm_rebin_xnorm->GetBinWidth(i);
        double y4 = PbPb_mupt_cent0to10_ynorm_rebin_xnorm->GetBinContent(i);
        double e4 = PbPb_mupt_cent0to10_ynorm_rebin_xnorm->GetBinError(i);
        if(w4!=0){
            PbPb_mupt_cent0to10_ynorm_rebin_xnorm->SetBinContent(i,y4/w4);
            PbPb_mupt_cent0to10_ynorm_rebin_xnorm->SetBinError(i,e4/w4);
        }
        double w5 = PbPb_jetpt_cent10to30_ynorm_rebin_xnorm->GetBinWidth(i);
        double y5 = PbPb_jetpt_cent10to30_ynorm_rebin_xnorm->GetBinContent(i);
        double e5 = PbPb_jetpt_cent10to30_ynorm_rebin_xnorm->GetBinError(i);
        if(w5!=0){
            PbPb_jetpt_cent10to30_ynorm_rebin_xnorm->SetBinContent(i,y5/w5);
            PbPb_jetpt_cent10to30_ynorm_rebin_xnorm->SetBinError(i,e5/w5);
        }
        double w6 = PbPb_mupt_cent10to30_ynorm_rebin_xnorm->GetBinWidth(i);
        double y6 = PbPb_mupt_cent10to30_ynorm_rebin_xnorm->GetBinContent(i);
        double e6 = PbPb_mupt_cent10to30_ynorm_rebin_xnorm->GetBinError(i);
        if(w6!=0){
            PbPb_mupt_cent10to30_ynorm_rebin_xnorm->SetBinContent(i,y6/w6);
            PbPb_mupt_cent10to30_ynorm_rebin_xnorm->SetBinError(i,e6/w6);
        }
        double w7 = PbPb_jetpt_cent30to50_ynorm_rebin_xnorm->GetBinWidth(i);
        double y7 = PbPb_jetpt_cent30to50_ynorm_rebin_xnorm->GetBinContent(i);
        double e7 = PbPb_jetpt_cent30to50_ynorm_rebin_xnorm->GetBinError(i);
        if(w7!=0){
            PbPb_jetpt_cent30to50_ynorm_rebin_xnorm->SetBinContent(i,y7/w7);
            PbPb_jetpt_cent30to50_ynorm_rebin_xnorm->SetBinError(i,e7/w7);
        }
        double w8 = PbPb_mupt_cent30to50_ynorm_rebin_xnorm->GetBinWidth(i);
        double y8 = PbPb_mupt_cent30to50_ynorm_rebin_xnorm->GetBinContent(i);
        double e8 = PbPb_mupt_cent30to50_ynorm_rebin_xnorm->GetBinError(i);
        if(w8!=0){
            PbPb_mupt_cent30to50_ynorm_rebin_xnorm->SetBinContent(i,y8/w8);
            PbPb_mupt_cent30to50_ynorm_rebin_xnorm->SetBinError(i,e8/w8);
        }
        double w9 = PbPb_jetpt_cent50to90_ynorm_rebin_xnorm->GetBinWidth(i);
        double y9 = PbPb_jetpt_cent50to90_ynorm_rebin_xnorm->GetBinContent(i);
        double e9 = PbPb_jetpt_cent50to90_ynorm_rebin_xnorm->GetBinError(i);
        if(w9!=0){
            PbPb_jetpt_cent50to90_ynorm_rebin_xnorm->SetBinContent(i,y9/w9);
            PbPb_jetpt_cent50to90_ynorm_rebin_xnorm->SetBinError(i,e9/w9);
        }
        double w10 = PbPb_mupt_cent50to90_ynorm_rebin_xnorm->GetBinWidth(i);
        double y10 = PbPb_mupt_cent50to90_ynorm_rebin_xnorm->GetBinContent(i);
        double e10 = PbPb_mupt_cent50to90_ynorm_rebin_xnorm->GetBinError(i);
        if(w10!=0){
            PbPb_mupt_cent50to90_ynorm_rebin_xnorm->SetBinContent(i,y10/w10);
            PbPb_mupt_cent50to90_ynorm_rebin_xnorm->SetBinError(i,e10/w10);
        }
    }

    // eta loop

    TH1D *pp_jeteta_ynorm_rebin_xnorm = (TH1D*) pp_jeteta_ynorm_rebin->Clone("pp_jeteta_ynorm_rebin_xnorm");
    TH1D *pp_mueta_ynorm_rebin_xnorm = (TH1D*) pp_mueta_ynorm_rebin->Clone("pp_mueta_ynorm_rebin_xnorm");
    TH1D *PbPb_jeteta_cent0to10_ynorm_rebin_xnorm = (TH1D*) PbPb_jeteta_cent0to10_ynorm_rebin->Clone("PbPb_jeteta_cent0to10_ynorm_rebin_xnorm");
    TH1D *PbPb_mueta_cent0to10_ynorm_rebin_xnorm = (TH1D*) PbPb_mueta_cent0to10_ynorm_rebin->Clone("PbPb_mueta_cent0to10_ynorm_rebin_xnorm");
    TH1D *PbPb_jeteta_cent10to30_ynorm_rebin_xnorm = (TH1D*) PbPb_jeteta_cent10to30_ynorm_rebin->Clone("PbPb_jeteta_cent10to30_ynorm_rebin_xnorm");
    TH1D *PbPb_mueta_cent10to30_ynorm_rebin_xnorm = (TH1D*) PbPb_mueta_cent10to30_ynorm_rebin->Clone("PbPb_mueta_cent10to30_ynorm_rebin_xnorm");
    TH1D *PbPb_jeteta_cent30to50_ynorm_rebin_xnorm = (TH1D*) PbPb_jeteta_cent30to50_ynorm_rebin->Clone("PbPb_jeteta_cent30to50_ynorm_rebin_xnorm");
    TH1D *PbPb_mueta_cent30to50_ynorm_rebin_xnorm = (TH1D*) PbPb_mueta_cent30to50_ynorm_rebin->Clone("PbPb_mueta_cent30to50_ynorm_rebin_xnorm");
    TH1D *PbPb_jeteta_cent50to90_ynorm_rebin_xnorm = (TH1D*) PbPb_jeteta_cent50to90_ynorm_rebin->Clone("PbPb_jeteta_cent50to90_ynorm_rebin_xnorm");
    TH1D *PbPb_mueta_cent50to90_ynorm_rebin_xnorm = (TH1D*) PbPb_mueta_cent50to90_ynorm_rebin->Clone("PbPb_mueta_cent50to90_ynorm_rebin_xnorm");

    int N_eta = pp_jeteta_ynorm_rebin_xnorm->GetSize();
    for(int i=0;i<N_eta;i++){
        double w1 = pp_jeteta_ynorm_rebin_xnorm->GetBinWidth(i);
        double y1 = pp_jeteta_ynorm_rebin_xnorm->GetBinContent(i);
        double e1 = pp_jeteta_ynorm_rebin_xnorm->GetBinError(i);
        if(w1!=0){
            pp_jeteta_ynorm_rebin_xnorm->SetBinContent(i,y1/w1);
            pp_jeteta_ynorm_rebin_xnorm->SetBinError(i,e1/w1);
        }
        double w2 = pp_mueta_ynorm_rebin_xnorm->GetBinWidth(i);
        double y2 = pp_mueta_ynorm_rebin_xnorm->GetBinContent(i);
        double e2 = pp_mueta_ynorm_rebin_xnorm->GetBinError(i);
        if(w2!=0){
            pp_mueta_ynorm_rebin_xnorm->SetBinContent(i,y2/w2);
            pp_mueta_ynorm_rebin_xnorm->SetBinError(i,e2/w2);
        }
        double w3 = PbPb_jeteta_cent0to10_ynorm_rebin_xnorm->GetBinWidth(i);
        double y3 = PbPb_jeteta_cent0to10_ynorm_rebin_xnorm->GetBinContent(i);
        double e3 = PbPb_jeteta_cent0to10_ynorm_rebin_xnorm->GetBinError(i);
        if(w3!=0){
            PbPb_jeteta_cent0to10_ynorm_rebin_xnorm->SetBinContent(i,y3/w3);
            PbPb_jeteta_cent0to10_ynorm_rebin_xnorm->SetBinError(i,e3/w3);
        }
        double w4 = PbPb_mueta_cent0to10_ynorm_rebin_xnorm->GetBinWidth(i);
        double y4 = PbPb_mueta_cent0to10_ynorm_rebin_xnorm->GetBinContent(i);
        double e4 = PbPb_mueta_cent0to10_ynorm_rebin_xnorm->GetBinError(i);
        if(w4!=0){
            PbPb_mueta_cent0to10_ynorm_rebin_xnorm->SetBinContent(i,y4/w4);
            PbPb_mueta_cent0to10_ynorm_rebin_xnorm->SetBinError(i,e4/w4);
        }
        double w5 = PbPb_jeteta_cent10to30_ynorm_rebin_xnorm->GetBinWidth(i);
        double y5 = PbPb_jeteta_cent10to30_ynorm_rebin_xnorm->GetBinContent(i);
        double e5 = PbPb_jeteta_cent10to30_ynorm_rebin_xnorm->GetBinError(i);
        if(w5!=0){
            PbPb_jeteta_cent10to30_ynorm_rebin_xnorm->SetBinContent(i,y5/w5);
            PbPb_jeteta_cent10to30_ynorm_rebin_xnorm->SetBinError(i,e5/w5);
        }
        double w6 = PbPb_mueta_cent10to30_ynorm_rebin_xnorm->GetBinWidth(i);
        double y6 = PbPb_mueta_cent10to30_ynorm_rebin_xnorm->GetBinContent(i);
        double e6 = PbPb_mueta_cent10to30_ynorm_rebin_xnorm->GetBinError(i);
        if(w6!=0){
            PbPb_mueta_cent10to30_ynorm_rebin_xnorm->SetBinContent(i,y6/w6);
            PbPb_mueta_cent10to30_ynorm_rebin_xnorm->SetBinError(i,e6/w6);
        }
        double w7 = PbPb_jeteta_cent30to50_ynorm_rebin_xnorm->GetBinWidth(i);
        double y7 = PbPb_jeteta_cent30to50_ynorm_rebin_xnorm->GetBinContent(i);
        double e7 = PbPb_jeteta_cent30to50_ynorm_rebin_xnorm->GetBinError(i);
        if(w7!=0){
            PbPb_jeteta_cent30to50_ynorm_rebin_xnorm->SetBinContent(i,y7/w7);
            PbPb_jeteta_cent30to50_ynorm_rebin_xnorm->SetBinError(i,e7/w7);
        }
        double w8 = PbPb_mueta_cent30to50_ynorm_rebin_xnorm->GetBinWidth(i);
        double y8 = PbPb_mueta_cent30to50_ynorm_rebin_xnorm->GetBinContent(i);
        double e8 = PbPb_mueta_cent30to50_ynorm_rebin_xnorm->GetBinError(i);
        if(w8!=0){
            PbPb_mueta_cent30to50_ynorm_rebin_xnorm->SetBinContent(i,y8/w8);
            PbPb_mueta_cent30to50_ynorm_rebin_xnorm->SetBinError(i,e8/w8);
        }
        double w9 = PbPb_jeteta_cent50to90_ynorm_rebin_xnorm->GetBinWidth(i);
        double y9 = PbPb_jeteta_cent50to90_ynorm_rebin_xnorm->GetBinContent(i);
        double e9 = PbPb_jeteta_cent50to90_ynorm_rebin_xnorm->GetBinError(i);
        if(w9!=0){
            PbPb_jeteta_cent50to90_ynorm_rebin_xnorm->SetBinContent(i,y9/w9);
            PbPb_jeteta_cent50to90_ynorm_rebin_xnorm->SetBinError(i,e9/w9);
        }
        double w10 = PbPb_mueta_cent50to90_ynorm_rebin_xnorm->GetBinWidth(i);
        double y10 = PbPb_mueta_cent50to90_ynorm_rebin_xnorm->GetBinContent(i);
        double e10 = PbPb_mueta_cent50to90_ynorm_rebin_xnorm->GetBinError(i);
        if(w10!=0){
            PbPb_mueta_cent50to90_ynorm_rebin_xnorm->SetBinContent(i,y10/w10);
            PbPb_mueta_cent50to90_ynorm_rebin_xnorm->SetBinError(i,e10/w10);
        }
    }

     // phi loop

    TH1D *pp_jetphi_ynorm_rebin_xnorm = (TH1D*) pp_jetphi_ynorm_rebin->Clone("pp_jetphi_ynorm_rebin_xnorm");
    TH1D *pp_muphi_ynorm_rebin_xnorm = (TH1D*) pp_muphi_ynorm_rebin->Clone("pp_muphi_ynorm_rebin_xnorm");
    TH1D *PbPb_jetphi_cent0to10_ynorm_rebin_xnorm = (TH1D*) PbPb_jetphi_cent0to10_ynorm_rebin->Clone("PbPb_jetphi_cent0to10_ynorm_rebin_xnorm");
    TH1D *PbPb_muphi_cent0to10_ynorm_rebin_xnorm = (TH1D*) PbPb_muphi_cent0to10_ynorm_rebin->Clone("PbPb_muphi_cent0to10_ynorm_rebin_xnorm");
    TH1D *PbPb_jetphi_cent10to30_ynorm_rebin_xnorm = (TH1D*) PbPb_jetphi_cent10to30_ynorm_rebin->Clone("PbPb_jetphi_cent10to30_ynorm_rebin_xnorm");
    TH1D *PbPb_muphi_cent10to30_ynorm_rebin_xnorm = (TH1D*) PbPb_muphi_cent10to30_ynorm_rebin->Clone("PbPb_muphi_cent10to30_ynorm_rebin_xnorm");
    TH1D *PbPb_jetphi_cent30to50_ynorm_rebin_xnorm = (TH1D*) PbPb_jetphi_cent30to50_ynorm_rebin->Clone("PbPb_jetphi_cent30to50_ynorm_rebin_xnorm");
    TH1D *PbPb_muphi_cent30to50_ynorm_rebin_xnorm = (TH1D*) PbPb_muphi_cent30to50_ynorm_rebin->Clone("PbPb_muphi_cent30to50_ynorm_rebin_xnorm");
    TH1D *PbPb_jetphi_cent50to90_ynorm_rebin_xnorm = (TH1D*) PbPb_jetphi_cent50to90_ynorm_rebin->Clone("PbPb_jetphi_cent50to90_ynorm_rebin_xnorm");
    TH1D *PbPb_muphi_cent50to90_ynorm_rebin_xnorm = (TH1D*) PbPb_muphi_cent50to90_ynorm_rebin->Clone("PbPb_muphi_cent50to90_ynorm_rebin_xnorm");

    int N_phi = pp_jetphi_ynorm_rebin_xnorm->GetSize();
    for(int i=0;i<N_phi;i++){
        double w1 = pp_jetphi_ynorm_rebin_xnorm->GetBinWidth(i);
        double y1 = pp_jetphi_ynorm_rebin_xnorm->GetBinContent(i);
        double e1 = pp_jetphi_ynorm_rebin_xnorm->GetBinError(i);
        if(w1!=0){
            pp_jetphi_ynorm_rebin_xnorm->SetBinContent(i,y1/w1);
            pp_jetphi_ynorm_rebin_xnorm->SetBinError(i,e1/w1);
        }
        double w2 = pp_muphi_ynorm_rebin_xnorm->GetBinWidth(i);
        double y2 = pp_muphi_ynorm_rebin_xnorm->GetBinContent(i);
        double e2 = pp_muphi_ynorm_rebin_xnorm->GetBinError(i);
        if(w2!=0){
            pp_muphi_ynorm_rebin_xnorm->SetBinContent(i,y2/w2);
            pp_muphi_ynorm_rebin_xnorm->SetBinError(i,e2/w2);
        }
        double w3 = PbPb_jetphi_cent0to10_ynorm_rebin_xnorm->GetBinWidth(i);
        double y3 = PbPb_jetphi_cent0to10_ynorm_rebin_xnorm->GetBinContent(i);
        double e3 = PbPb_jetphi_cent0to10_ynorm_rebin_xnorm->GetBinError(i);
        if(w3!=0){
            PbPb_jetphi_cent0to10_ynorm_rebin_xnorm->SetBinContent(i,y3/w3);
            PbPb_jetphi_cent0to10_ynorm_rebin_xnorm->SetBinError(i,e3/w3);
        }
        double w4 = PbPb_muphi_cent0to10_ynorm_rebin_xnorm->GetBinWidth(i);
        double y4 = PbPb_muphi_cent0to10_ynorm_rebin_xnorm->GetBinContent(i);
        double e4 = PbPb_muphi_cent0to10_ynorm_rebin_xnorm->GetBinError(i);
        if(w4!=0){
            PbPb_muphi_cent0to10_ynorm_rebin_xnorm->SetBinContent(i,y4/w4);
            PbPb_muphi_cent0to10_ynorm_rebin_xnorm->SetBinError(i,e4/w4);
        }
        double w5 = PbPb_jetphi_cent10to30_ynorm_rebin_xnorm->GetBinWidth(i);
        double y5 = PbPb_jetphi_cent10to30_ynorm_rebin_xnorm->GetBinContent(i);
        double e5 = PbPb_jetphi_cent10to30_ynorm_rebin_xnorm->GetBinError(i);
        if(w5!=0){
            PbPb_jetphi_cent10to30_ynorm_rebin_xnorm->SetBinContent(i,y5/w5);
            PbPb_jetphi_cent10to30_ynorm_rebin_xnorm->SetBinError(i,e5/w5);
        }
        double w6 = PbPb_muphi_cent10to30_ynorm_rebin_xnorm->GetBinWidth(i);
        double y6 = PbPb_muphi_cent10to30_ynorm_rebin_xnorm->GetBinContent(i);
        double e6 = PbPb_muphi_cent10to30_ynorm_rebin_xnorm->GetBinError(i);
        if(w6!=0){
            PbPb_muphi_cent10to30_ynorm_rebin_xnorm->SetBinContent(i,y6/w6);
            PbPb_muphi_cent10to30_ynorm_rebin_xnorm->SetBinError(i,e6/w6);
        }
        double w7 = PbPb_jetphi_cent30to50_ynorm_rebin_xnorm->GetBinWidth(i);
        double y7 = PbPb_jetphi_cent30to50_ynorm_rebin_xnorm->GetBinContent(i);
        double e7 = PbPb_jetphi_cent30to50_ynorm_rebin_xnorm->GetBinError(i);
        if(w7!=0){
            PbPb_jetphi_cent30to50_ynorm_rebin_xnorm->SetBinContent(i,y7/w7);
            PbPb_jetphi_cent30to50_ynorm_rebin_xnorm->SetBinError(i,e7/w7);
        }
        double w8 = PbPb_muphi_cent30to50_ynorm_rebin_xnorm->GetBinWidth(i);
        double y8 = PbPb_muphi_cent30to50_ynorm_rebin_xnorm->GetBinContent(i);
        double e8 = PbPb_muphi_cent30to50_ynorm_rebin_xnorm->GetBinError(i);
        if(w8!=0){
            PbPb_muphi_cent30to50_ynorm_rebin_xnorm->SetBinContent(i,y8/w8);
            PbPb_muphi_cent30to50_ynorm_rebin_xnorm->SetBinError(i,e8/w8);
        }
        double w9 = PbPb_jetphi_cent50to90_ynorm_rebin_xnorm->GetBinWidth(i);
        double y9 = PbPb_jetphi_cent50to90_ynorm_rebin_xnorm->GetBinContent(i);
        double e9 = PbPb_jetphi_cent50to90_ynorm_rebin_xnorm->GetBinError(i);
        if(w9!=0){
            PbPb_jetphi_cent50to90_ynorm_rebin_xnorm->SetBinContent(i,y9/w9);
            PbPb_jetphi_cent50to90_ynorm_rebin_xnorm->SetBinError(i,e9/w9);
        }
        double w10 = PbPb_muphi_cent50to90_ynorm_rebin_xnorm->GetBinWidth(i);
        double y10 = PbPb_muphi_cent50to90_ynorm_rebin_xnorm->GetBinContent(i);
        double e10 = PbPb_muphi_cent50to90_ynorm_rebin_xnorm->GetBinError(i);
        if(w10!=0){
            PbPb_muphi_cent50to90_ynorm_rebin_xnorm->SetBinContent(i,y10/w10);
            PbPb_muphi_cent50to90_ynorm_rebin_xnorm->SetBinError(i,e10/w10);
        }
    }




      TCanvas *c1 = new TCanvas("c2", "p_{T} spectra", 500, 500); 
      
      c1->cd();
      
      //c2->SetLogy();
      TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.3, 1., 1.);
      //pad1->SetGridx();         // Vertical grid
      pad2->SetLeftMargin(0.15);
      pad2->SetLogy();
      pad2->Draw();             // Draw the upper pad: pad1
      pad2->cd();
      

      // jet pT 
    
    
      Double_t msize = 0.9;

      
      pp_jetpt_ynorm_rebin_xnorm->SetMarkerStyle(8);
      pp_jetpt_ynorm_rebin_xnorm->SetMarkerColor(kBlack);
      pp_jetpt_ynorm_rebin_xnorm->SetMarkerSize(msize);
      pp_jetpt_ynorm_rebin_xnorm->SetFillColorAlpha(kBlack,0.7);
      pp_jetpt_ynorm_rebin_xnorm->SetStats(0);
      pp_jetpt_ynorm_rebin_xnorm->SetTitle("");

      PbPb_jetpt_cent50to90_ynorm_rebin_xnorm->SetMarkerStyle(8);
      PbPb_jetpt_cent50to90_ynorm_rebin_xnorm->SetMarkerColor(kBlue);
      PbPb_jetpt_cent50to90_ynorm_rebin_xnorm->SetMarkerSize(msize);
      PbPb_jetpt_cent50to90_ynorm_rebin_xnorm->SetFillColorAlpha(kBlue,0.7);
      PbPb_jetpt_cent50to90_ynorm_rebin_xnorm->SetStats(0);
      PbPb_jetpt_cent50to90_ynorm_rebin_xnorm->SetTitle("");

      PbPb_jetpt_cent30to50_ynorm_rebin_xnorm->SetMarkerStyle(8);
      PbPb_jetpt_cent30to50_ynorm_rebin_xnorm->SetMarkerColor(kGreen+2);
      PbPb_jetpt_cent30to50_ynorm_rebin_xnorm->SetMarkerSize(msize);
      PbPb_jetpt_cent30to50_ynorm_rebin_xnorm->SetFillColorAlpha(kGreen+2,0.7);
      PbPb_jetpt_cent30to50_ynorm_rebin_xnorm->SetStats(0);
      PbPb_jetpt_cent30to50_ynorm_rebin_xnorm->SetTitle("");

      PbPb_jetpt_cent10to30_ynorm_rebin_xnorm->SetMarkerStyle(8);
      PbPb_jetpt_cent10to30_ynorm_rebin_xnorm->SetMarkerColor(kRed);
      PbPb_jetpt_cent10to30_ynorm_rebin_xnorm->SetMarkerSize(msize);
      PbPb_jetpt_cent10to30_ynorm_rebin_xnorm->SetFillColorAlpha(kRed,0.7);
      PbPb_jetpt_cent10to30_ynorm_rebin_xnorm->SetStats(0);
      PbPb_jetpt_cent10to30_ynorm_rebin_xnorm->SetTitle(""); 

      PbPb_jetpt_cent0to10_ynorm_rebin_xnorm->SetMarkerStyle(8);
      PbPb_jetpt_cent0to10_ynorm_rebin_xnorm->SetMarkerColor(kBlue-10);
      PbPb_jetpt_cent0to10_ynorm_rebin_xnorm->SetMarkerSize(msize);
      PbPb_jetpt_cent0to10_ynorm_rebin_xnorm->SetFillColorAlpha(kBlue-10,0.7);
      PbPb_jetpt_cent0to10_ynorm_rebin_xnorm->SetStats(0);
      PbPb_jetpt_cent0to10_ynorm_rebin_xnorm->SetTitle("");

 






      //clayton_pt_rebin_xnorm_fitscaled->Draw("e2");
      pp_jetpt_ynorm_rebin_xnorm->Draw("e2");
      PbPb_jetpt_cent50to90_ynorm_rebin_xnorm->Draw("e2 same");
      PbPb_jetpt_cent30to50_ynorm_rebin_xnorm->Draw("e2 same");
      PbPb_jetpt_cent10to30_ynorm_rebin_xnorm->Draw("e2 same");
      PbPb_jetpt_cent0to10_ynorm_rebin_xnorm->Draw("e2 same");
      
      
      
      

      auto legend = new TLegend(0.7,0.5,0.85,0.7);
      legend->AddEntry(pp_jetpt_ynorm_rebin_xnorm,"pp","p");
      legend->AddEntry(PbPb_jetpt_cent50to90_ynorm_rebin_xnorm,"50-90 %","p");
      legend->AddEntry(PbPb_jetpt_cent30to50_ynorm_rebin_xnorm,"30-50 %","p");
      legend->AddEntry(PbPb_jetpt_cent10to30_ynorm_rebin_xnorm,"10-30 %","p");
      legend->AddEntry(PbPb_jetpt_cent0to10_ynorm_rebin_xnorm,"0-10 %","p");
      legend->SetBorderSize(0);
      legend->Draw();


      pp_jetpt_ynorm_rebin_xnorm->GetYaxis()->SetTitle("1/N^{jet} dN^{jet}/dp_{T}");
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
    TH1D *r1 = (TH1D*)PbPb_jetpt_cent0to10_ynorm_rebin_xnorm->Clone("r1");
    TH1D *r2 = (TH1D*)PbPb_jetpt_cent10to30_ynorm_rebin_xnorm->Clone("r2");
    TH1D *r3 = (TH1D*)PbPb_jetpt_cent30to50_ynorm_rebin_xnorm->Clone("r3");
    TH1D *r4 = (TH1D*)PbPb_jetpt_cent50to90_ynorm_rebin_xnorm->Clone("r4");
    r1->Divide(pp_jetpt_ynorm_rebin_xnorm);
    r2->Divide(pp_jetpt_ynorm_rebin_xnorm);
    r3->Divide(pp_jetpt_ynorm_rebin_xnorm);
    r4->Divide(pp_jetpt_ynorm_rebin_xnorm);
    r1->SetMarkerStyle(24);
    r2->SetMarkerStyle(24);
    r3->SetMarkerStyle(24);
    r4->SetMarkerStyle(24);
    r1->SetStats(0);
    r1->SetMinimum(0.0);
    r1->SetMaximum(4.0);
    r1->GetXaxis()->SetTitle("jet p_{T} [GeV/c]");
    r1->GetXaxis()->SetTitleSize(0.10);
    r1->GetXaxis()->SetLabelSize(0.08);
    r1->GetYaxis()->SetTitle("Centrality region / pp");
    r1->GetYaxis()->SetTitleSize(0.10);
    r1->GetYaxis()->SetLabelSize(0.08);
    r1->SetTitle("");
    r1->Draw("e2");
    r2->Draw("e2 same");
    r3->Draw("e2 same");
    r4->Draw("e2 same");

    c1->SaveAs("/home/clayton/Analysis/code/normalizeAndCompareSpectra/figures/DataPPvsPbPbCentralityRegionsPt.pdf");



} // end program