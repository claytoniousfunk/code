#include "myProcesses/hiforest/plugin/eventMap_hiForest.h"
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

void displayData(){


TString d;
double var = 9.05;
d = "A string";
d += var;
d += " in some units";
TLatex *l = new TLatex(0.2,0.7,d.Data());
// and so on

TCanvas *c1 = new TCanvas("c1","canvas",500,500);
//c1->cd();
//TPad *pad1 = new TPad("pad1", "pad1", 0., 0., 1., 1.);
//pad1->cd();
//pad1->Draw();
l->SetTextSize(0.02);
l->DrawLatex(0.1,0.5,d.Data());



}