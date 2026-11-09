#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include "TFile.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TList.h"
#include "TString.h"
#include "TH1F.h"
#include "TObject.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include <cmath>
#include <cstring>
#include <limits>
#include "BarlowTest.h"
#include "ReadBarlowTest.h"
#include "TProfile2D.h"

using namespace std;

// MC (LHC25f3_gp): 544886, 544868, 544122, 544116
// real (LHC23zzh_pass5): 544124, 544123, 544122, 544121, 544116, 544098, 544095, 544091
// with phi and pt weight: 544116, 544122
// with phi weight only: 544091, 544095, 544098, 544121, 544123, 544124
// with pt weight only: 544868, 544886

// observables vs. nTracks
// MCdata
// - <<2>> n = 1, 2, 3, 4
// - <<3>> (n1, n2, n3) = (2, 1, 1), (3, 2, 1), (4, 2, 2)
// - <<4>> n = 2, 3, 4
// - <<4>> (n1, n2, n3, n4) = (3, 2, 3, 2), (4, 2, 4, 2)

enum TwoPCorrName {eV11 = 0, eV22, eV33, eV44, n2};
enum ThreePCorrName {eV211 = 0, eV321, eV422, n3};
enum FourPCorrName {eV3232 = 0, eV4242, n4};
// size_t n2 = 4;
// size_t n3 = 3;
// size_t n4 = 5;
const size_t nbin = 20;

double TwoPCorr[n2][nbin];
double ThreePCorr[n3][nbin];
double FourPCorr[n4][nbin];
double nTracks[nbin];


int CalculateData_kpcorrVSntrack() {
    TFile* analysisFile = TFile::Open("AnalysisResults.root", "read");

    TString dirName = "multiparticle-cumulants";

    TDirectory* dir = dynamic_cast<TDirectory*>(analysisFile->Get(dirName.Data()));


    printf(" Open dir %s\n", dirName.Data());

    TList* baseList = nullptr;
    TList* histogramList = nullptr;
    dir->GetObject("Default list name", baseList);
    if(!baseList) {
        printf(" Failed to find list \"Default list name\"\n");
        return 0;
    }
    histogramList = dynamic_cast<TList*>(baseList->FindObject("MultiparticleCorrelationProfiles"));
    if(!histogramList) {
        printf(" Failed to find list \"MultiparticleCorrelationProfiles\"\n");
        return 0;
    }

    TList* runNumberList = dynamic_cast<TList*>(histogramList->FindObject("mcc_run544886"));

    if(!runNumberList) {
        printf(" Failed to find list \"mcc_run544886\"\n");
        return 0;
    }

    TProfile2D* prof2 = dynamic_cast<TProfile2D*>(runNumberList->FindObject("prof2NumTrackAfterCut"));
    TProfile2D* prof3 = dynamic_cast<TProfile2D*>(runNumberList->FindObject("prof3NumTrackAfterCut"));
    TProfile2D* prof4 = dynamic_cast<TProfile2D*>(runNumberList->FindObject("prof4NumTrackAfterCut"));

    for(int i = 0; i < nbin; i++) {
        TwoPCorr[eV11][i] = prof2->GetBinContent(prof2->GetXaxis()->FindBin("<v1^2>"), i + 1);
        TwoPCorr[eV22][i] = prof2->GetBinContent(prof2->GetXaxis()->FindBin("<v2^2>"), i + 1);
        TwoPCorr[eV33][i] = prof2->GetBinContent(prof2->GetXaxis()->FindBin("<v3^2>"), i + 1);
        TwoPCorr[eV44][i] = prof2->GetBinContent(prof2->GetXaxis()->FindBin("<v4^2>"), i + 1);
        ThreePCorr[eV211][i] = prof3->GetBinContent(prof3->GetXaxis()->FindBin("<v2v1v1>"), i + 1);
        ThreePCorr[eV321][i] = prof3->GetBinContent(prof3->GetXaxis()->FindBin("<v3v2v1>"), i + 1);
        ThreePCorr[eV422][i] = prof3->GetBinContent(prof3->GetXaxis()->FindBin("<v4v2v2>"), i + 1);
        FourPCorr[eV3232][i] = prof4->GetBinContent(prof4->GetXaxis()->FindBin("<v2^2v3^2>"), i + 1);
        FourPCorr[eV4242][i] = prof4->GetBinContent(prof4->GetXaxis()->FindBin("<v2^2v4^2>"), i + 1);
        nTracks[i] = prof2->GetYaxis()->GetBinCenter(i + 1);
    }


    TGraph* g2[n2];
    TGraph* g3[n3];
    TGraph* g4[n4];
    TCanvas* c = new TCanvas("c", "c", 800, 600);

    for(int i = 0; i < n2; i++) {
        g2[i] = new TGraph(nbin, nTracks, TwoPCorr[i]);
        g2[i]->Draw();
        c->Print(Form("result/temp_2p_%d.png", i));
        c->Clear();
    }
    for(int i = 0; i < n3; i++) {
        g3[i] = new TGraph(nbin, nTracks, ThreePCorr[i]);
        g3[i]->Draw();
        c->Print(Form("result/temp_3p_%d.png", i));
        c->Clear();
    }
    for(int i = 0; i < n4; i++) {
        g4[i] = new TGraph(nbin, nTracks, FourPCorr[i]);
        g4[i]->Draw();
        c->Print(Form("result/temp_4p_%d.png", i));
        c->Clear();
    }

    return 1;
}