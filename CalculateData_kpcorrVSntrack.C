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

enum enTwoPCorrName {eV11 = 0, eV22, eV33, eV44, n2};
enum enThreePCorrName {eV211 = 0, eV321, eV422, n3};
enum enFourPCorrName {eV2222 = 0, eV3333, eV4444, eV3232, eV4242, n4};

const int TwoPCorrName[n2][2] = {{1,-1}, {2,-2}, {3,-3}, {4,-4}};
const int ThreePCorrName[n3][3] = {{2,-1,-1}, {3,-2,-1}, {4,-2,-2}};
const int FourPCorrName[n4][4] = {{2,2,-2,-2}, {3,3,-3,-3}, {4,4,-4,-4}, {3,2,-3,-2}, {4,2,-4,-2}};

const size_t nRun = 4;
const int runNumber[nRun] = {544886, 544868, 544122, 544116};

// size_t n2 = 4;
// size_t n3 = 3;
// size_t n4 = 5;
const size_t nbin = 20;

double TwoPCorrByRun[nRun][n2][nbin];
double ThreePCorrByRun[nRun][n3][nbin];
double FourPCorrByRun[nRun][n4][nbin];

double TwoPCorr[n2][nbin];
double ThreePCorr[n3][nbin];
double FourPCorr[n4][nbin];

double nTracks[nbin];

const char* dataDir = "/scratch3/go52dab/O2tutorial/test20260922_analysis/";
const double nan = std::numeric_limits<double>::quiet_NaN();

int CalculateData_kpcorrVSntrack() {

    // fill array

    for (int i = 0; i < nRun; i++) {


        TFile* analysisFile = TFile::Open(Form("%sMCdata/%d/AnalysisResults.root", dataDir, runNumber[i]), "read"); // Form("%sMCdata/544886/AnalysisResults.root", dataDir)

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

        TList* runNumberList = dynamic_cast<TList*>(histogramList->FindObject(Form("mcc_run%d", runNumber[i])));

        if(!runNumberList) {
            printf(Form(" Failed to find list \"mcc_run%d\"\n", runNumber[i]));
            return 0;
        }

        TProfile2D* prof2 = dynamic_cast<TProfile2D*>(runNumberList->FindObject("prof2NumTrackAfterCut"));
        TProfile2D* prof3 = dynamic_cast<TProfile2D*>(runNumberList->FindObject("prof3NumTrackAfterCut"));
        TProfile2D* prof4 = dynamic_cast<TProfile2D*>(runNumberList->FindObject("prof4NumTrackAfterCut"));

        for (int j = 0; j < nbin; j++) {
            TwoPCorrByRun[i][eV11][j] = prof2->GetBinContent(prof2->GetXaxis()->FindBin("<v1^2>"), j + 1);
            TwoPCorrByRun[i][eV22][j] = prof2->GetBinContent(prof2->GetXaxis()->FindBin("<v2^2>"), j + 1);
            TwoPCorrByRun[i][eV33][j] = prof2->GetBinContent(prof2->GetXaxis()->FindBin("<v3^2>"), j + 1);
            TwoPCorrByRun[i][eV44][j] = prof2->GetBinContent(prof2->GetXaxis()->FindBin("<v4^2>"), j + 1);
            ThreePCorrByRun[i][eV211][j] = prof3->GetBinContent(prof3->GetXaxis()->FindBin("<v2v1v1>"), j + 1);
            ThreePCorrByRun[i][eV321][j] = prof3->GetBinContent(prof3->GetXaxis()->FindBin("<v3v2v1>"), j + 1);
            ThreePCorrByRun[i][eV422][j] = prof3->GetBinContent(prof3->GetXaxis()->FindBin("<v4v2v2>"), j + 1);
            FourPCorrByRun[i][eV2222][j] = prof4->GetBinContent(prof4->GetXaxis()->FindBin("<v2^2v2^2>"), j + 1);
            FourPCorrByRun[i][eV3333][j] = prof4->GetBinContent(prof4->GetXaxis()->FindBin("<v3^2v3^2>"), j + 1);
            FourPCorrByRun[i][eV4444][j] = prof4->GetBinContent(prof4->GetXaxis()->FindBin("<v4^2v4^2>"), j + 1);
            FourPCorrByRun[i][eV3232][j] = prof4->GetBinContent(prof4->GetXaxis()->FindBin("<v2^2v3^2>"), j + 1);
            FourPCorrByRun[i][eV4242][j] = prof4->GetBinContent(prof4->GetXaxis()->FindBin("<v2^2v4^2>"), j + 1);
            nTracks[j] = prof2->GetYaxis()->GetBinCenter(j + 1);
        }

    }

    // bootstrap
    int nValid[n2 + n3 + n4][nbin] = {};

    for (int run = 0; run < nRun; run++) {
        for (int bin = 0; bin < nbin; bin++) {
            for (int obs = 0; obs < n2; obs++) {
                double value = TwoPCorrByRun[run][obs][bin];
                if (value == 0.) {
                    continue;
                }
                TwoPCorr[obs][bin] += value;
                nValid[obs][bin]++;
            }
            for (int obs = 0; obs < n3; obs++) {
                double value = ThreePCorrByRun[run][obs][bin];
                if (value == 0.) {
                    continue;
                }
                ThreePCorr[obs][bin] += value;
                nValid[n2 + obs][bin]++;
            }
            for (int obs = 0; obs < n4; obs++) {
                double value = FourPCorrByRun[run][obs][bin];
                if (value == 0.) {
                    continue;
                }
                FourPCorr[obs][bin] += value;
                nValid[n2 + n3 + obs][bin]++;
            }
        }
    }

    for (int bin = 0; bin < nbin; bin++) {
        for (int obs = 0; obs < n2; obs++) {
            TwoPCorr[obs][bin] = (nValid[obs][bin] > 0) ? TwoPCorr[obs][bin] / nValid[obs][bin] : nan;
        }
        for (int obs = 0; obs < n3; obs++) {
            ThreePCorr[obs][bin] = (nValid[n2 + obs][bin] > 0) ? ThreePCorr[obs][bin] / nValid[n2 + obs][bin] : nan;
        }
        for (int obs = 0; obs < n4; obs++) {
            FourPCorr[obs][bin] = (nValid[n2 + n3 + obs][bin] > 0) ? FourPCorr[obs][bin] / nValid[n2 + n3 + obs][bin] : nan;
        }
    }

    // plot

    TGraph* g2[n2];
    TGraph* g3[n3];
    TGraph* g4[n4];
    TCanvas* c = new TCanvas("c", "c", 800, 600);

    for(int i = 0; i < n2; i++) {
        g2[i] = new TGraph(nbin, nTracks, TwoPCorr[i]);
        g2[i]->Draw();
        g2[i]->GetYaxis()->SetTitle(Form("#LT#LT2#GT#GT_{%d,%d}", TwoPCorrName[i][0], TwoPCorrName[i][1]));
        c->Print(Form("result/temp_2p_%d%d.png", TwoPCorrName[i][0], TwoPCorrName[i][1]));
        c->Clear();
        delete g2[i];
    }
    for(int i = 0; i < n3; i++) {
        g3[i] = new TGraph(nbin, nTracks, ThreePCorr[i]);
        g3[i]->Draw();
        g3[i]->GetYaxis()->SetTitle(Form("#LT#LT3#GT#GT_{%d,%d,%d}", ThreePCorrName[i][0], ThreePCorrName[i][1], ThreePCorrName[i][2]));
        c->Print(Form("result/temp_3p_%d%d%d.png", ThreePCorrName[i][0], ThreePCorrName[i][1], ThreePCorrName[i][2]));
        c->Clear();
        delete g3[i];
    }
    for(int i = 0; i < n4; i++) {
        g4[i] = new TGraph(nbin, nTracks, FourPCorr[i]);
        g4[i]->Draw();
        g4[i]->GetYaxis()->SetTitle(Form("#LT#LT4#GT#GT_{%d,%d,%d,%d}", FourPCorrName[i][0], FourPCorrName[i][1], FourPCorrName[i][2], FourPCorrName[i][3]));
        c->Print(Form("result/temp_4p_%d%d%d%d.png", FourPCorrName[i][0], FourPCorrName[i][1], FourPCorrName[i][2], FourPCorrName[i][3]));
        c->Clear();
        delete g4[i];
    }

    return 1;
}