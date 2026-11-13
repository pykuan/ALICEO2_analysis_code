#include <iostream>
#include <fstream> 
#include <complex>
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

using namespace std;

// try to combine two CalculateData but failed

// MC (LHC25f3_gp): 544886, 544868, 544122, 544116
// real (LHC23zzh_pass5): 544124, 544123, 544122, 544121, 544116, 544098, 544095, 544091
// with phi and pt weight: 544116, 544122
// with phi weight only: 544091, 544095, 544098, 544121, 544123, 544124
// with pt weight only: 544868, 544886

// observables vs. centrality
// realdata
// - vn{k}, NSC
// - vn{2} deta = 0, 0.4, 0.8, 1.0
// MCdata
// - vn{k}, NSC
// - vn{2} deta = 0, 0.4, 0.8, 1.0

// workflow:
// A. generate data (/result/barlow/result_a_realMC.csv, /result/barlow/result_b_realMC.csv, /result/barlow/result_c_realMC.csv) for Barlow test
// B. run Barlow test (BarlowTestC.C)
// C. read Barlow test result (/result/barlow/barlow_aVSc_realMC.csv, /result/barlow/barlow_bVSc_realMC.csv), generate data csv file (/result/data/result_vnscVScent_realMC.csv)

char* realMC = nullptr;
int nRunNumber;
std::vector<int> runNumber;
ofstream resultForBarlowFile;
int nBarlowTests;
std::vector<std::string> dirBarlowTests;
ofstream resultFile;
const int maxRunNumber = 8;

const std::vector<std::string> readBarlowTestTypeReal = {"phiVSall", "phiptVSall"};
const std::vector<std::string> readBarlowTestTypeMC = {"ptVSall", "phiptVSall"};

bool bSCnormalized = true;
const bool bLocal = false;

const int idReal = 0;
const int idMC = 0;
const char* dataDir = "/scratch3/go52dab/O2tutorial/test20260922_analysis/";

const int nCentBins = 8;
const char* centBinsName[nCentBins] = {"00-05", "05-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70"};
const double centBinsCenter[nCentBins] = {2.5, 7.5, 15., 25., 35., 45., 55., 65.};

const int nObservables = 18;
const char* observablesList[nObservables] = {"v2{2}", "v2{4}", "v3{2}", "v3{4}", "v4{2}", "v4{4}", "NSC(2,3)", "NSC(2,4)", "NSC(2,3,4)", "v2{2}gap1", "v2{2}gap2", "v2{2}gap3", "v3{2}gap1", "v3{2}gap2", "v3{2}gap3", "v4{2}gap1", "v4{2}gap2", "v4{2}gap3"};
enum observablesEnum {eV22 = 0, eV24, eV32, eV34, eV42, eV44, eNSC23, eNSC24, eNSC234, eV22gap1, eV22gap2, eV22gap3, eV32gap1, eV32gap2, eV32gap3, eV42gap1, eV42gap2, eV42gap3};
// observables:
//  vn{2} from <<2>>n
//  vn{4} from <<2>>n and <<4>>nn
//  NSC(n,m) from <<2>>n, <<2>>m and <<4>>nm
//  NSC(k,l,m) from <<2>>k, <<2>>l, <<2>>m, <<4>>kl, <<4>>km, <<4>>lm and <<6>>klm
//  vn{2} with 3 different eta gaps

bool bValidRun[maxRunNumber];

std::complex<double> MeanByRun[nObservables][maxRunNumber][nCentBins];
double Mean[nObservables][nCentBins];
double StatErr[nObservables][nCentBins];
double SysErr[nObservables][nCentBins];
double Err[nObservables][nCentBins];

template <size_t nObservables_, size_t nRunNumber_, size_t nCentBins_>
void FillArray(std::complex<double> (&MeanByRun_)[nObservables_][nRunNumber_][nCentBins_],
               bool (&bValidRun_)[nRunNumber_],
               const char* (&centBinsName_)[nCentBins_]){

    const double nan = std::numeric_limits<double>::quiet_NaN();

    for (int obs = 0; obs < nObservables_; ++obs) {
        for (int run = 0; run < nRunNumber_; ++run) {
            for (int cent = 0; cent < nCentBins_; ++cent) {
                MeanByRun_[obs][run][cent] = nan;
            }
        }
    }

    for(int i=0; i<nRunNumber_; i++) {
        bValidRun_[i] = true;
    }

    int id = 0;
    if (strcmp(realMC, "realdata") == 0) {
        id = idReal;
    }
    else if (strcmp(realMC, "MCdata") == 0) {
        id = idMC;
    }

    TFile* analysisFile = nullptr;

    for(int i=0; i<nCentBins_; i++){ //cent bins
        printf("cent = %s \n", centBinsName_[i]);

        if(bLocal){

            if (gSystem->AccessPathName(Form("%s%s/cent%s/AnalysisResults.root", dataDir, realMC, centBinsName_[i]), kFileExists)) {
                printf("%s/cent%s/AnalysisResults.root not exist\n", realMC, centBinsName_[i]);
                continue;
            }

            else {
                printf("found %s/cent%s/AnalysisResults.root\n", realMC, centBinsName_[i]);
                analysisFile = TFile::Open(Form("%s%s/cent%s/AnalysisResults.root", dataDir, realMC, centBinsName_[i]), "read");
            }
        }

        TList* runNumberList = nullptr;
        TProfile* prof2 = nullptr;
        TProfile* prof4 = nullptr;
        TProfile* prof6 = nullptr;
        TProfile* prof2gap1 = nullptr;
        TProfile* prof2gap2 = nullptr;
        TProfile* prof2gap3 = nullptr;
        for(int j=0; j<runNumber.size(); j++) { //run number

            printf(" Run %d: \n", runNumber[j]);

            if (!bLocal){

                if (gSystem->AccessPathName(Form("%s%s/%d/AnalysisResults.root", dataDir, realMC, runNumber[j]), kFileExists)) {
                    continue;

                }
                else {
                    analysisFile = TFile::Open(Form("%s%s/%d/AnalysisResults.root", dataDir, realMC, runNumber[j]), "read");
                    printf(" found %s/%d/AnalysisResults.root\n", realMC, runNumber[j]);
                }
            }

            TString dirName = Form("multiparticle-cumulants_cent%s_id%d", centBinsName_[i], id);
            TDirectory* dir = dynamic_cast<TDirectory*>(analysisFile->Get(dirName.Data()));
            if (!dir) {
                printf(" Failed to open dir %s\n", dirName.Data());
                dirName = Form("multiparticle-cumulants_cent%s", centBinsName_[i]);
                dir = dynamic_cast<TDirectory*>(analysisFile->Get(dirName.Data()));
                if (!dir) {
                    printf(" Failed to open dir %s\n", dirName.Data());
                    dirName = "multiparticle-cumulants";
                    dir = dynamic_cast<TDirectory*>(analysisFile->Get(dirName.Data()));        
                    if (!dir) {
                        printf(" Failed to open dir %s\n", dirName.Data());
                        bValidRun_[j] = false;
                        continue;
                    }
                    else {
                        printf(" Open dir %s\n", dirName.Data());
                    }
                }
                else {
                    printf(" Open dir %s\n", dirName.Data());
                }
            }
            else {
                printf(" Open dir %s\n", dirName.Data());
            }

            TList* baseList = nullptr;
            TList* histogramList = nullptr;
            dir->GetObject("Default list name", baseList);
            if(!baseList) {
                printf(" Failed to find list \"Default list name\"\n");
                bValidRun_[j] = false;
                continue;
            }
            histogramList = dynamic_cast<TList*>(baseList->FindObject("MultiparticleCorrelationProfiles"));
            if(!histogramList) {
                printf(" Failed to find list \"MultiparticleCorrelationProfiles\"\n");
                bValidRun_[j] = false;
                continue;
            }

            runNumberList = dynamic_cast<TList*>(histogramList->FindObject(Form("mcc_run%d", runNumber[j])));

            if(!runNumberList) {
                continue;
            }

            else {

                prof2 = dynamic_cast<TProfile*>(runNumberList->FindObject("prof2AfterCut"));
                prof4 = dynamic_cast<TProfile*>(runNumberList->FindObject("prof4AfterCut"));
                prof6 = dynamic_cast<TProfile*>(runNumberList->FindObject("prof6AfterCut"));
                prof2gap1 = dynamic_cast<TProfile*>(runNumberList->FindObject("prof2AfterCutGapped0.4"));
                prof2gap2 = dynamic_cast<TProfile*>(runNumberList->FindObject("prof2AfterCutGapped0.8"));
                prof2gap3 = dynamic_cast<TProfile*>(runNumberList->FindObject("prof2AfterCutGapped1.0"));

                if (!prof2 || !prof4 || !prof6 || !prof2gap1 || !prof2gap2 || !prof2gap3) {
                    continue;
                }

                // <<2>> and vn{2}
                std::complex<double> v2TwoPCu = std::complex<double>{prof2->GetBinContent(prof2->GetXaxis()->FindBin("<v2^2>")), 0}; // c2{2}
                std::complex<double> v3TwoPCu = std::complex<double>{prof2->GetBinContent(prof2->GetXaxis()->FindBin("<v3^2>")), 0};
                std::complex<double> v4TwoPCu = std::complex<double>{prof2->GetBinContent(prof2->GetXaxis()->FindBin("<v4^2>")), 0};
                
                std::complex<double> v2TwoPCorr = std::pow(v2TwoPCu, 0.5); // v2{2}
                std::complex<double> v3TwoPCorr = std::pow(v3TwoPCu, 0.5); // v3{2}
                std::complex<double> v4TwoPCorr = std::pow(v4TwoPCu, 0.5); // v4{2}

                // <<2>> and vn{2}, gapped
                std::complex<double> v2TwoPCuGap1 = std::complex<double>{prof2gap1->GetBinContent(prof2gap1->GetXaxis()->FindBin("<v2^2>")), 0};
                std::complex<double> v3TwoPCuGap1 = std::complex<double>{prof2gap1->GetBinContent(prof2gap1->GetXaxis()->FindBin("<v3^2>")), 0};
                std::complex<double> v4TwoPCuGap1 = std::complex<double>{prof2gap1->GetBinContent(prof2gap1->GetXaxis()->FindBin("<v4^2>")), 0};
                std::complex<double> v2TwoPCorrGap1 = std::pow(v2TwoPCuGap1, 0.5); // v2{2}
                std::complex<double> v3TwoPCorrGap1 = std::pow(v3TwoPCuGap1, 0.5); // v3{2}
                std::complex<double> v4TwoPCorrGap1 = std::pow(v4TwoPCuGap1, 0.5); // v4{2}

                std::complex<double> v2TwoPCuGap2 = std::complex<double>{prof2gap2->GetBinContent(prof2gap2->GetXaxis()->FindBin("<v2^2>")), 0};
                std::complex<double> v3TwoPCuGap2 = std::complex<double>{prof2gap2->GetBinContent(prof2gap2->GetXaxis()->FindBin("<v3^2>")), 0};
                std::complex<double> v4TwoPCuGap2 = std::complex<double>{prof2gap2->GetBinContent(prof2gap2->GetXaxis()->FindBin("<v4^2>")), 0};
                std::complex<double> v2TwoPCorrGap2 = std::pow(v2TwoPCuGap2, 0.5); // v2{2}
                std::complex<double> v3TwoPCorrGap2 = std::pow(v3TwoPCuGap2, 0.5); // v3{2}
                std::complex<double> v4TwoPCorrGap2 = std::pow(v4TwoPCuGap2, 0.5); // v4{2}

                std::complex<double> v2TwoPCuGap3 = std::complex<double>{prof2gap3->GetBinContent(prof2gap3->GetXaxis()->FindBin("<v2^2>")), 0};
                std::complex<double> v3TwoPCuGap3 = std::complex<double>{prof2gap3->GetBinContent(prof2gap3->GetXaxis()->FindBin("<v3^2>")), 0};
                std::complex<double> v4TwoPCuGap3 = std::complex<double>{prof2gap3->GetBinContent(prof2gap3->GetXaxis()->FindBin("<v4^2>")), 0};
                std::complex<double> v2TwoPCorrGap3 = std::pow(v2TwoPCuGap3, 0.5); // v2{2}
                std::complex<double> v3TwoPCorrGap3 = std::pow(v3TwoPCuGap3, 0.5); // v3{2}
                std::complex<double> v4TwoPCorrGap3 = std::pow(v4TwoPCuGap3, 0.5); // v4{2}

                // <<4>>
                double v2v2FourPCorr = prof4->GetBinContent(prof4->GetXaxis()->FindBin("<v2^2v2^2>")); // <v2^2 v2^2>
                double v2v3FourPCorr = prof4->GetBinContent(prof4->GetXaxis()->FindBin("<v2^2v3^2>")); // <v2^2 v3^2>
                double v2v4FourPCorr = prof4->GetBinContent(prof4->GetXaxis()->FindBin("<v2^2v4^2>")); // <v2^2 v4^2>
                double v3v3FourPCorr = prof4->GetBinContent(prof4->GetXaxis()->FindBin("<v3^2v3^2>")); // <v3^2 v3^2>
                double v3v4FourPCorr = prof4->GetBinContent(prof4->GetXaxis()->FindBin("<v3^2v4^2>")); // <v3^2 v4^2>
                double v4v4FourPCorr = prof4->GetBinContent(prof4->GetXaxis()->FindBin("<v4^2v4^2>")); // <v4^2 v4^2>

                // <<6>>
                double vnSixPCorr = prof6->GetBinContent(prof6->GetXaxis()->FindBin("<v2^2v3^2v4^2>")); // <v2^2 v3^2 v4^2>

                // vn{4}

                std::complex<double> v2FourPCu = 2. * std::pow(v2TwoPCorr, 4.) - v2v2FourPCorr; // -c2{4}
                std::complex<double> v3FourPCu = 2. * std::pow(v3TwoPCorr, 4.) - v3v3FourPCorr; // -c3{4}
                std::complex<double> v4FourPCu = 2. * std::pow(v4TwoPCorr, 4.) - v4v4FourPCorr; // -c4{4}

                std::complex<double> v2FourPCorr = std::pow(v2FourPCu, 0.25); // v2{4}
                std::complex<double> v3FourPCorr = std::pow(v3FourPCu, 0.25); // v3{4}
                std::complex<double> v4FourPCorr = std::pow(v4FourPCu, 0.25); // v4{4}

                // (N)SC
                double NSC23;
                double NSC24;
                double NSC234;

                if (bSCnormalized) {
                    NSC23 = (v2v3FourPCorr - v2TwoPCu.real() * v3TwoPCu.real()) / v2TwoPCuGap3.real() * v3TwoPCuGap3.real(); // NSC(2,3)
                    NSC24 = (v2v4FourPCorr - v2TwoPCu.real() * v4TwoPCu.real()) / v2TwoPCuGap3.real() * v4TwoPCuGap3.real(); // NSC(2,4)
                    NSC234 = (vnSixPCorr - v2v3FourPCorr * v4TwoPCu.real() - v2v4FourPCorr * v3TwoPCu.real() - v3v4FourPCorr * v2TwoPCu.real() + 2 * v2TwoPCu.real() * v3TwoPCu.real() * v4TwoPCu.real()) / (v2TwoPCuGap3.real() * v3TwoPCuGap3.real() * v4TwoPCuGap3.real()); // NSC(2,3,4)
                }
                else {
                    NSC23 = (v2v3FourPCorr - v2TwoPCu.real() * v3TwoPCu.real()); // SC(2,3)
                    NSC24 = (v2v4FourPCorr - v2TwoPCu.real() * v4TwoPCu.real()); // SC(2,4)
                    NSC234 = (vnSixPCorr - v2v3FourPCorr * v4TwoPCu.real() - v2v4FourPCorr * v3TwoPCu.real() - v3v4FourPCorr * v2TwoPCu.real() + 2 * v2TwoPCu.real() * v3TwoPCu.real() * v4TwoPCu.real()); // SC(2,3,4)
                }

                // printf("  v2{2} = %.12f \n", v2TwoPCorr);
                // printf("  v3{2} = %.12f \n", v3TwoPCorr);
                // printf("  v4{2} = %.12f \n", v4TwoPCorr);
                // printf("  v2{4} = %.12f \n", v2FourPCorr);
                // printf("  v3{4} = %.12f \n", v3FourPCorr);
                // printf("  v4{4} = %.12f \n", v4FourPCorr);
                // // printf("  <v2^2 v2^2> = %.12f \n", v2v2FourPCorr);
                // // printf("  <v2^2 v3^2> = %.12f \n", v2v3FourPCorr);
                // // printf("  <v2^2 v4^2> = %.12f \n", v2v4FourPCorr);
                // printf("  (N)SC(2,3) = %.12f \n", NSC23);
                // printf("  (N)SC(2,4) = %.12f \n", NSC24);
                // printf("  (N)SC(2,3,4) = %.12f \n", NSC234);

                MeanByRun_[eV22][j][i] = v2TwoPCorr;
                MeanByRun_[eV24][j][i] = v2FourPCorr;
                MeanByRun_[eV32][j][i] = v3TwoPCorr;
                MeanByRun_[eV34][j][i] = v3FourPCorr;
                MeanByRun_[eV42][j][i] = v4TwoPCorr;
                MeanByRun_[eV44][j][i] = v4FourPCorr;
                MeanByRun_[eNSC23][j][i] = NSC23;
                MeanByRun_[eNSC24][j][i] = NSC24;
                MeanByRun_[eNSC234][j][i] = NSC234;
                MeanByRun_[eV22gap1][j][i] = v2TwoPCorrGap1;
                MeanByRun_[eV32gap1][j][i] = v3TwoPCorrGap1;
                MeanByRun_[eV42gap1][j][i] = v4TwoPCorrGap1;
                MeanByRun_[eV22gap2][j][i] = v2TwoPCorrGap2;
                MeanByRun_[eV32gap2][j][i] = v3TwoPCorrGap2;
                MeanByRun_[eV42gap2][j][i] = v4TwoPCorrGap2;
                MeanByRun_[eV22gap3][j][i] = v2TwoPCorrGap3;
                MeanByRun_[eV32gap3][j][i] = v3TwoPCorrGap3;
                MeanByRun_[eV42gap3][j][i] = v4TwoPCorrGap3;
            }

            if (!bLocal && analysisFile) {
                analysisFile->Close();
                delete analysisFile;
                analysisFile = nullptr;
            }

        }

        if (bLocal && analysisFile) {
            analysisFile->Close();
            delete analysisFile;
            analysisFile = nullptr;
        }

    }

}

template <size_t nObservables_, size_t nRunNumber_, size_t nCentBins_>
void Bootstrap(std::complex<double> (&MeanByRun_)[nObservables_][nRunNumber_][nCentBins_],
               const char* (&centBinsName_)[nCentBins_],
               double (&Mean_)[nObservables_][nCentBins_],
               double (&StatErr_)[nObservables_][nCentBins_],
               double (&Err_)[nObservables_][nCentBins_],
               const char* (&observablesList_)[nObservables_])
{

    double nValid[nObservables_][nCentBins] = {}; // actually int

    for (int obs = 0; obs < nObservables_; ++obs) {
        for (int cent = 0; cent < nCentBins; ++cent) {
            Mean_[obs][cent] = 0.;
            StatErr_[obs][cent] = 0.;
        }
    }

    for (int run = 0; run < nRunNumber_; ++run) {
        for (int obs = 0; obs < nObservables_; ++obs) {
            for (int cent = 0; cent < nCentBins; ++cent) {

                std::complex<double> value = MeanByRun_[obs][run][cent];

                if (value != value) { // nan?
                    continue;
                }

                Mean_[obs][cent] += value.real();
                nValid[obs][cent]++;
            }
        }
    }

    for (int obs = 0; obs < nObservables_; ++obs) {
        for (int cent = 0; cent < nCentBins; ++cent) {

            if (nValid[obs][cent] > 0) {
                Mean_[obs][cent] /= nValid[obs][cent];
            } else {
                Mean_[obs][cent] = std::numeric_limits<double>::quiet_NaN();
            }
        }
    }

    for (int run = 0; run < nRunNumber_; ++run) {
        for (int obs = 0; obs < nObservables_; ++obs) {
            for (int cent = 0; cent < nCentBins; ++cent) {

                std::complex<double> value = MeanByRun_[obs][run][cent];

                if (value != value || nValid[obs][cent] < 2) { // !std::isfinite(value) || ...
                    continue;
                }

                StatErr_[obs][cent] += std::pow(value.real() - Mean_[obs][cent], 2);
            }
        }
    }

    for (int obs = 0; obs < nObservables_; ++obs) {
        for (int cent = 0; cent < nCentBins; ++cent) {

            const int n = nValid[obs][cent];

            if (n >= 2) {
                StatErr_[obs][cent] = std::sqrt(StatErr_[obs][cent] /(static_cast<double>(n) * (n - 1)));
            } else {
                StatErr_[obs][cent] = std::numeric_limits<double>::quiet_NaN(); //0
                if (!bLocal) {
                    Mean_[obs][cent] = std::numeric_limits<double>::quiet_NaN();
                }
            }

            Err_[obs][cent] = StatErr_[obs][cent];

            printf("%s, cent %s: %d valid runs\n", observablesList_[obs], centBinsName_[cent], n);
        }
    }
}



template <size_t nObservables_, size_t nCentBins_>
void Barlow(
               double (&Mean_)[nObservables_][nCentBins_],
               double (&StatErr_)[nObservables_][nCentBins_],
               double (&SysErr_)[nObservables_][nCentBins_],
               double (&Err_)[nObservables_][nCentBins_],
               const char* (&observablesList_)[nObservables_],
               const std::vector<std::string>& dirBarlowTests_
){

    for (int i=0; i<nCentBins_; i++) {
        for(int j=0; j<nObservables_; j++) {
            SysErr_[j][i] = 0.;
        }
    }

    double a[nObservables_][nCentBins_] = {};
    for (int m=0; m<nObservables_; m++) {
        for (const auto& path : dirBarlowTests_) {
            auto array = ReadBarlowTest(path.c_str(), observablesList_[m]);
            for (int p=0; p<nCentBins_; p++) {
                a[m][p] += std::pow(array[p], 2);
            }
        }

        for (int q=0; q<nCentBins_; q++) {
            SysErr_[m][q] = a[m][q];
        }

    }

    for (int p=0; p<nCentBins_; p++) {
        for (int m=0; m<nObservables_; m++) {
            SysErr_[m][p] = std::pow(SysErr_[m][p], 0.5);
            SysErr_[m][p] *= Mean_[m][p];
            Err_[m][p] = StatErr_[m][p] + SysErr_[m][p];
        }
    }

    for (int i=0; i<nCentBins_; i++) {
        for(int j=0; j<nObservables_; j++) {
            if (SysErr_[j][i] != SysErr_[j][i]) { // check for NaN
                SysErr_[j][i] = 0.;
            }
        }
    }

}

template <size_t nObservables_, size_t nRunNumber_, size_t nCentBins_>
void Print(std::complex<double> (&MeanByRun_)[nObservables_][nRunNumber_][nCentBins_],
           bool (&bValidRun_)[nRunNumber_],
           const char* (&centBinsName_)[nCentBins_],
           double (&Mean_)[nObservables_][nCentBins_],
           double (&StatErr_)[nObservables_][nCentBins_],
           double (&Err_)[nObservables_][nCentBins_],
           const char* (&observablesList_)[nObservables_],
           const double (&centBinsCenter_)[nCentBins_],
           double (&SysErr_)[nObservables_][nCentBins_],
           const std::vector<std::string>& dirBarlowTests_){

    if (!bLocal) {
        resultForBarlowFile << "#X,Y,staterror" << endl;
    }

    FillArray(MeanByRun_, bValidRun_, centBinsName_);
    Bootstrap(MeanByRun_, centBinsName_, Mean_, StatErr_, Err_, observablesList_);
    if (!bLocal) {
        Barlow(Mean_, StatErr_, SysErr_, Err_, observablesList_, dirBarlowTests_);
    }

    for(int m=0; m<nObservables_; m++) {
        if (!bLocal) {
            resultForBarlowFile << "#:" << observablesList_[m] << endl;
        }
        resultFile << "#:" << observablesList_[m] << endl;
        resultFile << "X,DUMMY,DUMMY," << observablesList_[m] << ",stat +,stat -,sys +,sys -" << endl;
        for (int i=0; i<nCentBins_; i++) {
            if (!bLocal) {
                resultForBarlowFile << centBinsCenter_[i] << "," << Mean_[m][i] << "," << StatErr_[m][i] << endl;
            }
            resultFile << centBinsCenter_[i] << "," // CENTRALITY [PCT]
                        << 0. << "," // CENTRALITY [PCT] LOW
                        << 0. << "," // CENTRALITY [PCT] HIGH
                        << Mean_[m][i] << "," // V2{2}
                        << StatErr_[m][i] << "," // stat +
                        << - StatErr_[m][i] << "," // stat -
                        << SysErr_[m][i] << "," // sys +
                        << - SysErr_[m][i] // sys -
                        << endl;
        }
        if (!bLocal) {
            resultForBarlowFile << "#end" << endl;
        }
        resultFile << endl;
    }

    if (!bLocal) {
        resultForBarlowFile.close();
    }
    resultFile.close();

}

template <size_t nObservables_, size_t nRunNumber_temp, size_t nCentBins_>
void CalculateData(const char* realMC_,
                   const char* generateBarlowTestType,
                   
                   const char* (&centBinsName_)[nCentBins_],
                   const double (&centBinsCenter_)[nCentBins_],
                   bool (&bValidRun_)[nRunNumber_temp],
                   std::complex<double> (&MeanByRun_)[nObservables_][nRunNumber_temp][nCentBins_],
                   
                   double (&Mean_)[nObservables_][nCentBins_],
                   double (&StatErr_)[nObservables_][nCentBins_],
                   double (&SysErr_)[nObservables_][nCentBins_],
                   double (&Err_)[nObservables_][nCentBins_],
                   const char* (&observablesList_)[nObservables_],
                   
                   const std::vector<std::string>& dirBarlowTests_
                
                ) {


    std::vector<std::string> readBarlowTestType;

    if (strcmp(realMC_, "realdata") == 0) {
        realMC = "realdata";
        readBarlowTestType = readBarlowTestTypeReal;
        if (strcmp(generateBarlowTestType, "all") == 0) {
            runNumber = {544124, 544123, 544122, 544121, 544116, 544098, 544095, 544091};
        }
        else if (strcmp(generateBarlowTestType, "phipt") == 0) {
            runNumber = {544116, 544122};
        }
        else if (strcmp(generateBarlowTestType, "phi") == 0) {
            runNumber = {544091, 544095, 544098, 544121, 544123, 544124};
        }
    }
    else if (strcmp(realMC_, "MCdata") == 0) {
        realMC = "MCdata";
        readBarlowTestType = readBarlowTestTypeMC;
        if (strcmp(generateBarlowTestType, "all") == 0) {
            runNumber = {544886, 544868, 544122, 544116};
        }
        else if (strcmp(generateBarlowTestType, "phipt") == 0) {
            runNumber = {544116, 544122};
        }
        else if (strcmp(generateBarlowTestType, "pt") == 0) {
            runNumber = {544886, 544868};
        }
    }
    nRunNumber = runNumber.size();
    resultForBarlowFile.open(Form("result/barlow/result_%s_%s.csv", generateBarlowTestType, realMC_));
    nBarlowTests = readBarlowTestType.size();
    for (const auto& barlowTest : readBarlowTestType){
        dirBarlowTests.emplace_back(Form("result/barlow/barlow_%s_%s.csv", barlowTest.c_str(), realMC_));
    }
    resultFile.open(Form("result/data/result_vnscVScent_%s.csv", realMC_));

    Print(MeanByRun_, bValidRun_, centBinsName_, Mean_, StatErr_, Err_, observablesList_, centBinsCenter_, SysErr_, dirBarlowTests_);
    dirBarlowTests.clear();
}

int CalculateData_vnscVScent(bool isNormalized = true) {

    bSCnormalized = isNormalized;

    if (!bLocal) {
        // generate data for Barlow test

        CalculateData("realdata", "all", centBinsName, centBinsCenter, bValidRun, MeanByRun, Mean, StatErr, SysErr, Err, observablesList, dirBarlowTests);
        CalculateData("realdata", "phipt", centBinsName, centBinsCenter, bValidRun, MeanByRun, Mean, StatErr, SysErr, Err, observablesList, dirBarlowTests);
        CalculateData("realdata", "phi", centBinsName, centBinsCenter, bValidRun, MeanByRun, Mean, StatErr, SysErr, Err, observablesList, dirBarlowTests);
        CalculateData("MCdata", "all", centBinsName, centBinsCenter, bValidRun, MeanByRun, Mean, StatErr, SysErr, Err, observablesList, dirBarlowTests);
        CalculateData("MCdata", "phipt", centBinsName, centBinsCenter, bValidRun, MeanByRun, Mean, StatErr, SysErr, Err, observablesList, dirBarlowTests);
        CalculateData("MCdata", "pt", centBinsName, centBinsCenter, bValidRun, MeanByRun, Mean, StatErr, SysErr, Err, observablesList, dirBarlowTests);

        // generate Barlow test result
        for (const auto& barlowTest : readBarlowTestTypeReal) {
            std::size_t pos = barlowTest.find("VS");
            std::string strTrial = barlowTest.substr(0, pos);
            std::string strDefault = barlowTest.substr(pos + 2);
            // trial, default, sys err data path
            BarlowTest(Form("result/barlow/result_%s_realdata.csv", strTrial.c_str()), Form("result/barlow/result_%s_realdata.csv", strDefault.c_str()), Form("result/barlow/barlow_%s_realdata.csv", barlowTest.c_str()));
        }
        for (const auto& barlowTest : readBarlowTestTypeMC) {
            std::size_t pos = barlowTest.find("VS");
            std::string strTrial = barlowTest.substr(0, pos);
            std::string strDefault = barlowTest.substr(pos + 2);
            BarlowTest(Form("result/barlow/result_%s_MCdata.csv", strTrial.c_str()), Form("result/barlow/result_%s_MCdata.csv", strDefault.c_str()), Form("result/barlow/barlow_%s_MCdata.csv", barlowTest.c_str()));
        }
    }

    // apply systematic error and print total result
    CalculateData("realdata", "all", centBinsName, centBinsCenter, bValidRun, MeanByRun, Mean, StatErr, SysErr, Err, observablesList, dirBarlowTests);
    if (!bLocal) {
        CalculateData("MCdata", "all", centBinsName, centBinsCenter, bValidRun, MeanByRun, Mean, StatErr, SysErr, Err, observablesList, dirBarlowTests);
    }

    return 0;
}
