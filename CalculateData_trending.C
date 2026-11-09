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

using namespace std;

// generate data csv file

// observables vs. runnumber
// realdata
// - multiplicity
// - nParticles
// - pt
// - <<2>> for one or two cent bin

const char* dataDir = "/scratch3/go52dab/O2tutorial/test20260827_analysis/";

const int nCentBins = 8;
const char* centBinsName[nCentBins] = {"00-05", "05-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70"};
const double centBinsCenter[nCentBins] = {2.5, 7.5, 15., 25., 35., 45., 55., 65.};

const int nObservables = 9;
const char* observablesList[nObservables] = {"multiplicity", "Nparticles", "pt", "<<2>>2cent10-20", "<<2>>3cent10-20", "<<2>>4cent10-20", "<<2>>2cent50-60", "<<2>>3cent50-60", "<<2>>4cent50-60"};
enum observablesEnum {eMult = 0, eNParticles, ePt, eTwoP2Cent15, eTwoP3Cent15, eTwoP4Cent15, eTwoP2Cent55, eTwoP3Cent55, eTwoP4Cent55};
// observables:
//  <reference multiplicity>
//  <number of selected particles>
//  <pT>
//  <<2>> for cent 10-20
//  <<2>> for cent 50-60

const int nRunNumber = 8;
const int runNumber[nRunNumber] = {544124, 544123, 544122, 544121, 544116, 544098, 544095, 544091};
// const int nRunNumber = 1;
// const int runNumber[nRunNumber] = {544123};
// MC (LHC25f3_gp): 544886, 544868, 544122, 544116
// real (LHC23zzh_pass5): 544124, 544123, 544122, 544121, 544116, 544098, 544095, 544091
// with phi and pt weight: 544116, 544122
// with phi weight only: 544091, 544095, 544098, 544121, 544123, 544124
// with pt weight only: 544868, 544886


ofstream resultFile("result/data/result_trending_realdata.csv");

double Mean[nObservables][nRunNumber];
double Err[nObservables][nRunNumber];

double dummy = 9999.;

void FillArray(){


    TFile* analysisFile = nullptr;

    for(int j=0; j<nRunNumber; j++) { //run number

        TH1F* multHist = nullptr;
        TH1F* ptHist = nullptr;
        TProfile* prof2 = nullptr;

        // run number
        printf(" Run %d: \n", runNumber[j]);
        if (gSystem->AccessPathName(Form("%srealdata/%d/AnalysisResults.root", dataDir, runNumber[j]), kFileExists)) {
            printf("realdata/%d/AnalysisResults.root not exist\n", runNumber[j]);

        }
        else {
            analysisFile = TFile::Open(Form("%srealdata/%d/AnalysisResults.root", dataDir, runNumber[j]), "read");
            printf("found realdata/%d/AnalysisResults.root\n", runNumber[j]);

            // "multiparticle-cumulants"
            TString dirName = "multiparticle-cumulants";
            TDirectory* dir = dynamic_cast<TDirectory*>(analysisFile->Get(dirName.Data()));
            if (!dir) {
                printf(" Failed to open dir %s\n", dirName.Data());
            }
            else {
                printf(" Open dir %s\n", dirName.Data());

                // "Default list name"
                TList* baseList = nullptr;
                dir->GetObject("Default list name", baseList);
                if(!baseList) {
                    printf("  Failed to find list \"Default list name\"\n");
                }
                else {
                    printf("  Open list \"Default list name\"\n");

                    // "EventHistograms"
                    TList* EVhistogramList = dynamic_cast<TList*>(baseList->FindObject("EventHistograms"));
                    if(!EVhistogramList) {
                        printf("  Failed to find list \"EventHistograms\"\n");
                    }
                    else {
                        printf("  Open list \"EventHistograms\"\n");

                        // "fHistMultiplicity[eRec][after cut]"
                        multHist = dynamic_cast<TH1F*>(EVhistogramList->FindObject("fHistMultiplicity[eRec][after cut]"));
                        if(!multHist) {
                            printf("  Failed to find histogram \"fHistMultiplicity[eRec][after cut]\"\n");
                        }
                        else {
                            printf("  Found histogram \"fHistMultiplicity[eRec][after cut]\"\n");
                            double mult_mean = multHist->GetMean();
                            double mult_err = multHist->GetStdDev();
                            Mean[eMult][j] = mult_mean;
                            Err[eMult][j] = mult_err;

                        }
                    }

                    // "ParticleHistograms"
                    TList* PThistogramList = dynamic_cast<TList*>(baseList->FindObject("ParticleHistograms"));
                    if(!PThistogramList) {
                        printf("  Failed to find list \"ParticleHistograms\"\n");
                    }
                    else {
                        printf("  Open list \"ParticleHistograms\"\n");

                        // "fHistPt[eRec][after cut]"
                        ptHist = dynamic_cast<TH1F*>(PThistogramList->FindObject("fHistPt[eRec][after cut]"));
                        if(!ptHist) {
                            printf("  Failed to find histogram \"fHistPt[eRec][after cut]\"\n");
                        }
                        else {
                            printf("  Found histogram \"fHistPt[eRec][after cut]\"\n");
                            double nparticle = static_cast<double>(ptHist->GetEntries());
                            double pt_mean = ptHist->GetMean();
                            double pt_err = ptHist->GetStdDev();
                            Mean[eNParticles][j] = nparticle;
                            Err[eNParticles][j] = 0.;
                            Mean[ePt][j] = pt_mean;
                            Err[ePt][j] = pt_err;

                        }
                    }

                }

            }

            // "multiparticle-cumulants_cent10-20"
            TString dirName15 = "multiparticle-cumulants_cent10-20";
            TDirectory* dir15 = dynamic_cast<TDirectory*>(analysisFile->Get(dirName15.Data()));
            if (!dir15) {
                printf(" Failed to open dir %s\n", dirName15.Data());
            }
            else {
                printf(" Open dir %s\n", dirName15.Data());

                // "Default list name"
                TList* baseList = nullptr;
                dir15->GetObject("Default list name", baseList);
                if(!baseList) {
                    printf("  Failed to find list \"Default list name\"\n");
                }
                else {
                    printf("  Open list \"Default list name\"\n");

                    // "MultiparticleCorrelationProfiles"
                    TList* MChistogramList = dynamic_cast<TList*>(baseList->FindObject("MultiparticleCorrelationProfiles"));
                    if(!MChistogramList) {
                        printf("  Failed to find list \"MultiparticleCorrelationProfiles\"\n");
                    } 
                    else {
                        printf("  Found list \"MultiparticleCorrelationProfiles\"\n");

                        // "MultiparticleCorrelationProfiles" -> "mcc_run%d"
                        TList* runNumberList = dynamic_cast<TList*>(MChistogramList->FindObject(Form("mcc_run%d", runNumber[j])));
                        if(!runNumberList) {
                            printf("   Failed to find list \"mcc_run%d\"\n", runNumber[j]);
                        }
                        else {
                            printf("   Open list \"mcc_run%d\"\n", runNumber[j]);

                            prof2 = dynamic_cast<TProfile*>(runNumberList->FindObject("prof2AfterCut"));
                            if (!prof2) {
                                printf("    Failed to find prof2");
                            }
                            else {
                                printf("    Found prof2");

                                double v2TwoPCorr = prof2->GetBinContent(prof2->GetXaxis()->FindBin("<v2^2>"));
                                double v3TwoPCorr = prof2->GetBinContent(prof2->GetXaxis()->FindBin("<v3^2>"));
                                double v4TwoPCorr = prof2->GetBinContent(prof2->GetXaxis()->FindBin("<v4^2>"));
                                double v2TwoPCorr_err = prof2->GetStdDev(prof2->GetXaxis()->FindBin("<v2^2>")); //getmeanerror
                                double v3TwoPCorr_err = prof2->GetStdDev(prof2->GetXaxis()->FindBin("<v3^2>"));
                                double v4TwoPCorr_err = prof2->GetStdDev(prof2->GetXaxis()->FindBin("<v4^2>"));

                                Mean[eTwoP2Cent15][j] = v2TwoPCorr;
                                Mean[eTwoP3Cent15][j] = v3TwoPCorr;
                                Mean[eTwoP4Cent15][j] = v4TwoPCorr;
                                Err[eTwoP2Cent15][j] = v2TwoPCorr_err;
                                Err[eTwoP3Cent15][j] = v3TwoPCorr_err;
                                Err[eTwoP4Cent15][j] = v4TwoPCorr_err;

                            }


                        }



                    }

                }

            }


            // "multiparticle-cumulants_cent50-60"
            TString dirName55 = "multiparticle-cumulants_cent50-60";
            TDirectory* dir55 = dynamic_cast<TDirectory*>(analysisFile->Get(dirName55.Data()));
            if (!dir55) {
                printf(" Failed to open dir %s\n", dirName55.Data());
            }
            else {
                printf(" Open dir %s\n", dirName55.Data());

                // "Default list name"
                TList* baseList = nullptr;
                dir55->GetObject("Default list name", baseList);
                if(!baseList) {
                    printf("  Failed to find list \"Default list name\"\n");
                }
                else {
                    printf("  Open list \"Default list name\"\n");

                    // "MultiparticleCorrelationProfiles"
                    TList* MChistogramList = dynamic_cast<TList*>(baseList->FindObject("MultiparticleCorrelationProfiles"));
                    if(!MChistogramList) {
                        printf("  Failed to find list \"MultiparticleCorrelationProfiles\"\n");
                    } 
                    else {
                        printf("  Found list \"MultiparticleCorrelationProfiles\"\n");

                        // "MultiparticleCorrelationProfiles" -> "mcc_run%d"
                        TList* runNumberList = dynamic_cast<TList*>(MChistogramList->FindObject(Form("mcc_run%d", runNumber[j])));
                        if(!runNumberList) {
                            printf("   Failed to find list \"mcc_run%d\"\n", runNumber[j]);
                        }
                        else {
                            printf("   Open list \"mcc_run%d\"\n", runNumber[j]);

                            prof2 = dynamic_cast<TProfile*>(runNumberList->FindObject("prof2AfterCut"));
                            if (!prof2) {
                                printf("    Failed to find prof2");
                            }
                            else {
                                printf("    Found prof2");

                                double v2TwoPCorr = prof2->GetBinContent(prof2->GetXaxis()->FindBin("<v2^2>"));
                                double v3TwoPCorr = prof2->GetBinContent(prof2->GetXaxis()->FindBin("<v3^2>"));
                                double v4TwoPCorr = prof2->GetBinContent(prof2->GetXaxis()->FindBin("<v4^2>"));
                                double v2TwoPCorr_err = prof2->GetStdDev(prof2->GetXaxis()->FindBin("<v2^2>"));
                                double v3TwoPCorr_err = prof2->GetStdDev(prof2->GetXaxis()->FindBin("<v3^2>"));
                                double v4TwoPCorr_err = prof2->GetStdDev(prof2->GetXaxis()->FindBin("<v4^2>"));

                                Mean[eTwoP2Cent55][j] = v2TwoPCorr;
                                Mean[eTwoP3Cent55][j] = v3TwoPCorr;
                                Mean[eTwoP4Cent55][j] = v4TwoPCorr;
                                Err[eTwoP2Cent55][j] = v2TwoPCorr_err;
                                Err[eTwoP3Cent55][j] = v3TwoPCorr_err;
                                Err[eTwoP4Cent55][j] = v4TwoPCorr_err;

                            }

                        }

                    }

                }

            }

        }

    }

}

int CalculateData_trending(){

    FillArray();

    for(int m=0; m<nObservables; m++) {
        resultFile << "#:" << observablesList[m] << endl;
        resultFile << "RUN NUMBER,DUMMY 1,DUMMY 2," << observablesList[m] << ",stat +,stat -,sys +,sys -" << endl;
        for (int i=0; i<nRunNumber; i++) {
            resultFile << runNumber[i] << "," // RUN NUMBER
                        << 0. << "," // DUMMY 1
                        << 0. << "," // DUMMY 2
                        << Mean[m][i] << "," // observable
                        << Err[m][i] << "," // stat +
                        << - Err[m][i] << "," // stat -
                        << 0 << "," // sys +
                        << 0 // sys -
                        << endl;
        }
        resultFile << endl;
    }

    resultFile.close();
    return 0;
}