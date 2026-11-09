// Get phi histograms from AnalysisResults.root of real data LHC23zzh_pass5
// Get pt histograms from AnalysisResults.root of MC data LHC25f3_gp
// Invert the histograms
// Write to weights.root
// Upload weight.root to CCDB
// Use weight.root in CCDB as input weight histograms in analysis

// MC rec data pt histograms are valid if processMonteCarlo
// MC sim data pt histograms are valid if processMonteCarlo -> MC rec pt vs. MC sim pt
// Phi histograms are valid if processRealData

#include <TList.h>
#include <TH1F.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TString.h>
#include <TObject.h>
#include <TIterator.h>

#include <map>
#include <cmath>
#include <cstdio>

using namespace std;
TList* getHistogramList(TFile* file_weight, int id) {

    TString dirName = Form("multiparticle-cumulants_id%d", id);
    TDirectory* dir = dynamic_cast<TDirectory*>(file_weight->Get(dirName.Data()));
    if (!dir) {
        // printf(" Failed to open dir %s\n", dirName.Data());
        dir = dynamic_cast<TDirectory*>(file_weight->Get("multiparticle-cumulants"));
        if (!dir) {
            // printf(" Failed to open dir \"multiparticle-cumulants\"\n");
        }
        else {
            printf(" Found \"multiparticle-cumulants\"\n");
        }
    }
    else {
        printf(" Found %s\n", dirName.Data());
    }

    TList* baseList = nullptr;
    TList* histogramList = nullptr;
    dir->GetObject("Default list name", baseList);
    if(!baseList) {
        printf(" failed to find list \"Default list name\"\n");
    }
    histogramList = dynamic_cast<TList*>(baseList->FindObject("WeightHistograms"));
    if(!histogramList) {
        printf(" failed to find list \"WeightHistograms\"\n");
    }

    return histogramList;

}

int makeWeightHistogram(){

    TFile* file = new TFile("weights.root", "recreate");
    TList* ccdb_object = new TList;
    ccdb_object->SetOwner(kTRUE);
    ccdb_object->SetName("ccdb_object");

    int runNumberMC[4] = {544886, 544868, 544122, 544116}; //544886, 544868, 544122, 544116
    int runNumberReal[8] = {544124, 544123, 544122, 544121, 544116, 544098, 544095, 544091};//544124, 544123, 544122, 544121, 544116, 544098, 544095, 544091
    int idMC = 58663;
    int idReal = 58664;

    std::map<TString, TList*> HistExistMap;
    TFile* file_real_weight = nullptr;
    TFile* file_MC_weight = nullptr;

    //phi weight

    for(int runNumber : runNumberReal) {//test20260706_make_weight_hist

        printf("processing real data:\n");

        file_real_weight = TFile::Open(Form("/scratch3/go52dab/O2tutorial/test20261030_make_weight_hist/realdata/%d/AnalysisResults.root", runNumber), "READ");

        if (!file_real_weight || file_real_weight->IsZombie()) {
            printf(" failed to open real %d/AnalysisResults.root\n", runNumber);
            return 1;
        }

        TList* histogramRealList = getHistogramList(file_real_weight, idReal);
        TList* runListReal = nullptr;
        TIter nextReal(histogramRealList);
        TObject* objectReal = nullptr;
        while(true) {

            objectReal = nextReal();
            if(!objectReal) {
                break;
            }

            TH1F* histReal = nullptr;
            histReal = dynamic_cast<TH1F*>(objectReal);
            if (!histReal) {
                printf("  failed to open real data histogram\n");
                continue;
            }

            TString histNameReal = histReal->GetName();
            printf("  processing %s \n", histNameReal.Data());

            if(histNameReal.BeginsWith("hPhi_run")) { 
                printf("  is phi real histogram\n");  

                TString runNumberName = histNameReal;
                runNumberName.ReplaceAll("hPhi_run", "");

                TH1F* histClone = dynamic_cast<TH1F*>(histReal->Clone());
                TString histCloneName = histNameReal;
                histCloneName.ReplaceAll("hPhi_run", "wPhi_run");
                histClone->SetName(histCloneName);
                histClone->SetTitle(Form("phi weight for run %s", runNumberName.Data()));
                histClone->SetDirectory(nullptr);

                int nBins = histClone->GetNbinsX();
                for(int i=1; i<=nBins; i++){
                    if(histClone->GetBinContent(i)>0){
                        histClone->SetBinContent(i, 1/histClone->GetBinContent(i) * pow(10, 6));
                    }
                    else{
                        histClone->SetBinContent(i, 0);
                    }

                }
                histClone->SetBinContent(0, 0.);
                histClone->SetBinContent(nBins+1, 0.);

                if (HistExistMap.find(runNumberName) == HistExistMap.end()) {
                    runListReal = nullptr;
                    runListReal = new TList;
                    runListReal->SetOwner(kTRUE);
                    runListReal->SetName(Form("%09d", runNumberName.Atoi()));
                    HistExistMap[runNumberName] = runListReal;
                }
                HistExistMap[runNumberName]->Add(histClone);

            }
            else continue;

        }
    }

    // pt weight

    for(int runNumber : runNumberMC) {

        printf("processing MC data:\n");

        file_MC_weight = TFile::Open(Form("/scratch3/go52dab/O2tutorial/test20261030_make_weight_hist/MCdata/%d/AnalysisResults.root", runNumber), "READ");
        
        if (!file_MC_weight || file_MC_weight->IsZombie()) {
            printf(" failed to open MC %d/AnalysisResults.root\n", runNumber);
            return 1;
        }

        TList* histogramMCList = getHistogramList(file_MC_weight, idMC);
        TList* runListMC = nullptr;
        TIter nextMC(histogramMCList);
        TObject* objectMC = nullptr;
        while(true) {

            objectMC = nextMC();
            if(!objectMC) {
                break;
            }

            TH1F* histMC = nullptr;
            histMC = dynamic_cast<TH1F*>(objectMC);
            if (!histMC) {
                printf("  failed to open MC data histogram\n");
                continue;
            }

            TString histNameMC = histMC->GetName();
            printf("  processing %s \n", histNameMC.Data());

            if(histNameMC.BeginsWith("hPtReal_run")) { //MC rec
                printf("  is pt MC rec histogram\n");  
            
                TString runNumberNameMC = histNameMC;
                runNumberNameMC.ReplaceAll("hPtReal_run", ""); //

                TH1F* histCloneMC = dynamic_cast<TH1F*>(histMC->Clone());
                TString histCloneMCName = histNameMC;
                histCloneMCName.ReplaceAll("hPtReal_run", "wPtReal_run");
                histCloneMC->SetName(histCloneMCName);
                histCloneMC->SetTitle(Form("pt weight for run %s", runNumberNameMC.Data()));
                histCloneMC->SetDirectory(nullptr);

                TIter nextMCsim(histogramMCList);
                TObject* objectMCsim = nullptr;
                while(true) {

                    objectMCsim = nextMCsim();
                    if(!objectMCsim) {
                        break;
                    }

                    TH1F* histMCsim = nullptr;
                    histMCsim = dynamic_cast<TH1F*>(objectMCsim);
                    if (!histMCsim) {
                        printf("  failed to open MC data histogram\n");
                        continue;
                    }

                    TString histNameMCsim = histMCsim->GetName();

                    if(histNameMCsim.BeginsWith("hPtMC_run")) {   //MC sim

                        // MC/Real
                        TString runNumberNameMCsim = histNameMCsim;
                        runNumberNameMCsim.ReplaceAll("hPtMC_run", "");

                        if(runNumberNameMCsim == runNumberNameMC) {
                            printf("  paired pt MC sim histogram found\n");  
                            int nBins = histCloneMC->GetNbinsX();
                            for(int i=1; i<=nBins; i++){
                                if(histCloneMC->GetBinContent(i)>0){
                                    double content = histMCsim->GetBinContent(i) / histCloneMC->GetBinContent(i);
                                    double scale = (histCloneMC->GetEntries() / histCloneMC->GetNbinsX()) / (histMCsim->GetEntries() / histMCsim->GetNbinsX()); 
                                    histCloneMC->SetBinContent(i, content); //*
                                }
                                else{
                                    histCloneMC->SetBinContent(i, 0);
                                }
                                // histCloneMC->SetBinError(i, 0.);
                            }
                            histCloneMC->SetBinContent(0, 0.);
                            histCloneMC->SetBinContent(nBins+1, 0.);
                            // histCloneMC->Scale(1./histCloneMC->Integral());
                        }
                        else continue;

                        if (HistExistMap.find(runNumberNameMC) == HistExistMap.end()) {
                            runListMC = nullptr;
                            runListMC = new TList;
                            runListMC->SetOwner(kTRUE);
                            runListMC->SetName(Form("%09d", runNumberNameMC.Atoi()));
                            HistExistMap[runNumberNameMC] = runListMC;
                        }
                        HistExistMap[runNumberNameMC]->Add(histCloneMC);

                    }
                    else continue;
                }

            }

        }
    }

    file->cd();

    for (auto it = HistExistMap.begin(); it != HistExistMap.end(); it++){
        ccdb_object->Add(it->second);
    }

    ccdb_object->Write(ccdb_object->GetName(), TObject::kSingleKey + TObject::kOverwrite);
    file->Write();

    file->Close();
    delete file;
    file = nullptr;

    file_real_weight->Close();
    file_MC_weight->Close();
    delete file_real_weight;
    delete file_MC_weight;

    return 0;
}