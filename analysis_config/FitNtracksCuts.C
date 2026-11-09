#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLine.h"
#include "TList.h"
#include "TString.h"
#include "TSystem.h"

struct FitParam {
  double A;
  double mean;
  double sig;
  double fitMin;
  double fitMax;
  double chi2;
  int ndf;
};

FitParam Fit(TH1F* nTracksHist) {
  if (!nTracksHist || nTracksHist->GetEntries() < 100) {
    throw std::runtime_error("Missing or insufficient NumTrack entries");
  }

  const double peak = nTracksHist->GetBinCenter(nTracksHist->GetMaximumBin());
  const double rms = nTracksHist->GetRMS();
  if (!std::isfinite(rms) || rms <= 0.) {
    throw std::runtime_error("Invalid NumTrack RMS");
  }

  // Fit the core; tails outside this range do not determine mu and sigma.
  const double fitMin = std::max(nTracksHist->GetXaxis()->GetXmin(), peak - 2. * rms);
  const double fitMax = std::min(nTracksHist->GetXaxis()->GetXmax(), peak + 2. * rms);
  TF1 nTracksFit("nTracksFit", "gaus", fitMin, fitMax);
  nTracksFit.SetParameters(nTracksHist->GetMaximum(), peak, rms);

  // N: do not store a second TF1 in the histogram; S: inspect fit status.
  const auto status = nTracksHist->Fit(&nTracksFit, "RQSN");
  if (int(status) != 0 || !status.Get() || !status->IsValid() ||
      !std::isfinite(nTracksFit.GetParameter(1)) ||
      !std::isfinite(nTracksFit.GetParameter(2)) ||
      nTracksFit.GetParameter(2) <= 0.) {
    throw std::runtime_error("Gaussian fit failed");
  }

  return {nTracksFit.GetParameter(0), nTracksFit.GetParameter(1),
          nTracksFit.GetParameter(2), fitMin, fitMax,
          nTracksFit.GetChisquare(), nTracksFit.GetNDF()};
}

void RunAll(const char* realMC) {

  constexpr int nCentBins = 8;
  constexpr double nSigma = 3.;
  const char* centBinsNames[nCentBins] = {
      "00-05", "05-10", "10-20", "20-30",
      "30-40", "40-50", "50-60", "60-70"};

  const TString filePath = TString::Format(
      "/scratch3/go52dab/O2tutorial/test20260922_analysis/%s/AnalysisResults.root",
      realMC);
  TFile* analysisFile = TFile::Open(filePath.Data(), "READ");
  if (!analysisFile || analysisFile->IsZombie()) {
    std::cerr << "Cannot open " << filePath << '\n';
    delete analysisFile;
    return;
  }

  for (const char* centBinsName : centBinsNames) {
    TDirectory* dir = analysisFile->GetDirectory(
        TString::Format("multiparticle-cumulants_cent%s", centBinsName));
    if (!dir) {
      std::cerr << "Cannot find multiparticle-cumulants" << centBinsName << '\n';
      continue;
    }

    TList* baseList = nullptr;
    dir->GetObject("Default list name", baseList);
    if (!baseList) {
      std::cerr << "Cannot find Default list name for " << centBinsName << '\n';
      continue;
    }

    auto* eventHistograms =
        dynamic_cast<TList*>(baseList->FindObject("EventHistograms"));
    if (!eventHistograms) {
      std::cerr << "Cannot find EventHistograms for " << centBinsName << '\n';
      continue;
    }

    auto* nTracksHist = dynamic_cast<TH1F*>(
        eventHistograms->FindObject("fHistNumTrack[eRec][after cut]"));
    if (!nTracksHist) {
      std::cerr << "Cannot find fHistNumTrack[eRec][after cut] for "
                << centBinsName << '\n';
      continue;
    }

    FitParam fit;
    try {
      fit = Fit(nTracksHist);
    } catch (const std::exception& error) {
      std::cerr << "Skipping " << centBinsName << ": " << error.what() << '\n';
      continue;
    }

    const double lower = std::max(0., fit.mean - nSigma * fit.sig);
    const double upper = fit.mean + nSigma * fit.sig;
    const TString outputDir = TString::Format(
        "NtracksCutConfig/%s/%s",
        realMC, centBinsName);
    if (gSystem->AccessPathName(outputDir.Data()) &&
        gSystem->mkdir(outputDir.Data(), true) != 0) {
      std::cerr << "Cannot create " << outputDir << '\n';
      continue;
    }

    const TString csvPath = outputDir + "/CorrelationCuts.csv";
    const TString jsonPath = outputDir + "/CorrelationCuts.json";
    const TString plotPath = outputDir + "/CorrelationCuts.png";

    std::ofstream csv(csvPath.Data());
    std::ofstream json(jsonPath.Data());
    if (!csv || !json) {
      std::cerr << "Cannot write CSV or JSON in " << outputDir << '\n';
      continue;
    }

    csv << "centrality,entries,A,mean,sigma,nSigma,lower,upper,chi2,ndf,fitMin,fitMax\n";
    csv << centBinsName << ',' << nTracksHist->GetEntries() << ','
        << std::fixed << std::setprecision(6)
        << fit.A << ',' << fit.mean << ',' << fit.sig << ',' << nSigma
        << ',' << lower << ',' << upper << ',' << fit.chi2 << ','
        << fit.ndf << ',' << fit.fitMin << ',' << fit.fitMax << '\n';

    json << "{\n"
         << "  \"multiparticle-cumulants\": {\n"
         << "    \"cfNumTrackCut\": {\"values\": [\""
         << std::fixed << std::setprecision(3) << lower << "\", \""
         << upper << "\"]}\n"
         << "  }\n"
         << "}\n";
    csv.close();
    json.close();

    gStyle->SetOptStat(0);
    TCanvas canvas(TString::Format("cNumTrack_%s", centBinsName),
                   "NumTrack cut", 1100, 800);
    // nTracksHist->SetTitle(TString::Format(
    //     "NumTrack distribution, centrality %s%%;NumTrack;Events", centBinsName));
    nTracksHist->SetTitle("");
    nTracksHist->Draw("HIST");

    TF1 drawFit("drawFit", "gaus", fit.fitMin, fit.fitMax);
    drawFit.SetParameters(fit.A, fit.mean, fit.sig);
    drawFit.SetLineColor(kRed + 1);
    drawFit.SetLineWidth(2);
    drawFit.Draw("SAME");

    const double yMax = nTracksHist->GetMaximum() * 1.05;
    TLine lowerLine(lower, 0., lower, yMax);
    TLine upperLine(upper, 0., upper, yMax);
    lowerLine.SetLineColor(kBlue + 1);
    upperLine.SetLineColor(kBlue + 1);
    lowerLine.SetLineStyle(2);
    upperLine.SetLineStyle(2);
    lowerLine.SetLineWidth(2);
    upperLine.SetLineWidth(2);
    lowerLine.Draw("SAME");
    upperLine.Draw("SAME");

    TLegend legend(0.60, 0.73, 0.89, 0.89);
    legend.AddEntry(nTracksHist, "Events", "l");
    legend.AddEntry(&drawFit, "Gaussian fit", "l");
    legend.AddEntry(&lowerLine, Form("#pm %d#sigma cut", static_cast<int>(nSigma)), "l");
    legend.Draw();
    canvas.SaveAs(plotPath.Data());

    std::cout << centBinsName << "%: mean=" << fit.mean
              << ", sigma=" << fit.sig << ", cut=["
              << lower << ", " << upper << "]\n";
  }

  analysisFile->Close();
  delete analysisFile;
}


void FitNtracksCuts() {
  RunAll("realdata");
  RunAll("MCdata");
}
