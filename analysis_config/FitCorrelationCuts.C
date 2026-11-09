#include <TCanvas.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResultPtr.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TList.h>
#include <TPad.h>
#include <TString.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>



struct SliceInfo {
  double x = 0.;
  double median = 0.;
  double sigma = 0.;
};

struct CutParameters {
  double a = 0.; // upper slope
  double b = 0.; // upper intercept
  double c = 0.; // lower slope
  double d = 0.; // lower intercept
  int nSlices = 0;
  bool valid = false;
};

// For every x bin, project the TH2F onto y and calculate
// median = Q50 and sigma = (Q84 - Q16) / 2.
std::vector<SliceInfo> ProjectAndCalculate(TH2F* histogram, double minEntriesPerSlice = 1000.)
{
  std::vector<SliceInfo> slices;
  if (!histogram) {
    return slices;
  }

  static int projectionCounter = 0;
  const double probabilities[3] = {0.16, 0.50, 0.84};
  const float cutTail = 0.75;

  for (int xBin = 1; xBin <= cutTail * histogram->GetNbinsX(); ++xBin) {
    const TString projectionName =
      TString::Format("projectionY_%d", projectionCounter++);
    TH1D* projection =
      histogram->ProjectionY(projectionName.Data(), xBin, xBin);
    projection->SetDirectory(nullptr);

    if (projection->Integral() < minEntriesPerSlice) {
      delete projection;
      continue;
    }

    double quantiles[3] = {0., 0., 0.};
    projection->GetQuantiles(3, quantiles, probabilities);

    SliceInfo slice;
    slice.x = histogram->GetXaxis()->GetBinCenter(xBin);
    slice.median = quantiles[1];
    slice.sigma = 0.5 * (quantiles[2] - quantiles[0]);

    if (slice.sigma > 0.) {
      slices.push_back(slice);
      // printf("%f, ", quantiles[2]);
    }
    delete projection;
  }
  // printf("\n\n");

  return slices;
}

// Construct median +/- nSigma*sigma for every slice and fit the upper and
// lower points separately with y = slope*x + intercept.
CutParameters FindCutParameters(const TH2F* histogram, const std::vector<SliceInfo>& slices,
                                double nSigma = 4.)
{
  CutParameters result;
  result.nSlices = static_cast<int>(slices.size());
  if (slices.size() < 2 || nSigma <= 0.) {
    return result;
  }

  TGraph upperGraph(result.nSlices);
  TGraph lowerGraph(result.nSlices);

  for (int i = 0; i < result.nSlices; ++i) {
    const double upper = slices[i].median + nSigma * slices[i].sigma;
    const double lower = slices[i].median - nSigma * slices[i].sigma;
    upperGraph.SetPoint(i, slices[i].x, upper);
    lowerGraph.SetPoint(i, slices[i].x, lower);
  }

  static int fitCounter = 0;
  const int currentFit = fitCounter++;
  const TString upperFitName = TString::Format("upperFit_%d", currentFit);
  const TString lowerFitName = TString::Format("lowerFit_%d", currentFit);
  TF1 upperFit(upperFitName.Data(), "pol1", slices.front().x, slices.back().x);
  TF1 lowerFit(lowerFitName.Data(), "pol1", slices.front().x, slices.back().x);

  const int upperStatus = upperGraph.Fit(&upperFit, "QNR");
  const int lowerStatus = lowerGraph.Fit(&lowerFit, "QNR");
  if (upperStatus != 0 || lowerStatus != 0) {
    return result;
  }

  // For ROOT pol1: parameter 0 is the intercept and parameter 1 is the slope.
  result.a = upperFit.GetParameter(1);
  result.b = upperFit.GetParameter(0);
  result.c = lowerFit.GetParameter(1);
  result.d = lowerFit.GetParameter(0);

  // const double xMin = histogram->GetXaxis()->GetXmin();
  // const double xMax = histogram->GetXaxis()->GetXmax();
  const double xMin = slices.front().x;
  const double xMax = slices.back().x;
  const bool correctOrderAtMin = result.a * xMin + result.b >
                                 result.c * xMin + result.d;
  const bool correctOrderAtMax = result.a * xMax + result.b >
                                 result.c * xMax + result.d;
  result.valid = correctOrderAtMin && correctOrderAtMax;
  return result;
}

void RunAll(const char* realMC, const char* ba, const char* estm, const char* id) {

  const char* filePath = Form("/scratch3/go52dab/O2tutorial/test20260906_generate_corr_cut/%s/AnalysisResults.root", realMC);
  const char* outputPath = Form("CentMultCorrCutConfig/%s/%s/CorrelationCuts.csv", realMC, estm);
  const char* jsonPath = Form("CentMultCorrCutConfig/%s/%s/CorrelationCuts.json", realMC, estm);
  const char* plotPath = Form("CentMultCorrCutConfig/%s/%s/CorrelationCuts.pdf", realMC, estm);

  const double minEntriesPerSlice = 1000.;
  const double nSigma = 4.;

  TFile* analysisFile = TFile::Open(filePath, "READ");
  if (!analysisFile || analysisFile->IsZombie()) {
    std::cerr << "Cannot open " << filePath << '\n';
    return;
  }

  TDirectory* dir =
    dynamic_cast<TDirectory*>(analysisFile->Get(Form("multiparticle-cumulants%s", id)));
  if (!dir) {
    std::cerr << "Cannot find multiparticle-cumulants\n";
    analysisFile->Close();
    delete analysisFile;
    return;
  }

  TList* baseList = nullptr;
  dir->GetObject("Default list name", baseList);
  if (!baseList) {
    std::cerr << "Cannot find Default list name\n";
    analysisFile->Close();
    delete analysisFile;
    return;
  }

  TList* QAHistogramsList =
    dynamic_cast<TList*>(baseList->FindObject("QualityAssuranceHistograms"));
  TList* CorrHistogramsList =
    dynamic_cast<TList*>(baseList->FindObject("CorrelationHistograms"));

  if (!QAHistogramsList) {
    std::cerr << "Cannot find CorrelationHistograms\n";
    analysisFile->Close();
    delete analysisFile;
    return;
  }

  if (!CorrHistogramsList) {
    std::cerr << "Cannot find QualityAssuranceHistograms\n";
    analysisFile->Close();
    delete analysisFile;
    return;
  }

  TH2F* centAB = dynamic_cast<TH2F*>(CorrHistogramsList->FindObject(Form("fHistCentrality[FT0C vs. FT0M][%s cut]", ba)));
  TH2F* centAC = dynamic_cast<TH2F*>(CorrHistogramsList->FindObject(Form("fHistCentrality[FT0C vs. FV0A][%s cut]", ba)));
  TH2F* centBC = dynamic_cast<TH2F*>(CorrHistogramsList->FindObject(Form("fHistCentrality[FT0M vs. FV0A][%s cut]", ba)));
  TH2F* multAB = dynamic_cast<TH2F*>(CorrHistogramsList->FindObject(Form("fHistMultiplicity[FT0C vs. FT0M][%s cut]", ba)));
  TH2F* multAC = dynamic_cast<TH2F*>(CorrHistogramsList->FindObject(Form("fHistMultiplicity[FT0C vs. FV0A][%s cut]", ba)));
  TH2F* multBC = dynamic_cast<TH2F*>(CorrHistogramsList->FindObject(Form("fHistMultiplicity[FT0M vs. FV0A][%s cut]", ba)));
  TH2F* multNcontrib = dynamic_cast<TH2F*>(QAHistogramsList->FindObject(Form("fHistMultiplicity[%s cut]", ba)));

  struct HistogramInput {
    const char* label;
    const char* configKey;
    TH2F* histogram;
  };

  const HistogramInput inputs[] = {
    {"centAB", "cfCentABCorrCut", centAB},
    {"centAC", "cfCentACCorrCut", centAC},
    {"centBC", "cfCentBCCorrCut", centBC},
    {"multAB", "cfMultABCorrCut", multAB},
    {"multAC", "cfMultACCorrCut", multAC},
    {"multBC", "cfMultBCCorrCut", multBC},
    {"multNcontrib", "cfNumContribMultCorrCut", multNcontrib},
  };
  constexpr int nInputs = sizeof(inputs) / sizeof(inputs[0]);
  CutParameters cutResults[nInputs];

  std::ofstream output(outputPath);
  if (!output.is_open()) {
    std::cerr << "Cannot create " << outputPath << '\n';
    analysisFile->Close();
    delete analysisFile;
    return;
  }

  output << "label,a,b,c,d,nSlices,nSigma\n";
  output << std::setprecision(10);

  for (int i = 0; i < nInputs; ++i) {
    const auto& input = inputs[i];
    if (!input.histogram) {
      std::cerr << "Cannot find histogram: " << input.label << '\n';
      continue;
    }

    const std::vector<SliceInfo> slices =
      ProjectAndCalculate(input.histogram, minEntriesPerSlice);
    const CutParameters cuts = FindCutParameters(input.histogram, slices, nSigma);
    cutResults[i] = cuts;

    if (!cuts.valid) {
      std::cerr << "Cannot determine valid cuts for " << input.label << '\n';
      continue;
    }

    output << input.label << ','
           << cuts.a << ',' << cuts.b << ','
           << cuts.c << ',' << cuts.d << ','
           << cuts.nSlices << ',' << nSigma << '\n';

    std::cout << input.label << ": "
              << "a=" << cuts.a << ", b=" << cuts.b
              << ", c=" << cuts.c << ", d=" << cuts.d << '\n';
  }

  output.close();

  bool allCutsValid = true;
  for (const auto& cuts : cutResults) {
    if (!cuts.valid) {
      allCutsValid = false;
      break;
    }
  }

  if (allCutsValid && strcmp(ba, "before") == 0) {
    std::ofstream jsonOutput(jsonPath);
    if (!jsonOutput.is_open()) {
      std::cerr << "Cannot create " << jsonPath << '\n';
      analysisFile->Close();
      delete analysisFile;
      return;
    }

    jsonOutput << std::setprecision(10);
    jsonOutput << "{\n";
    jsonOutput << "  \"multiparticle-cumulants\": {\n";
    for (int i = 0; i < nInputs; ++i) {
      const auto& cuts = cutResults[i];
      jsonOutput << "    \"" << inputs[i].configKey
                 << "\" : {\"values\": [\""
                 << cuts.a << "\", \"" << cuts.b << "\", \""
                 << cuts.c << "\", \"" << cuts.d << "\"]}";
      if (i + 1 < nInputs) {
        jsonOutput << ',';
      }
      jsonOutput << '\n';
    }
    jsonOutput << "  }\n";
    jsonOutput << "}\n";
    jsonOutput.close();
  } else {
    std::cerr << "JSON was not written\n"; // because one or more cut fits are invalid
  }

  // Draw one histogram per page and use that TH2F's full x/y ranges.
  TCanvas canvas("correlationCutsCanvas", "Correlation cuts", 1000, 800);
  const TString openPdf = TString::Format("%s[", plotPath);
  canvas.Print(openPdf.Data());

  for (int i = 0; i < nInputs; ++i) {
    const auto& input = inputs[i];
    const auto& cuts = cutResults[i];
    if (!input.histogram || !cuts.valid) {
      continue;
    }

    canvas.Clear();
    canvas.cd();
    gPad->SetRightMargin(0.15);
    gPad->SetLogz();

    input.histogram->SetStats(false);
    input.histogram->GetXaxis()->SetRange(1, input.histogram->GetNbinsX());
    input.histogram->GetYaxis()->SetRange(1, input.histogram->GetNbinsY());
    input.histogram->Draw("COLZ");
    input.histogram->SetTitle("");

    const double xMin = input.histogram->GetXaxis()->GetXmin();
    const double xMax = input.histogram->GetXaxis()->GetXmax();

    const TString upperName = TString::Format("upper_%d", i);
    const TString lowerName = TString::Format("lower_%d", i);
    TF1 upperLine(upperName.Data(), "[0]*x+[1]", xMin, xMax);
    TF1 lowerLine(lowerName.Data(), "[0]*x+[1]", xMin, xMax);
    upperLine.SetParameters(cuts.a, cuts.b);
    lowerLine.SetParameters(cuts.c, cuts.d);
    upperLine.SetLineColor(kRed + 1);
    lowerLine.SetLineColor(kBlack);
    upperLine.SetLineWidth(3);
    lowerLine.SetLineWidth(3);
    upperLine.Draw("SAME");
    lowerLine.Draw("SAME");

    TLegend legend(0.15, 0.75, 0.42, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(&upperLine, "upper: a x + b", "l");
    legend.AddEntry(&lowerLine, "lower: c x + d", "l");
    // legend.Draw();

    canvas.Modified();
    canvas.Update();
    canvas.Print(plotPath);

    const TString pngPath = TString::Format("CentMultCorrCutConfig/%s/%s/chap3_CorrelationCuts_%s_%s.png", realMC, estm, input.label, ba);
    canvas.SaveAs(pngPath.Data());
  }

  const TString closePdf = TString::Format("%s]", plotPath);
  canvas.Print(closePdf.Data());
  analysisFile->Close();
  delete analysisFile;

  std::cout << "Results saved to " << outputPath << " and " << plotPath;
  if (allCutsValid) {
    std::cout << " and " << jsonPath;
  }
  std::cout << '\n';

}

int FitCorrelationCuts()
{
  // RunAll("realdata", "before", "multFT0M", "");
  // RunAll("realdata", "before", "multFV0A", "_id62521");
  // RunAll("realdata", "before", "multFT0C", "_id62520");

  RunAll("MCdata", "before", "multFT0M", "_id62802");
  RunAll("MCdata", "before", "multFV0A", "_id62521");
  RunAll("MCdata", "before", "multFT0C", "_id62520");

  return 0;
}
