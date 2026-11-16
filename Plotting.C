#include <cmath>
#include <iostream>
#include <vector>
#include <limits>

#include "TCanvas.h"
#include "TColor.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TObject.h"
#include "TString.h"
#include "TSystem.h"

#include "ShiftAlongXaxis.h"
#include "loadtable.h"

using namespace std;

const bool bLocal = false;
bool bNormalized = true;
bool bCumulant = false;
bool bShowMC = false;

const char* resultPath_vnscVScent_real = "result/data/result_vnscVScent_realdata.csv";
const char* resultPath_vnscVScent_MC = "result/data/result_vnscVScent_MCdata.csv";
const char* resultPath_trending = "result/data/result_trending_realdata.csv";
const char* resultPath_root_cutchecks = "/scratch3/go52dab/O2tutorial/test20260922_analysis/realdata/AnalysisResults.root";
const char* resultPath_weights = "/scratch3/go52dab/O2tutorial/analysis_code/weights.root";

const int nRunNumberMC = 4;
const int nRunNumberReal = 8;
const int runNumberMC[nRunNumberMC] = {544886, 544868, 544122, 544116};
const int runNumberReal[nRunNumberReal] = {544124, 544123, 544122, 544121, 544116, 544098, 544095, 544091};

// tools

double getYMax(TMultiGraph* mg, bool errbar = true){

    // get maximum
    double yMax = -std::numeric_limits<double>::infinity();

    TIter next(mg->GetListOfGraphs());
    TObject* obj = nullptr;

    while ((obj = next())) {
      auto* gr = dynamic_cast<TGraphErrors*>(obj);
      auto gr_arr = gr->GetY();
      auto grerr_arr = gr->GetEY();

      for (int i = 0; i < gr->GetN(); i++) {
        double element_pos = 0.;
        if (errbar) {
          element_pos = gr_arr[i] + grerr_arr[i];
        } else {
          element_pos = gr_arr[i];
        }
        if (element_pos > yMax) {
            yMax = element_pos;
        }
      }
    }
    return yMax;
}

double getYMin(TMultiGraph* mg, bool errbar = true){

    // get minimum
    double yMin = std::numeric_limits<double>::infinity();

    TIter next(mg->GetListOfGraphs());
    TObject* obj = nullptr;

    while ((obj = next())) {
      auto* gr = dynamic_cast<TGraphErrors*>(obj);
      auto gr_arr = gr->GetY();
      auto grerr_arr = gr->GetEY();

      for (int i = 0; i < gr->GetN(); i++) {
        double element_neg = 0.;
        if (errbar) {
          element_neg = gr_arr[i] - grerr_arr[i];
        } else {
          element_neg = gr_arr[i];
        }
        if (element_neg < yMin) {
            yMin = element_neg;
        }
      }
    }
    return yMin;
}

void setRunPositions(TGraphErrors* graph, int nRunNumber)
{
    for (int i = 0; i < nRunNumber; i++) {
        double oldX = 0.;
        double y = 0.;
        graph->GetPoint(i, oldX, y);

        const double errorY = graph->GetErrorY(i);

        graph->SetPoint(i, i + 1, y);
        graph->SetPointError(i, 0., errorY);
    }
}

void setRunLabels(TMultiGraph* multigraph, const int runNumber[], int nRunNumber)
{
    TAxis* axis = multigraph->GetXaxis();

    axis->Set(nRunNumber, 0.5, nRunNumber + 0.5);

    for (int i = 0; i < nRunNumber; i++) {
        axis->SetBinLabel(i + 1, Form("%d", runNumber[i]));
    }

    axis->LabelsOption("v");
    axis->SetTitle("Run number");
}

void convertToRelativeDeviation(TGraphErrors* graph)
{
    const int nPoints = graph->GetN();

    double mean = 0.;
    for (int i = 0; i < nPoints; i++) {
        double x = 0.;
        double y = 0.;
        graph->GetPoint(i, x, y);
        mean += y;
    }
    mean /= nPoints;

    for (int i = 0; i < nPoints; i++) {
        double x = 0.;
        double y = 0.;
        graph->GetPoint(i, x, y);

        const double errorY = graph->GetErrorY(i);
        const double relativeValue = 100. * (y - mean) / mean;
        const double relativeError = 100. * errorY / std::abs(mean);

        graph->SetPoint(i, x, relativeValue);
        graph->SetPointError(i, 0., relativeError);
    }
}

// plots functions

void plotFlow(int harmonic, TGraphErrors* twoParticleGraphReal_stat, TGraphErrors* twoParticleGraphMC_stat, TGraphErrors* fourParticleGraphReal_stat, TGraphErrors* fourParticleGraphMC_stat, const std::vector<TGraphErrors*>& dataV2Compare_stat, const std::vector<TGraphErrors*>& dataV3V4Compare_stat,
                            TGraphErrors* twoParticleGraphReal_sys, TGraphErrors* twoParticleGraphMC_sys, TGraphErrors* fourParticleGraphReal_sys, TGraphErrors* fourParticleGraphMC_sys, const std::vector<TGraphErrors*>& dataV2Compare_sys, const std::vector<TGraphErrors*>& dataV3V4Compare_sys, const bool bShowCompare) {


  std::vector<TGraphErrors*> twoParticleGraphReal = {twoParticleGraphReal_sys, twoParticleGraphReal_stat};
  std::vector<TGraphErrors*> twoParticleGraphMC = {twoParticleGraphMC_sys, twoParticleGraphMC_stat};
  std::vector<TGraphErrors*> fourParticleGraphReal = {fourParticleGraphReal_sys, fourParticleGraphReal_stat};
  std::vector<TGraphErrors*> fourParticleGraphMC = {fourParticleGraphMC_sys, fourParticleGraphMC_stat};

  std::vector<std::vector<TGraphErrors*>> dataV2Compare = {dataV2Compare_sys, dataV2Compare_stat};
  std::vector<std::vector<TGraphErrors*>> dataV3V4Compare = {dataV3V4Compare_sys, dataV3V4Compare_stat};

  auto* canvas = new TCanvas(Form("c_v%d", harmonic), "", 1500, 900);
  canvas->SetGrid();
  const float alpha = 0.5;

  for (size_t s = 0; s < 2; s++) {
    twoParticleGraphReal[s]->SetMarkerStyle(20);
    twoParticleGraphReal[s]->SetMarkerColor(kRed);
    twoParticleGraphMC[s]->SetMarkerStyle(20);
    twoParticleGraphMC[s]->SetMarkerColor(kGreen + 2);
    fourParticleGraphReal[s]->SetMarkerStyle(harmonic == 2 ? 22 : 26);
    fourParticleGraphReal[s]->SetMarkerColor(kRed);
    fourParticleGraphMC[s]->SetMarkerStyle(harmonic == 2 ? 22 : 26);
    fourParticleGraphMC[s]->SetMarkerColor(kGreen + 2);
    if (s == 0) {
      twoParticleGraphReal[s]->SetFillColorAlpha(kRed - 9, alpha);
      twoParticleGraphMC[s]->SetFillColorAlpha(kGreen - 9, alpha);
      fourParticleGraphReal[s]->SetFillColorAlpha(kRed - 9, alpha);
      fourParticleGraphMC[s]->SetFillColorAlpha(kGreen - 9, alpha);
      twoParticleGraphReal[s]->SetLineWidth(0); 
      twoParticleGraphMC[s]->SetLineWidth(0); 
      fourParticleGraphReal[s]->SetLineWidth(0); 
      fourParticleGraphMC[s]->SetLineWidth(0); 
    }

    if (bShowCompare) {
      if (harmonic == 2) {
        const int markers[4] = {20, 22, 25, 32};
        for (int i = 0; i < 4; ++i) {
          dataV2Compare[s][i]->SetMarkerStyle(markers[i]);
          dataV2Compare[s][i]->SetMarkerColor(kBlue);
          if (s == 0) {
            dataV2Compare[s][i]->SetFillColorAlpha(kBlue - 9, alpha);
            dataV2Compare[s][i]->SetLineWidth(0); 
          }
        }
      } else {
        const int index = harmonic - 3;
        dataV3V4Compare[s][index]->SetMarkerStyle(20);
        dataV3V4Compare[s][index]->SetMarkerColor(kBlue);
        if (s == 0) {
          dataV3V4Compare[s][index]->SetFillColorAlpha(kBlue - 9, alpha);
          dataV3V4Compare[s][index]->SetLineWidth(0); 
        } 
      }
    }

  }

  auto* multigraph = new TMultiGraph();
  multigraph->SetTitle(" ");

  const char* drawOp[2] = {"2", "P"};
  for (int s = 0; s < 2; s++) {
    multigraph->Add(twoParticleGraphReal[s], drawOp[s]);
    multigraph->Add(fourParticleGraphReal[s], drawOp[s]);
    multigraph->Add(twoParticleGraphMC[s], drawOp[s]);
    multigraph->Add(fourParticleGraphMC[s], drawOp[s]);

    if (bShowCompare) {
      if (harmonic == 2) {
        for (auto* graph : dataV2Compare[s]) {
          multigraph->Add(graph, drawOp[s]);
        }
      } else {
        multigraph->Add(dataV3V4Compare[s][harmonic - 3], drawOp[s]);
      }
    }
  }
  
  TIter next(multigraph->GetListOfGraphs());
  int multigraphSize = multigraph->GetListOfGraphs()->GetSize();
  int graphIndex = 0;
  int graphIndexTrue = 0;
  if (harmonic != 4) {
    while (auto* object = next()) {
      auto* graph = dynamic_cast<TGraphErrors*>(object);
      if (graphIndexTrue >= multigraphSize / 2) {
        graphIndexTrue = 0;
        graphIndex = 0;
      }
      // if (!graph) {
      //   continue;
      // }
      if (graphIndexTrue == 2) {
        ShiftAlongXaxis(graph, 1);
        ++graphIndexTrue;
        continue;
      }
      else if (graphIndexTrue == 3) {
        ShiftAlongXaxis(graph, 2);
        ++graphIndexTrue;
        continue;
      }
      graph->SetLineColor(kBlack);
      ShiftAlongXaxis(graph, - 0.6 + 0.7 * graphIndex);
      ++graphIndexTrue;
      ++graphIndex;
    }
  }
  else {
    while (auto* object = next()) {
      auto* graph = dynamic_cast<TGraphErrors*>(object);
      if (graphIndexTrue >= multigraphSize / 2) {
        graphIndexTrue = 0;
        graphIndex = 0;
      }
      // if (!graph) {
      //   continue;
      // }
      graph->SetLineColor(kBlack);
      ShiftAlongXaxis(graph, - 0.6 + 0.7 * graphIndex);
      ++graphIndexTrue;
      ++graphIndex;
    }

  }

  TLegend* legend = nullptr;
  const double yMax = getYMax(multigraph);
  const double yMin = getYMin(multigraph);
  if (bShowCompare) {
    if (yMax > 0. && yMin >= 0.) {
      legend = harmonic == 2 ? new TLegend(0.40, 0.65, 0.85, 0.85)
                             : new TLegend(0.40, 0.73, 0.85, 0.85);
      multigraph->SetMinimum(0.);
      multigraph->SetMaximum(1.7 * yMax);
    } else if (yMax > 0. && yMin < 0. && std::abs(yMax) < std::abs(yMin)) {
      legend = harmonic == 2 ? new TLegend(0.15, 0.15, 0.60, 0.15 + 0.04 * graphIndex)
                             : new TLegend(0.15, 0.15, 0.60, 0.27);
      multigraph->SetMinimum(1.7 * yMax);
    } else if (yMax > 0. && yMin < 0. && std::abs(yMax) > std::abs(yMin)) {
      legend = harmonic == 2 ? new TLegend(0.40, 0.65, 0.85, 0.85)
                             : new TLegend(0.40, 0.73, 0.85, 0.85);
      multigraph->SetMinimum(1.7 * yMin);
      multigraph->SetMaximum(1.7 * yMax);
    } else {
      legend = harmonic == 2 ? new TLegend(0.15, 0.15, 0.60, 0.15 + 0.04 * graphIndex)
                             : new TLegend(0.15, 0.15, 0.60, 0.27);
      multigraph->SetMinimum(1.7 * yMin);
      multigraph->SetMaximum(0);
    }


  } else {
    legend = new TLegend(0.55, 0.73, 0.85, 0.85);
    multigraph->SetMinimum(0.);
    multigraph->SetMaximum(1.7 * yMax);
  }

  multigraph->Draw("A");
  multigraph->GetXaxis()->SetLimits(0., harmonic == 2 && bShowCompare ? 80. : 70.);
  multigraph->GetXaxis()->SetTitle("Centrality (\%)");
  multigraph->GetYaxis()->SetTitle(Form("v_{%d}", harmonic));
  multigraph->GetYaxis()->SetTitleOffset(1);

  if (harmonic == 2) {
    if (bShowCompare) {
      legend->SetNColumns(3);

      legend->AddEntry((TObject*)0, "PYTHIA8", "");
      legend->AddEntry((TObject*)0, "ALICE Pb-Pb", "");
      legend->AddEntry((TObject*)0, "", "");

      legend->AddEntry((TObject*)0, "5.36 TeV", "");
      legend->AddEntry((TObject*)0, "5.36 TeV", "");
      legend->AddEntry((TObject*)0, "5.02 TeV", "");

      legend->AddEntry(twoParticleGraphMC[0], "v_{2}{2, |#Delta#eta|>1}");
      legend->AddEntry(twoParticleGraphReal[0], "v_{2}{2, |#Delta#eta|>1}");
      legend->AddEntry(dataV2Compare[0][0], "v_{2}{2, |#Delta#eta|>1}");

      legend->AddEntry(fourParticleGraphMC[0], "v_{2}{4}");
      legend->AddEntry(fourParticleGraphReal[0], "v_{2}{4}");
      legend->AddEntry(dataV2Compare[0][1], "v_{2}{4}");

      legend->AddEntry((TObject*)0, "", "");
      legend->AddEntry((TObject*)0, "", "");
      legend->AddEntry(dataV2Compare[0][2], "v_{2}{6}");

      legend->AddEntry((TObject*)0, "", "");
      legend->AddEntry((TObject*)0, "", "");
      legend->AddEntry(dataV2Compare[0][3], "v_{2}{8}");
    }
    else {
      legend->SetNColumns(2);

      legend->AddEntry((TObject*)0, "PYTHIA8", "");
      legend->AddEntry((TObject*)0, "ALICE Pb-Pb", "");

      legend->AddEntry((TObject*)0, "5.36 TeV", "");
      legend->AddEntry((TObject*)0, "5.36 TeV", "");

      legend->AddEntry(twoParticleGraphMC[0], "v_{2}{2}");
      legend->AddEntry(twoParticleGraphReal[0], "v_{2}{2}");

      legend->AddEntry(fourParticleGraphMC[0], "v_{2}{4}");
      legend->AddEntry(fourParticleGraphReal[0], "v_{2}{4}");

    }

  } else {

    if (bShowCompare) {
      legend->SetNColumns(3);

      legend->AddEntry((TObject*)0, "PYTHIA8", "");
      legend->AddEntry((TObject*)0, "ALICE Pb-Pb", "");
      legend->AddEntry((TObject*)0, "", "");

      legend->AddEntry((TObject*)0, "5.36 TeV", "");
      legend->AddEntry((TObject*)0, "5.36 TeV", "");
      legend->AddEntry((TObject*)0, "5.02 TeV", "");

      legend->AddEntry(twoParticleGraphMC[0], Form("v_{%d}{2, |#Delta#eta|>1}", harmonic));
      legend->AddEntry(twoParticleGraphReal[0], Form("v_{%d}{2, |#Delta#eta|>1}", harmonic));
      legend->AddEntry(dataV3V4Compare[0][harmonic - 3], Form("v_{%d}{2, |#Delta#eta|>1}", harmonic));

      legend->AddEntry(fourParticleGraphMC[0], Form("v_{%d}{4}", harmonic));
      legend->AddEntry(fourParticleGraphReal[0], Form("v_{%d}{4}", harmonic));
      legend->AddEntry((TObject*)0, "", "");

    } else {
      legend->SetNColumns(2);

      legend->AddEntry((TObject*)0, "PYTHIA8", "");
      legend->AddEntry((TObject*)0, "ALICE Pb-Pb", "");

      legend->AddEntry((TObject*)0, "5.36 TeV", "");
      legend->AddEntry((TObject*)0, "5.36 TeV", "");

      legend->AddEntry(twoParticleGraphMC[0], Form("v_{%d}{2}", harmonic));
      legend->AddEntry(twoParticleGraphReal[0], Form("v_{%d}{2}", harmonic));

      legend->AddEntry(fourParticleGraphMC[0], Form("v_{%d}{4}", harmonic));
      legend->AddEntry(fourParticleGraphReal[0], Form("v_{%d}{4}", harmonic));
    }

  }

  legend->SetBorderSize(0);
  legend->Draw();
  if (bShowCompare) {
    canvas->Print(Form("result/plots/chap4_result_v%d.pdf", harmonic));
    canvas->Print(Form("result/plots/chap4_result_v%d.png", harmonic));
  }
  else {
    canvas->Print(Form("result/plots/chap4_result_v%d_wocompare.pdf", harmonic));
    canvas->Print(Form("result/plots/chap4_result_v%d_wocompare.png", harmonic));
  }

  delete canvas;
  // delete multigraph;
}


void plotFlow2pGap(int harmonic, TGraphErrors* twoParticleGap0, TGraphErrors* twoParticleGap1, TGraphErrors* twoParticleGap2, TGraphErrors* twoParticleGap3) {

  auto* canvas = new TCanvas(Form("c_v%d", harmonic), "", 1500, 900);
  canvas->SetGrid();

  twoParticleGap0->SetMarkerStyle(20);
  twoParticleGap0->SetMarkerColor(kRed);
  twoParticleGap1->SetMarkerStyle(20);
  twoParticleGap1->SetMarkerColor(kBlue);
  twoParticleGap2->SetMarkerStyle(20);
  twoParticleGap2->SetMarkerColor(kGreen + 2);
  twoParticleGap3->SetMarkerStyle(20);
  twoParticleGap3->SetMarkerColor(kYellow + 2);

  auto* multigraph = new TMultiGraph();
  multigraph->SetTitle(" ");
  multigraph->Add(twoParticleGap0);
  multigraph->Add(twoParticleGap1);
  multigraph->Add(twoParticleGap2);
  multigraph->Add(twoParticleGap3);

  TIter next(multigraph->GetListOfGraphs());
  int graphIndex = 0;
  while (auto* object = next()) {
    auto* graph = dynamic_cast<TGraphErrors*>(object);
    if (!graph) {
      continue;
    }
    graph->SetLineColor(kBlack);
    ShiftAlongXaxis(graph, 0.7*graphIndex);
    // printf("index:%d\n", graphIndex);
    ++graphIndex;
  }

  TLegend* legend = nullptr;
  const double yMax = getYMax(multigraph);
  legend = new TLegend(0.65, 0.67, 0.85, 0.85);
  multigraph->SetMinimum(0.);
  multigraph->SetMaximum(2. * yMax);

  multigraph->Draw("ap");
  multigraph->GetXaxis()->SetLimits(0., 70.);
  multigraph->GetXaxis()->SetTitle("Centrality (\%)");
  multigraph->GetYaxis()->SetTitle(Form("v_{%d}", harmonic));

  legend->AddEntry((TObject*)0, "ALICE Pb-Pb", "");
  legend->AddEntry((TObject*)0, "5.36 TeV", "");

  legend->AddEntry(twoParticleGap0, Form("v_{%d}{2}", harmonic));
  legend->AddEntry(twoParticleGap1, Form("v_{%d}{2}, |#Delta#eta|>0.4", harmonic));
  legend->AddEntry(twoParticleGap2, Form("v_{%d}{2}, |#Delta#eta|>0.8", harmonic));
  legend->AddEntry(twoParticleGap3, Form("v_{%d}{2}, |#Delta#eta|>1.0", harmonic));

  legend->SetBorderSize(0);
  legend->Draw();
  canvas->Print(Form("result/plots/chap4_result_v%d_etagap.pdf", harmonic));
  canvas->Print(Form("result/plots/chap4_result_v%d_etagap.png", harmonic));
  delete canvas;
  // delete multigraph;
}

void plotNSC(int index, TGraphErrors* resultGraphReal, TGraphErrors* resultGraphMC, TGraphErrors* comparisonGraph, const bool bShowCompare, const bool bShowMC, const bool bNormalized)
{

  auto* canvas = new TCanvas(Form("c_NSC%d", index), "", 1500, 900);
  canvas->SetGrid();

  resultGraphReal->SetMarkerStyle(20);
  resultGraphReal->SetMarkerColor(kRed);

  if (bShowMC) {
  resultGraphMC->SetMarkerStyle(20);
  resultGraphMC->SetMarkerColor(kGreen + 2);
  }

  if (bShowCompare) {
    comparisonGraph->SetMarkerStyle(20);
    comparisonGraph->SetMarkerColor(kBlue);
  }

  auto* multigraph = new TMultiGraph();
  multigraph->SetTitle(" ");
  multigraph->Add(resultGraphReal);
  if (bShowMC) {
    multigraph->Add(resultGraphMC);
  }
  if (bShowCompare) {
    multigraph->Add(comparisonGraph);
  }

  TIter next(multigraph->GetListOfGraphs());
  int graphIndex = 0;
  while (auto* object = next()) {
    auto* graph = dynamic_cast<TGraphErrors*>(object);
    if (!graph) {
      continue;
    }
    graph->SetLineColor(kBlack);
    ShiftAlongXaxis(graph, 0.7 * graphIndex);
    ++graphIndex;
  }

  const double yMax = getYMax(multigraph);
  const double yMin = getYMin(multigraph);
  TLegend* legend = nullptr;

  if (yMax > 0. && yMin > 0.) {
    legend = new TLegend(0.53, 0.73, 0.83, 0.85);
    multigraph->SetMinimum(0.);
    multigraph->SetMaximum(1.5 * yMax);
  } else if (yMax > 0. && yMin < 0. && std::abs(yMax) < std::abs(yMin)) {
    legend = new TLegend(0.15, 0.15, 0.45, 0.27);
    multigraph->SetMinimum(1.5 * yMin);
    multigraph->SetMaximum(1.5 * yMax);
  } else if (yMax > 0. && yMin < 0. && std::abs(yMax) > std::abs(yMin)) {
    legend = new TLegend(0.53, 0.73, 0.83, 0.85);
    multigraph->SetMinimum(1.5 * yMin);
    multigraph->SetMaximum(1.5 * yMax);
  } else {
    legend = new TLegend(0.15, 0.15, 0.45, 0.27);
    multigraph->SetMinimum(1.5 * yMin);
    multigraph->SetMaximum(0);
  }

  multigraph->Draw("ap");
  multigraph->GetXaxis()->SetLimits(0., 70.);
  multigraph->GetXaxis()->SetTitle("Centrality (\%)");

  TString label;
  if (bNormalized) {
    label = index == 5 ? "NSC(2,3,4)" : Form("NSC(2,%d)", index);
  }
  else {
    label = index == 5 ? "SC(2,3,4)" : Form("SC(2,%d)", index);
  }
  multigraph->GetYaxis()->SetTitle(label.Data());

  if (bShowCompare && bShowMC) {
    legend->SetNColumns(3);

    legend->AddEntry((TObject*)0, "PYTHIA8", "");
    legend->AddEntry((TObject*)0, "ALICE Pb-Pb", "");
    legend->AddEntry((TObject*)0, "", "");

    legend->AddEntry((TObject*)0, "5.36 TeV", "");
    legend->AddEntry((TObject*)0, "5.36 TeV", "");
    legend->AddEntry((TObject*)0, "2.76 TeV", "");

    legend->AddEntry(resultGraphMC, label.Data());
    legend->AddEntry(resultGraphReal, label.Data());
    legend->AddEntry(comparisonGraph, label.Data());

  } else if (!bShowCompare && bShowMC){
    legend->SetNColumns(2);

    legend->AddEntry((TObject*)0, "ALICE Pb-Pb", "");
    legend->AddEntry((TObject*)0, "", "");

    legend->AddEntry((TObject*)0, "5.36 TeV", "");
    legend->AddEntry((TObject*)0, "", "");

    legend->AddEntry(resultGraphMC, label.Data());
    legend->AddEntry(resultGraphReal, label.Data());

  } else if (bShowCompare && !bShowMC){
    legend->SetNColumns(2);

    legend->AddEntry((TObject*)0, "ALICE Pb-Pb", "");
    legend->AddEntry((TObject*)0, "", "");

    legend->AddEntry((TObject*)0, "5.36 TeV", "");
    legend->AddEntry((TObject*)0, "2.76 TeV", "");

    legend->AddEntry(resultGraphReal, label.Data());
    legend->AddEntry(comparisonGraph, label.Data());
  }

  legend->SetBorderSize(0);
  legend->Draw();

  if (bNormalized) {
    if (index == 5) {
      canvas->Print("result/plots/chap4_result_NSC234.pdf");
      canvas->Print("result/plots/chap4_result_NSC234.png");
    } else {
      canvas->Print(Form("result/plots/chap4_result_NSC2%d.pdf", index));
      canvas->Print(Form("result/plots/chap4_result_NSC2%d.png", index));
    }
  }
  else {
    if (index == 5) {
      canvas->Print("result/plots/chap4_result_SC234.pdf");
      canvas->Print("result/plots/chap4_result_SC234.png");
    } else {
      canvas->Print(Form("result/plots/chap4_result_SC2%d.pdf", index));
      canvas->Print(Form("result/plots/chap4_result_SC2%d.png", index));
    }
  }

  delete canvas;
  // delete multigraph;
}


void plotTrending(std::vector<TGraphErrors*> TrendingGraphs)
{
  // const int nRunNumber = 8;
  // const int runNumber[nRunNumber] = {544124, 544123, 544122, 544121, 544116, 544098, 544095, 544091};
  const int nObservables = 9;
  const char* observablesList[nObservables] = {"multiplicity", "Nparticles", "pT", "<<2>>2cent10-20", "<<2>>3cent10-20", "<<2>>4cent10-20", "<<2>>2cent50-60", "<<2>>3cent50-60", "<<2>>4cent50-60"};
  const char* observablesUnitList[nObservables] = {"", "", "[GeV]", "", "", "", "", "", ""};

  TCanvas* canvas = nullptr;
  TMultiGraph* multigraph = nullptr;

  for (int i = 0; i < 4; i++) {

    canvas = new TCanvas("c", "", 1500, 900);
    canvas->SetBottomMargin(0.20);
    canvas->SetGrid();

    if (i != 3) {

      TGraphErrors* Graph = TrendingGraphs[i];
      setRunPositions(Graph, nRunNumberReal);
      // convertToRelativeDeviation(Graph);
      Graph->SetMarkerStyle(20);
      Graph->SetMarkerColor(kRed);

      multigraph = new TMultiGraph();
      multigraph->SetTitle(" ");
      multigraph->Add(Graph);

      const double yRangeMax = std::abs(getYMax(multigraph, false));
      const double yRangeMin = std::abs(getYMin(multigraph, false));
      multigraph->SetMinimum(yRangeMin - 0.2 * std::abs(yRangeMax));
      multigraph->SetMaximum(yRangeMax + 0.2 * std::abs(yRangeMax));

      TLegend* legend = new TLegend(0.55, 0.73, 0.85, 0.85);
      legend->SetBorderSize(0);
      legend->Draw();

      multigraph->Draw("AP");
      setRunLabels(multigraph, runNumberReal, nRunNumberReal);
      multigraph->GetXaxis()->SetTitle("Run number");
      multigraph->GetYaxis()->SetTitle(Form("%s %s", observablesList[i], observablesUnitList[i]));
      multigraph->GetXaxis()->SetTitleOffset(2);


      canvas->Print(Form("result/plots/chap4_result_trending_%s.pdf", observablesList[i]));
      canvas->Print(Form("result/plots/chap4_result_trending_%s.png", observablesList[i]));
      delete canvas;
      canvas = nullptr;
      // delete multigraph;
      multigraph = nullptr;
    }
    else { // 2p corr

      const int nGraphs = nObservables - 3; // 6
      canvas = new TCanvas("c", "", 1500, 900);
      canvas->SetBottomMargin(0.20);
      canvas->SetGrid();

      auto* multigraph = new TMultiGraph();
      multigraph->SetTitle(" ");
      const int markerStyle[nGraphs] = {20, 20, 20, 24, 24, 24};
      const int markerColor[nGraphs] = {kRed + 1, kBlue + 1, kGreen + 2, kRed + 1, kBlue + 1, kGreen + 2};
      for (int j = 0; j < nGraphs; j++) {

        auto* Graph = TrendingGraphs[i + j];
        setRunPositions(Graph, nRunNumberReal);
        // convertToRelativeDeviation(Graph);
        Graph->SetMarkerStyle(markerStyle[j]);
        Graph->SetMarkerColor(markerColor[j]);
        multigraph->Add(Graph);

      }

      const double shiftStep = 0.08;

      TIter next(multigraph->GetListOfGraphs());
      int graphIndex = 0;

      while (auto* object = next()) {
          auto* graph = dynamic_cast<TGraphErrors*>(object);
          if (!graph) {
              continue;
          }

          const double shift = (graphIndex - 0.5 * (nGraphs - 1)) * shiftStep;

          ShiftAlongXaxis(graph, shift);

          ++graphIndex;
      }

      const double yRangeMax = std::abs(getYMax(multigraph, false));
      const double yRangeMin = std::abs(getYMin(multigraph, false));
      multigraph->SetMinimum(yRangeMin - 0.2 * std::abs(yRangeMax));
      multigraph->SetMaximum(yRangeMax + 0.6 * std::abs(yRangeMax));


      multigraph->Draw("AP");//APX
      setRunLabels(multigraph, runNumberReal, nRunNumberReal);
      multigraph->GetXaxis()->SetTitle("Run number");
      multigraph->GetYaxis()->SetTitle("#LT#LT2#GT#GT_{n}");
      multigraph->GetXaxis()->SetTitleOffset(2);


      TLegend* legend = new TLegend(0.55, 0.70, 0.8, 0.85);
      legend->SetNColumns(2);
      legend->AddEntry((TObject*)0, "10\%-20\%", "");
      legend->AddEntry((TObject*)0, "50\%-60\%", "");
      legend->AddEntry(TrendingGraphs[i], "#LT#LT2#GT#GT_{2}");
      legend->AddEntry(TrendingGraphs[i + 3], "#LT#LT2#GT#GT_{2}");
      legend->AddEntry(TrendingGraphs[i + 1], "#LT#LT2#GT#GT_{3}");
      legend->AddEntry(TrendingGraphs[i + 4], "#LT#LT2#GT#GT_{3}");
      legend->AddEntry(TrendingGraphs[i + 2], "#LT#LT2#GT#GT_{4}");
      legend->AddEntry(TrendingGraphs[i + 5], "#LT#LT2#GT#GT_{4}");

      legend->SetBorderSize(0);
      legend->Draw();
      canvas->Print("result/plots/chap4_result_trending_2pcorr.pdf");
      canvas->Print("result/plots/chap4_result_trending_2pcorr.png");


      delete canvas;
      canvas = nullptr;
      // delete multigraph;
      multigraph = nullptr;

    }

  }

}

void PlotCutCheck(TFile* resultFile, const char* quantityName, const char* xTitle, const char* histListName, const char* centRange, const bool bLogY) {


  TCanvas* canvas = nullptr;
  TH1F* histBefore = nullptr;
  TH1F* histAfter = nullptr;

  // "multiparticle-cumulants"
  TString dirName = Form("multiparticle-cumulants%s", centRange);
  TDirectory* dir = dynamic_cast<TDirectory*>(resultFile->Get(dirName.Data()));
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

      // "___Histograms"
      TList* EVhistogramList = dynamic_cast<TList*>(baseList->FindObject(histListName));
      if(!EVhistogramList) {
          printf("  Failed to find list \"%s\"\n", histListName);
      }
      else {
        printf("  Open list \"%s\"\n", histListName);

        // "fHist___[eRec][before cut]"
        histBefore = dynamic_cast<TH1F*>(EVhistogramList->FindObject(Form("fHist%s[eRec][before cut]", quantityName)));
        if(!histBefore) {
            printf("  Failed to find histogram \"fHist%s[eRec][before cut]\"\n", quantityName);
        }
        else {
            printf("  Found histogram \"fHist%s[eRec][before cut]\"\n", quantityName);
        }

        // "fHist___[eRec][after cut]"
        histAfter = dynamic_cast<TH1F*>(EVhistogramList->FindObject(Form("fHist%s[eRec][after cut]", quantityName)));
        if(!histAfter) {
            printf("  Failed to find histogram \"fHist%s[eRec][after cut]\"\n", quantityName);
        }
        else {
            printf("  Found histogram \"fHist%s[eRec][after cut]\"\n", quantityName);
        }

        canvas = new TCanvas("c", "", 1500, 900);
        canvas->SetBottomMargin(0.20);
        canvas->SetGrid();
        canvas->SetLogy(bLogY);

        histBefore->GetXaxis()->SetTitle(xTitle);
        histBefore->GetYaxis()->SetTitle("");
        histBefore->SetTitle("");
        if (!bLogY) {
          histBefore->SetMinimum(0);
        } else {
          histBefore->SetMinimum(1e-1);
        }
        histBefore->SetLineColor(2);
        histAfter->SetLineColor(4);
        histBefore->SetLineWidth(2);
        histAfter->SetLineWidth(2);
        histBefore->Draw();
        histAfter->Draw("same");

        TLegend* legend = new TLegend(0.7, 0.68, 0.85, 0.78);
        legend->AddEntry(histBefore, "before cut");
        legend->AddEntry(histAfter, "after cut");
        legend->SetBorderSize(0);
        legend->Draw();

        canvas->Print(Form("result/plots/chap3_%sCutCheck%s.png", quantityName, centRange));
        canvas->Print(Form("result/plots/chap3_%sCutCheck%s.pdf", quantityName, centRange));
        delete canvas;
        canvas = nullptr;
        delete legend;
        legend = nullptr;
      }
    }
  }
}

void plotPhiWeights(TFile* weightFile, const int runNumber)
{

  TList* ccdbObject = nullptr;
  weightFile->GetObject("ccdb_object", ccdbObject);

  auto* run544116 = dynamic_cast<TList*>(ccdbObject->FindObject(Form("000%d", runNumber)));
  auto* h544116_0 = dynamic_cast<TH1F*>(run544116->FindObject(Form("wPhi_run%d", runNumber)));
  auto* h544116 = static_cast<TH1*>(h544116_0->Clone("hPhiWeight544116"));

  h544116->Scale(1.0 / h544116->Integral());

  h544116->SetDirectory(nullptr);
  // weightFile->Close();

  gStyle->SetOptStat(0);
  h544116->SetTitle("");
  h544116->GetXaxis()->SetTitle("#varphi");
  h544116->GetYaxis()->SetTitle("Weight");
  h544116->SetLineColor(kBlue + 1);
  h544116->SetLineWidth(2);
  const double yMin = h544116->GetMinimum();
  const double yMax = h544116->GetMaximum();
  const double margin = 0.08 * (yMax - yMin);
  h544116->SetMinimum(yMin - margin);
  h544116->SetMaximum(yMax + margin);

  auto* canvas = new TCanvas("cPhiWeights", "", 1000, 700);

  h544116->Draw("HIST");

  // auto* legend = new TLegend(0.72, 0.79, 0.89, 0.90);
  // legend->SetBorderSize(0);
  // legend->SetFillStyle(0);
  // legend->AddEntry(h544116, "run 544116", "l");
  // legend->Draw();
  canvas->Print(Form("result/plots/chap3_phiWeights_%d.png", runNumber));
}

void plotPtWeights(TFile* weightFile, const int runNumber)
{

  TList* ccdbObject = nullptr;
  weightFile->GetObject("ccdb_object", ccdbObject);

  auto* run544116 = dynamic_cast<TList*>(ccdbObject->FindObject(Form("000%d", runNumber)));
  auto* h544116_0 = dynamic_cast<TH1F*>(run544116->FindObject(Form("wPtReal_run%d", runNumber)));
  auto* h544116 = static_cast<TH1*>(h544116_0->Clone("hPtWeight544116"));

  h544116->Scale(1.0 / h544116->Integral());

  h544116->SetDirectory(nullptr);
  // weightFile->Close();

  gStyle->SetOptStat(0);
  h544116->SetTitle("");
  h544116->GetXaxis()->SetTitle("p_{T}");
  h544116->GetXaxis()->SetRangeUser(0.2, 5.);
  h544116->GetYaxis()->SetTitle("Weight");
  h544116->SetLineColor(kBlue + 1);
  h544116->SetLineWidth(2);
  const double yMin = h544116->GetMinimum();
  const double yMax = h544116->GetMaximum();
  const double margin = 0.08 * (yMax - yMin);
  h544116->SetMinimum(yMin - margin);
  h544116->SetMaximum(yMax + margin);

  auto* canvas = new TCanvas("cPtWeights", "", 1000, 700);

  h544116->Draw("HIST");

  // auto* legend = new TLegend(0.72, 0.79, 0.89, 0.90);
  // legend->SetBorderSize(0);
  // legend->SetFillStyle(0);
  // legend->AddEntry(h544116, "run 544116", "l");
  // legend->Draw();
  canvas->Print(Form("result/plots/chap3_ptWeights_%d.png", runNumber));
}


int Plotting(bool isNormalized = true, bool isCumulant = false)
{
  bNormalized = isNormalized;
  if (!isNormalized) {
    bShowMC = true;
  }
  bCumulant = isCumulant;

  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.045, "X");
  gStyle->SetTitleSize(0.05, "Y");
  gStyle->SetMarkerSize(1.5);

  Color_t semiTransWhite = TColor::GetColorTransparent(kWhite, 0.8);
  gStyle->SetLegendFillColor(semiTransWhite);
  gStyle->SetLegendBorderSize(1); 


  // data + stat err
  const auto resultData_vnscVScent_real = loadtable(resultPath_vnscVScent_real, true, false);
  const auto resultData_vnscVScent_MC = loadtable(resultPath_vnscVScent_MC, true, false);
  const auto dataV2Compare = loadtable("data/HEPData-1602.01119v2-csv/Table1.csv", true, false); // V2{2, |DELTAETA|>1}, V2{4}, V2{6}, V2{8}
  const auto dataV3V4Compare = loadtable("data/HEPData-1602.01119v2-csv/Table2.csv", true, false); // V3{2, |DELTAETA|>1}, V4{2, |DELTAETA|>1}
  const auto dataNSC23Compare = loadtable("data/HEPData-PhysRevLett.117.182301-csv/Table2.csv", true, false); // NSC(2,3)
  const auto dataNSC24Compare = loadtable("data/HEPData-PhysRevLett.117.182301-csv/Table2.csv", true, false); // NSC(2,4)
  const auto dataNSC234Compare = loadtable("data/HEPData-PhysRevLett.127.092302-csv/Table5.csv", true, false); // NSC(2,3,4)
  const auto dataSC23Compare = loadtable("data/HEPData-PhysRevLett.117.182301-csv/Table1.csv", true, false); // SC(2,3)
  const auto dataSC24Compare = loadtable("data/HEPData-PhysRevLett.117.182301-csv/Table1.csv", true, false); // SC(2,4)
  const auto dataSC234Compare = loadtable("data/HEPData-PhysRevLett.127.092302-csv/Table1.csv", true, false); // SC(2,3,4)

  // data + sys err
  const auto resultData_vnscVScent_real_sys = loadtable(resultPath_vnscVScent_real, false, true);
  const auto resultData_vnscVScent_MC_sys = loadtable(resultPath_vnscVScent_MC, false, true);
  const auto dataV2Compare_sys = loadtable("data/HEPData-1602.01119v2-csv/Table1.csv", false, true); // V2{2, |DELTAETA|>1}, V2{4}, V2{6}, V2{8}
  const auto dataV3V4Compare_sys = loadtable("data/HEPData-1602.01119v2-csv/Table2.csv", false, true); // V3{2, |DELTAETA|>1}, V4{2, |DELTAETA|>1}
  const auto dataNSC23Compare_sys = loadtable("data/HEPData-PhysRevLett.117.182301-csv/Table2.csv", false, true); // NSC(2,3)
  const auto dataNSC24Compare_sys = loadtable("data/HEPData-PhysRevLett.117.182301-csv/Table2.csv", false, true); // NSC(2,4)
  const auto dataNSC234Compare_sys = loadtable("data/HEPData-PhysRevLett.127.092302-csv/Table5.csv", false, true); // NSC(2,3,4)
  const auto dataSC23Compare_sys = loadtable("data/HEPData-PhysRevLett.117.182301-csv/Table1.csv", false, true); // SC(2,3)
  const auto dataSC24Compare_sys = loadtable("data/HEPData-PhysRevLett.117.182301-csv/Table1.csv", false, true); // SC(2,4)
  const auto dataSC234Compare_sys = loadtable("data/HEPData-PhysRevLett.127.092302-csv/Table1.csv", false, true); // SC(2,3,4)

  // others
  const auto resultData_trending = loadtable(resultPath_trending);
  const auto resultData_root_cutchecks = TFile::Open(resultPath_root_cutchecks, "read");
  const auto resultData_root_weights = TFile::Open(resultPath_weights, "read");


  if (!bCumulant) {

    // vn without comparison
    // plotFlow(2, dynamic_cast<TGraphErrors*>(resultData_vnscVScent_real[0]->Clone()), dynamic_cast<TGraphErrors*>(resultData_vnscVScent_MC[0]->Clone()), dynamic_cast<TGraphErrors*>(resultData_vnscVScent_real[1]->Clone()), dynamic_cast<TGraphErrors*>(resultData_vnscVScent_MC[1]->Clone()), dataV2Compare, dataV3V4Compare, false);
    // plotFlow(3, dynamic_cast<TGraphErrors*>(resultData_vnscVScent_real[2]->Clone()), dynamic_cast<TGraphErrors*>(resultData_vnscVScent_MC[2]->Clone()), dynamic_cast<TGraphErrors*>(resultData_vnscVScent_real[3]->Clone()), dynamic_cast<TGraphErrors*>(resultData_vnscVScent_MC[3]->Clone()), dataV2Compare, dataV3V4Compare, false);
    // plotFlow(4, dynamic_cast<TGraphErrors*>(resultData_vnscVScent_real[4]->Clone()), dynamic_cast<TGraphErrors*>(resultData_vnscVScent_MC[4]->Clone()), dynamic_cast<TGraphErrors*>(resultData_vnscVScent_real[5]->Clone()), dynamic_cast<TGraphErrors*>(resultData_vnscVScent_MC[5]->Clone()), dataV2Compare, dataV3V4Compare, false);

    // vn with comparison
    plotFlow(2, dynamic_cast<TGraphErrors*>(resultData_vnscVScent_real[11]->Clone()), dynamic_cast<TGraphErrors*>(resultData_vnscVScent_MC[11]->Clone()), resultData_vnscVScent_real[1], resultData_vnscVScent_MC[1], dataV2Compare, dataV3V4Compare, 
                dynamic_cast<TGraphErrors*>(resultData_vnscVScent_real_sys[11]->Clone()), dynamic_cast<TGraphErrors*>(resultData_vnscVScent_MC_sys[11]->Clone()), resultData_vnscVScent_real_sys[1], resultData_vnscVScent_MC_sys[1], dataV2Compare_sys, dataV3V4Compare_sys, true);
    plotFlow(3, dynamic_cast<TGraphErrors*>(resultData_vnscVScent_real[14]->Clone()), dynamic_cast<TGraphErrors*>(resultData_vnscVScent_MC[14]->Clone()), resultData_vnscVScent_real[3], resultData_vnscVScent_MC[3], dataV2Compare, dataV3V4Compare, 
                dynamic_cast<TGraphErrors*>(resultData_vnscVScent_real_sys[14]->Clone()), dynamic_cast<TGraphErrors*>(resultData_vnscVScent_MC_sys[14]->Clone()), resultData_vnscVScent_real_sys[3], resultData_vnscVScent_MC_sys[3], dataV2Compare_sys, dataV3V4Compare_sys, true);
    plotFlow(4, dynamic_cast<TGraphErrors*>(resultData_vnscVScent_real[17]->Clone()), dynamic_cast<TGraphErrors*>(resultData_vnscVScent_MC[17]->Clone()), resultData_vnscVScent_real[5], resultData_vnscVScent_MC[5], dataV2Compare, dataV3V4Compare, 
                dynamic_cast<TGraphErrors*>(resultData_vnscVScent_real_sys[17]->Clone()), dynamic_cast<TGraphErrors*>(resultData_vnscVScent_MC_sys[17]->Clone()), resultData_vnscVScent_real_sys[5], resultData_vnscVScent_MC_sys[5], dataV2Compare_sys, dataV3V4Compare_sys, true);
    
    plotFlow2pGap(2, resultData_vnscVScent_real[0], resultData_vnscVScent_real[9], resultData_vnscVScent_real[10], dynamic_cast<TGraphErrors*>(resultData_vnscVScent_real[11]->Clone()));
    plotFlow2pGap(3, resultData_vnscVScent_real[2], resultData_vnscVScent_real[12], resultData_vnscVScent_real[13], dynamic_cast<TGraphErrors*>(resultData_vnscVScent_real[14]->Clone()));
    plotFlow2pGap(4, resultData_vnscVScent_real[4], resultData_vnscVScent_real[15], resultData_vnscVScent_real[16], dynamic_cast<TGraphErrors*>(resultData_vnscVScent_real[17]->Clone()));

    if (bNormalized) {
      plotNSC(3, resultData_vnscVScent_real[6], resultData_vnscVScent_MC[6], dataNSC23Compare[0], true, bShowMC, bNormalized);
      plotNSC(4, resultData_vnscVScent_real[7], resultData_vnscVScent_MC[7], dataNSC24Compare[1], true, bShowMC, bNormalized);
      plotNSC(5, resultData_vnscVScent_real[8], resultData_vnscVScent_MC[8], dataNSC234Compare[0], true, bShowMC, bNormalized);
    }
    else {
      plotNSC(3, resultData_vnscVScent_real[6], resultData_vnscVScent_MC[6], dataSC23Compare[0], true, bShowMC, bNormalized);
      plotNSC(4, resultData_vnscVScent_real[7], resultData_vnscVScent_MC[7], dataSC24Compare[1], true, bShowMC, bNormalized);
      plotNSC(5, resultData_vnscVScent_real[8], resultData_vnscVScent_MC[8], dataSC234Compare[0], true, bShowMC, bNormalized);
    }

    if (!bLocal) {
      plotTrending(resultData_trending);
      PlotCutCheck(resultData_root_cutchecks, "Multiplicity", "Multiplicity", "EventHistograms", "_cent60-70", true);
      PlotCutCheck(resultData_root_cutchecks, "Multiplicity", "Multiplicity", "EventHistograms", "_cent50-60", true);
      PlotCutCheck(resultData_root_cutchecks, "Multiplicity", "Multiplicity", "EventHistograms", "_cent40-50", true);
      PlotCutCheck(resultData_root_cutchecks, "Multiplicity", "Multiplicity", "EventHistograms", "_cent30-40", true);
      PlotCutCheck(resultData_root_cutchecks, "Multiplicity", "Multiplicity", "EventHistograms", "_cent20-30", true);
      PlotCutCheck(resultData_root_cutchecks, "Multiplicity", "Multiplicity", "EventHistograms", "_cent10-20", true);
      PlotCutCheck(resultData_root_cutchecks, "Multiplicity", "Multiplicity", "EventHistograms", "_cent05-10", true);
      PlotCutCheck(resultData_root_cutchecks, "Multiplicity", "Multiplicity", "EventHistograms", "_cent00-05", true);
      PlotCutCheck(resultData_root_cutchecks, "Pt", "p_{T} [GeV]", "ParticleHistograms", "", true);
      PlotCutCheck(resultData_root_cutchecks, "Eta", "#eta", "ParticleHistograms", "", false);
      PlotCutCheck(resultData_root_cutchecks, "Phi", "#phi", "ParticleHistograms", "", false);
      PlotCutCheck(resultData_root_cutchecks, "DCAXY", "DCA_{XY} [cm]", "ParticleHistograms", "", false);
      PlotCutCheck(resultData_root_cutchecks, "DCAZ", "DCA_{Z} [cm]", "ParticleHistograms", "", false);
      PlotCutCheck(resultData_root_cutchecks, "VertexZ", "Vertex Z [cm]", "EventHistograms", "", false);

      for (int const& i : runNumberMC) {
        plotPtWeights(resultData_root_weights, i);
      }
      for (int const& i : runNumberReal) {
        plotPhiWeights(resultData_root_weights, i);
      }

    }

  } else {
    return 1;
    // PlotCumulant()
  }

  return 0;
}
