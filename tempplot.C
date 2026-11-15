void tempplot()
{

TFile* weightFile = TFile::Open("weights.root", "read");
TFile* weightFile_ = TFile::Open("weights_.root", "read");
const int runNumber = 544122;

  TList* ccdbObject = nullptr;
  weightFile->GetObject("ccdb_object", ccdbObject);

  auto* run544122 = dynamic_cast<TList*>(ccdbObject->FindObject(Form("000%d", runNumber)));
  auto* h544122_0 = dynamic_cast<TH1F*>(run544122->FindObject(Form("wPhi_run%d", runNumber)));
  auto* h544122 = static_cast<TH1*>(h544122_0->Clone("hPhiWeight544122"));

  h544122->SetDirectory(nullptr);
  // weightFile->Close();


  TList* ccdbObject_ = nullptr;
  weightFile_->GetObject("ccdb_object", ccdbObject_);

  auto* run544122_ = dynamic_cast<TList*>(ccdbObject_->FindObject(Form("000%d", runNumber)));
  auto* h544122_0_ = dynamic_cast<TH1F*>(run544122_->FindObject(Form("wPhi_run%d", runNumber)));
  auto* h544122_ = static_cast<TH1*>(h544122_0_->Clone("hPhiWeight544122_"));
  h544122_->SetDirectory(nullptr);

  h544122->Divide(h544122_);

  gStyle->SetOptStat(0);
  h544122->SetTitle("");
  h544122->GetXaxis()->SetTitle("#varphi");
  h544122->GetYaxis()->SetTitle("Weight");
  h544122->SetLineColor(kBlue + 1);
  h544122->SetLineWidth(2);
  const double yMin = h544122->GetMinimum();
  const double yMax = h544122->GetMaximum();
  const double margin = 0.08 * (yMax - yMin);
  h544122->SetMinimum(yMin - margin);
  h544122->SetMaximum(yMax + margin);

  auto* canvas = new TCanvas("cPhiWeights", "", 1000, 700);

  h544122->Draw("HIST");

  // auto* legend = new TLegend(0.72, 0.79, 0.89, 0.90);
  // legend->SetBorderSize(0);
  // legend->SetFillStyle(0);
  // legend->AddEntry(h544122, "run 544122", "l");
  // legend->Draw();
  canvas->Print(Form("temp_checkphiweight_%d.png", runNumber));
}

