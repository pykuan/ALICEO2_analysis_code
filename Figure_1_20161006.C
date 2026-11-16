// Plotting macro for Fig. 1 in "Correlated event-by-event fluctuations of flow harmonics in Pb-Pb collisions at $\sqrt{s_{_{\rm NN}}}=2.76$~TeV" 
// Version: 20161006 (Remark: If you upgrade the version, replace ALL occurences in the macro)

// Save the figures in pdf, eps and png formats (in the directory from which the macro was launched):
Bool_t bSaveFigures = kTRUE; 

// Functions:
void SetStyle(Bool_t graypalette=kFALSE);
void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15);
void DrawLogo (Int_t logo=0, Double_t xmin =  0.28, Double_t ymin= 0.68);
void FakeHistosOnlyForExample(TH1*&hstat, TH1*&hsyst, TH1*&hsystCorr);
void LoadLibs();

// Main:
void Figure_1_20161006() 
{
 // Load necessary libraries
 LoadLibs();
 // Set the default style
 SetStyle();

 // Prepare Figure, please stick to the default canvas size(s) unless absolutely necessary in your case
 // a) Rectangular:
 TCanvas *cfig = new TCanvas("cfig", "Alice Figure Template", 800, 600); 
 // b) Square:
 // TCanvas *cfig = new TCanvas("cfig", "Alice Figure Template", 800, 800); 
 // cfig->SetLogy();

 // Set titles etc..
 TH1 *h = cfig->DrawFrame(0,-2.69e-6,70,3.59e-6);
 // Choose titles:
 const char *xTitle = "Centrality percentile";
 const char *yTitle = "SC(m,n)";
 // Set titles:
 h->SetXTitle(xTitle);
 h->SetYTitle(yTitle);
 // Set title offsets:
 h->GetXaxis()->SetTitleOffset(0.98844);
 h->GetYaxis()->SetTitleOffset(0.744);

 //===================================================================================================

 // SC(3,2) in Hijing:
 Double_t x_SC32_Hijing[9] = {2.5, 7.5, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0};
 Double_t y_SC32_Hijing[9] = {-1.139098e-09, -2.878181e-09, -5.721163e-09, -5.849574e-09, 9.227799e-09, -4.069350e-08, -8.059965e-08, 1.194018e-07, 1.223952e-06};
 Double_t xErr_SC32_Hijing[9] = {0.};
 Double_t yErr_SC32_Hijing[9] = {2.141186e-09, 2.887575e-09, 3.528635e-09, 7.114722e-09, 1.565906e-08, 3.874696e-08, 1.140743e-07, 4.394331e-07, 2.206863e-06};
 Int_t nPoints_SC32_Hijing = sizeof(x_SC32_Hijing)/sizeof(Double_t);
 TGraphErrors *ge_SC32_Hijing = new TGraphErrors(nPoints_SC32_Hijing-1,x_SC32_Hijing,y_SC32_Hijing,xErr_SC32_Hijing,yErr_SC32_Hijing);
 ge_SC32_Hijing->SetMarkerStyle(kOpenCircle);
 ge_SC32_Hijing->SetMarkerColor(kBlue+1);
 ge_SC32_Hijing->SetMarkerSize(1.2);
 ge_SC32_Hijing->SetLineColor(kBlue+1);
 ge_SC32_Hijing->SetFillColor(kBlue-8);
 ge_SC32_Hijing->SetFillStyle(3012);
 ShiftAlongXaxis(ge_SC32_Hijing,-0.5);
    
 //===================================================================================================

 // SC(4,2) in Hijing:
 Double_t x_SC42_Hijing[9] = {2.5, 7.5, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0};
 Double_t y_SC42_Hijing[9] = {-1.006451e-09, -2.088190e-09, -2.288267e-09, -1.243485e-08, -1.476088e-08, 2.919063e-10, -2.135256e-07, 1.557772e-07, -5.178595e-07};
 Double_t xErr_SC42_Hijing[9] = {0.};
 Double_t yErr_SC42_Hijing[9] = {1.815147e-09, 2.527433e-09, 3.180418e-09, 6.386532e-09, 1.441397e-08, 3.641242e-08, 1.068634e-07, 4.213652e-07, 2.079544e-06};
 Int_t nPoints_SC42_Hijing = sizeof(x_SC42_Hijing)/sizeof(Double_t);
 TGraphErrors *ge_SC42_Hijing = new TGraphErrors(nPoints_SC42_Hijing-1,x_SC42_Hijing,y_SC42_Hijing,xErr_SC42_Hijing,yErr_SC42_Hijing);
 ge_SC42_Hijing->SetMarkerStyle(kOpenSquare);
 ge_SC42_Hijing->SetMarkerColor(kRed);
 ge_SC42_Hijing->SetMarkerSize(1.);
 ge_SC42_Hijing->SetLineColor(kRed);
 ge_SC42_Hijing->SetFillColor(kRed-10);
 ge_SC42_Hijing->SetFillStyle(3001);
 ge_SC42_Hijing->SetFillColor(kRed-8);
 ge_SC42_Hijing->SetFillStyle(3008);
 ShiftAlongXaxis(ge_SC42_Hijing,0.5);
    
 //===================================================================================================

 // SC(3,2) values and statistical erros in data:
 Double_t x_SC32[9] = {2.5, 7.5, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0};
 Double_t y_SC32[9] = {-1.474034e-08, -6.189796e-08, -3.210893e-07, -6.376161e-07, -1.026290e-06, -1.254218e-06, -1.077982e-06, -1.093781e-06, -4.616095e-07};
 Double_t xErr_SC32[9] = {0.};
 Double_t yErr_SC32[9] = {1.856762e-09, 4.132604e-09, 7.362652e-09, 1.535372e-08, 2.897027e-08, 5.509401e-08, 1.154303e-07, 3.056604e-07, 1.284629e-06};
 Int_t nPoints_SC32 = sizeof(x_SC32)/sizeof(Double_t);
 TGraphErrors *ge_SC32 = new TGraphErrors(nPoints_SC32,x_SC32,y_SC32,xErr_SC32,yErr_SC32);
 ge_SC32->SetMarkerStyle(kFullCircle);
 ge_SC32->SetMarkerColor(kBlue+1);
 ge_SC32->SetMarkerSize(1.2);
 ge_SC32->SetLineColor(kBlue+1);
 //ShiftAlongXaxis(ge_SC32,0.0);
    
 // SC(3,2) systematical error:
 Double_t sysError_SC32 = 0.098;
 Double_t x_SysError_SC32[9] = {2.5, 7.5, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0};
 Double_t y_SysError_SC32[9] = {-1.474034e-08, -6.189796e-08, -3.210893e-07, -6.376161e-07, -1.026290e-06, -1.254218e-06, -1.077982e-06, -1.093781e-06, -4.616095e-07};
 Double_t xErr_SysError_SC32[9] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};
 Double_t yErr_SysError_SC32[9] = {0.0};
 for (Int_t p=0; p<9; p++)
 {
  yErr_SysError_SC32[p] = sysError_SC32 * y_SC32[p];
 }
 Int_t nPoints_SysError_SC32 = sizeof(x_SysError_SC32)/sizeof(Double_t);
 TGraphErrors *ge_SysError_SC32 = new TGraphErrors(nPoints_SysError_SC32,x_SysError_SC32,y_SysError_SC32,xErr_SysError_SC32,yErr_SysError_SC32);
 ge_SysError_SC32->SetMarkerStyle(kFullCircle);
 ge_SysError_SC32->SetMarkerColor(kBlue+1);
 ge_SysError_SC32->SetMarkerSize(1.2);
 ge_SysError_SC32->SetLineColor(kBlue+1);
 ge_SysError_SC32->SetFillColor(kBlue-9);
 //ShiftAlongXaxis(ge_SysError_SC32,0.0);

 //===================================================================================================

 // SC(4,2) values and statistical erros in data:
 Double_t x_SC42[9] = {2.5, 7.5, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0};
 Double_t y_SC42[9] = {1.382924e-09, 1.739164e-08, 8.423141e-08, 4.483835e-07, 1.126789e-06, 1.861415e-06, 2.147265e-06, 2.207977e-06, 3.535638e-06};
 Double_t xErr_SC42[9] = {0.};
 Double_t yErr_SC42[9] = {1.195944e-09, 2.724491e-09, 5.124662e-09, 1.194504e-08, 2.510275e-08, 5.113180e-08, 1.106269e-07, 3.004486e-07, 1.253560e-06};
 Int_t nPoints_SC42 = sizeof(x_SC42)/sizeof(Double_t);
 TGraphErrors *ge_SC42 = new TGraphErrors(nPoints_SC42,x_SC42,y_SC42,xErr_SC42,yErr_SC42);
 ge_SC42->SetMarkerStyle(kFullSquare);
 ge_SC42->SetMarkerColor(kRed);
 ge_SC42->SetMarkerSize(1.);
 ge_SC42->SetLineColor(kRed);
 //ShiftAlongXaxis(ge_SC42,0.0);
    
 // SC(4,2) systematical error:
 Double_t sysError_SC42 = 0.169;
 Double_t x_SysError_SC42[9] = {2.5, 7.5, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0};
 Double_t y_SysError_SC42[9] = {1.382924e-09, 1.739164e-08, 8.423141e-08, 4.483835e-07, 1.126789e-06, 1.861415e-06, 2.147265e-06, 2.207977e-06, 3.535638e-06};
 Double_t xErr_SysError_SC42[9] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};
 Double_t yErr_SysError_SC42[9] = {0.0};
 for (Int_t p=0; p<9; p++)
 {
  yErr_SysError_SC42[p] = sysError_SC42 * y_SysError_SC42[p];
 }
 Int_t nPoints_SysError_SC42 = sizeof(x_SysError_SC42)/sizeof(Double_t);
 TGraphErrors *ge_SysError_SC42 = new TGraphErrors(nPoints_SysError_SC42,x_SysError_SC42,y_SysError_SC42,xErr_SysError_SC42,yErr_SysError_SC42);
 ge_SysError_SC42->SetMarkerStyle(kFullSquare);
 ge_SysError_SC42->SetMarkerColor(kRed);
 ge_SysError_SC42->SetMarkerSize(1.);
 ge_SysError_SC42->SetLineColor(kRed);
 ge_SysError_SC42->SetFillColor(kRed-10);
 //ShiftAlongXaxis(ge_SysError_SC42,0.0);
    
 //===================================================================================================    

 // Draw the logo:   
 //  0: Just "ALICE" (for final data), to be added only if ALICE does not appear otherwise (e.g. in the legend)
 //  >0: ALICE Preliminary
 // DrawLogo(1, 0.59, 0.81);

 // You should always specify the colliding system:
 // NOTATION: pp, p-Pb, Pb-Pb. 
 // Don't forget to use #sqrt{s_{NN}} for p-Pb and Pb-Pb
 // You can change the position of this with
 /* 
 TLatex *text = new TLatex (45, 2.89e-6,"Pb-Pb #sqrt{#it{s}}_{NN} = 2.76 TeV");
 text->SetTextSizePixels(24);
 text->Draw();
 */
    
 // Draw line at y=0:
 TLine *line0 = new TLine(0., 0., 70., 0.);
 line0->SetLineColor(1);
 line0->SetLineWidth(2);
 line0->SetLineStyle(1);
 line0->Draw("lsame");
    
 // Legend for ALICE:
 TLegend *legend_ALICE = new TLegend(0.19,0.61,0.465,0.87);
 legend_ALICE->SetFillStyle(0); // white legend background
 legend_ALICE->SetTextSize(0.035);
 //legend_ALICE->SetBorderSize(0.0);
 legend_ALICE->SetTextFont(42);
 legend_ALICE->SetMargin(0.3);
 legend_ALICE->SetHeader("ALICE Pb-Pb #sqrt{#it{s}}_{NN} = 2.76 TeV");
 legend_ALICE->AddEntry(ge_SysError_SC42,"SC(4,2)","pf");
 legend_ALICE->AddEntry(ge_SysError_SC32,"SC(3,2)","pf");
  
 // Legend for Hijing:
 TLegend *legend_Hijing = new TLegend(0.18,0.19,0.49,0.38);
 legend_Hijing->SetFillStyle(0); // white legend background
 legend_Hijing->SetTextSize(0.035);
 //legend_Hijing->SetBorderSize(0.0);
 //legend_Hijing->SetTextFont(22);
 legend_Hijing->SetMargin(0.3);
 legend_Hijing->SetHeader("HIJING");
 legend_Hijing->AddEntry(ge_SC42_Hijing,"SC(4,2)","p");
 legend_Hijing->AddEntry(ge_SC32_Hijing,"SC(3,2)","p");
    
 // Final drawing (ordering matters!):   
 legend_ALICE->Draw("same");
 legend_Hijing->Draw("same");
 ge_SC32_Hijing->Draw("psame,Z");
 ge_SC42_Hijing->Draw("psame,Z");
 ge_SysError_SC32->Draw("2same");
 ge_SC32->Draw("psame,Z");
 ge_SysError_SC42->Draw("2same");
 ge_SC42->Draw("psame,Z");
    
 // Save or die:   
 if(bSaveFigures)
 {
  cfig->SaveAs("Figure_1_20161006.pdf");
  cfig->SaveAs("Figure_1_20161006.eps");
  cfig->SaveAs("Figure_1_20161006.png");
 } // if(bSaveFigures)

} // Figure_1_20161006()

//===================================================================================================    

void ShiftAlongXaxis(TGraphErrors *ge, Double_t shift)
{
 // Shift original TGraphErrors along x-axis by Double_t shift.
 
 if(!ge){cout<<"!ge"<<endl;exit(0);}
 
 Int_t nPoints = ge->GetN();
 Double_t x = 0.;
 Double_t y = 0.;
 for(Int_t p=0;p<nPoints;p++)
 { 
  ge->GetPoint(p,x,y);
  x+=shift;
  ge->SetPoint(p,x,y);
 } // end of for(Int_t p=0;p<nPoints;p++)

} // end of void ShiftAlongXaxis(TGraphErrors *ge, Double_t shift)

//===================================================================================================    

// Functions from ALICE Figures Template macro:

void LoadLibs() {
  gSystem->Load("libCore.so");  
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
}

void DrawLogo (Int_t logo, Double_t xmin, Double_t ymin) {

  // Logo is not needed anymore, now we only write alice preliminary
  // Logo:
  // 0: Justr writes "ALICE" (for final data)
  // Anything eles: writes "ALICE Preliminary"

  TLatex *   tex = new TLatex(xmin,ymin, logo ? "ALICE Preliminary" : "ALICE");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->Draw();

  // OLD logo
  //  TPad * currentPad = gPad;
  // Double_t AliLogo_LowX =xmin;
  // Double_t AliLogo_LowY = ymin;
  // Double_t AliLogo_Height = size;
  // //ALICE logo is a  file that is 821x798 pixels->should be wider than a square
  // Double_t AliLogo_Width  = (821./798.) * AliLogo_Height * gPad->GetWh() / gPad->GetWw();
  
  // TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",AliLogo_LowX,AliLogo_LowY,AliLogo_LowX+AliLogo_Width,AliLogo_LowY+AliLogo_Height);
  // myPadSetUp(myPadLogo,0,0,0,0);
  // //myPadLogo->SetFixedAspectRatio(1);
  // myPadLogo->Draw();
  // myPadLogo->cd();
  // if (logo == 0) {
  //   myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
  // } else if (logo == 1){
  //   TASImage *myAliceLogo = new TASImage(performanceLogoPath);
  //   myAliceLogo->Draw();
  // } else if (logo == 2) {
  //   TASImage *myAliceLogo = new TASImage(preliminaryLogoPath);
  //   myAliceLogo->Draw();
  // }
  // // go back to the old pad
  // currentPad->cd();

}

void myPadSetUp(TPad *currentPad, float currentLeft, float currentTop, float currentRight, float currentBottom){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}

void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;
  
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);

}

void FakeHistosOnlyForExample(TH1* &hstat, TH1* &hsyst, TH1*&hsystCorr) {

  TF1 * fExpo = new TF1 ("fExpo", "expo");
  fExpo->SetParameters(10, -0.3);
  hstat     = new TH1F("hstat", "hstat", 100, 0, 10);
  hsyst     = new TH1F("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1F("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2(); 
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2(); 
  hstat->FillRandom("fExpo",20000);
  Int_t nbinx = hstat->GetNbinsX();
  
  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
			  hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
			  hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
			  hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }
  
}
