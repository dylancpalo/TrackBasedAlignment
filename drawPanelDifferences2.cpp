#include <math.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include "TSystem.h"
#include "TProfile2D.h"

Double_t doubleGaussianConv(Double_t *x, Double_t *par){
   Float_t xx = x[0];

   Double_t sig1Sq = par[2]*par[2];
   Double_t sig2Sq = par[5]*par[5];
   Double_t mu1 = par[1];
   Double_t mu2 = par[1];
   Double_t gaus1 = par[0]*par[0]*exp(-0.5*((xx-mu1*2)*(xx-mu1*2)/(2*sig1Sq)));
   Double_t gaus2 = par[0]*par[0]*par[3]*par[3]*exp(-0.5*((xx-mu2*2)*(xx-mu2*2)/(2*sig2Sq)));
   Double_t gaus3 = 2*par[0]*par[0]*par[3]*exp(-0.5*((xx-(mu1+mu2))*(xx-(mu1+mu2))/(sig1Sq+sig2Sq)));

   return gaus1+gaus2+gaus3;
}

void drawPanelDifferences2(const char *badWireFileName = "CombinedTrkStrawStatus_08_21.txt"){

   int drawWireHistograms = 0;  // Set to 1 to draw the ~20k individual-wire histograms and fits.

   gStyle->SetOptStat(0);
   gStyle->SetStatFontSize(0.8);
   gStyle->SetTitleFontSize(0.07);
   gStyle->SetPadBottomMargin(0.13);
   gStyle->SetPadRightMargin(0.05);
   gStyle->SetPadLeftMargin(0.1);
   char name[600];
   char name2[600];

   snprintf(name, sizeof(name), "09_17_DT_DX_Planev7_ZMinus");
   snprintf(name2, sizeof(name2), "6Panel");

   char fileName[600];
   snprintf(fileName, sizeof(fileName), "TrackerROOT6Panel_09_21_DX_DT_Planev7_ZMinus.root");
   cout << fileName << endl;
   // Keep all CSV output for this configuration/run together.
   char csvFolder[600];
   snprintf(csvFolder, sizeof(csvFolder), "CSV_%s_%s", name, name2);
   gSystem->mkdir(csvFolder, true);

   char csvName[600];
   snprintf(csvName, sizeof(csvName), "%s/TrackShiftsWire_%s_%s.csv", csvFolder, name, name2);

   TFile *f1 = TFile::Open(fileName);

   TCanvas *cZPOCAErrorY = new TCanvas("cZPOCAErrorY", "cZPOCAErrorY", 0, 0, 1350, 1500);
   cZPOCAErrorY->Divide(1);
   cZPOCAErrorY->cd(1);
   TProfile *hZPOCAErrorY; f1->GetObject("a1/_hZPOCAErrorY", hZPOCAErrorY);
   hZPOCAErrorY->GetYaxis()->SetRangeUser(-20, 20);
   hZPOCAErrorY->Draw();
   cZPOCAErrorY->SaveAs("cZPOCAErrorY.pdf");

   TF1 *fDOCA = new TF1("fDOCA","[0]*exp(-0.5*(x-[1])*(x-[1])/([2]*[2])) + [0]*[5]*exp(-0.5*(x-[3])*(x-[3])/([4]*[4]))",-750, 750);
   fDOCA->SetParameter(0,50000);
   fDOCA->SetParameter(1,0.);
   fDOCA->SetParameter(2,200);
   fDOCA->SetParameter(3,-150.);
   fDOCA->SetParameter(4,400.);
   fDOCA->SetParameter(5,0.2);
   fDOCA->SetParLimits(5,0.25, 0.25);
   fDOCA->SetParLimits(2,150, 500);
   fDOCA->SetParLimits(4,300, 1300);
   fDOCA->SetParLimits(1,-250, 250);
   fDOCA->SetParLimits(3,-500, 500);

   TCanvas *cDOCAResolution = new TCanvas("cDOCAResolution", "cDOCAResolution", 0, 0, 1350, 1500);
   char title[600];

   const Int_t nStations = 18;
   const Int_t nPanels = 12;
   const Int_t nWiresPerSector = 12;
   const Int_t nSectorsPerPanel = 8;
   const Int_t nWiresPerPanel = nWiresPerSector * nSectorsPerPanel;
   const Int_t nPanelHistograms = nStations * nPanels;
   const Int_t nWireHistograms = nStations * nPanels * nWiresPerPanel;

   // MN-ID lookup from TrkPanelMap.toml for run range [121700, 200000].
   // Indexing is [frame/slot][panel-in-station].
   // panel-in-station 0-5  -> even plane 2*slot, TOML panel 0-5
   // panel-in-station 6-11 -> odd  plane 2*slot+1, TOML panel 0-5
   const Int_t panelMNID[nStations][nPanels] = {
      {120, 244, 126, 210, 269, 191,  73, 242,  74,  65,  77,  38},
      { 56, 239,  26,  99,  59,  60, 127,  90, 130, 119, 108,  88},
      { 42,  83,  70, 209, 279, 256,  30, 215, 185,  81, 189, 115},
      { 40, 136, 168,  96, 134,  94, 135, 131, 176, 181,  79, 179},
      {100, 196, 159, 152, 245, 216, 192, 184, 227, 194, 183, 195},
      { 93, 118,  34,  85, 111, 113, 243, 175, 257, 260, 259, 252},
      {170, 187, 172, 186, 171, 211, 221, 188, 282, 277, 177, 205},
      {161, 166, 148, 163, 154, 167, 246,  46, 217,  36, 238,  44},
      {213, 235, 247, 101, 253, 219, 224, 273, 276, 261, 248, 262},
      { 45, 106, 137,  78, 129,  98, 228, 204, 230, 208, 223, 226},
      {150, 199, 278,  80, 122, 133, 202, 200, 270, 112,  72, 193},
      {169,  39,  43,  53,  52,  35, 229, 275, 274, 254, 151, 225},
      {207, 139, 155, 145, 132,  82,  64,  63,  67,  55,  69,  62},
      {233,  66, 117, 138, 263, 264, 158, 250, 272, 280, 153, 220},
      {103, 147,  84, 258, 203, 104,  57,  47,  48, 128,  51, 105},
      {251, 231,  68, 142, 109, 182, 174, 236,  41, 206, 164, 214},
      {281, 240, 265, 271, 190, 255, 241,  86, 178,  61,  54, 266},
      { 97, 123,  92, 121,  89, 125, 160,  91,  71, 107, 149, 124}
   };

   // Bad-wire masks used only for the panel-level occupancy overlays.
   bool badWire[nStations][nPanels][nWiresPerPanel] = {};
   bool badWireShort[nStations][nPanels][nWiresPerPanel] = {};

   {
      std::ifstream badWireFile(badWireFileName);
      if (!badWireFile.is_open()){
         std::cerr << "WARNING: Could not open bad-wire file: " << badWireFileName << std::endl;
      } else {
         std::string line;
         Int_t nBadLoaded = 0;

         while (std::getline(badWireFile, line)){
            if (line.empty() || line[0] == '#') continue;

            std::istringstream iss(line);
            Int_t station = -1;
            Int_t panel = -1;
            Int_t wire = -1;
            Int_t bad = 0;
            if (!(iss >> station >> panel >> wire >> bad)) continue;

            if (station < 0 || station >= nStations) continue;
            if (panel < 0 || panel >= nPanels) continue;
            if (wire < 0 || wire >= nWiresPerPanel) continue;

            if (bad == 0){
               badWire[station][panel][wire] = true;
            } else if (bad == 1){
               badWireShort[station][panel][wire] = true;
            }
            nBadLoaded++;
         }

         std::cout << "Loaded " << nBadLoaded
                   << " bad-wire entries from " << badWireFileName << std::endl;
      }
   }

   cDOCAResolution->Divide(1);
   TH1F *h1Y; f1->GetObject("a1/_hDOCAErrorCleaned", h1Y);
   h1Y->SetLineColor(kBlack);
   cDOCAResolution->cd(1);
   h1Y->GetXaxis()->SetTitle("RDrift - UDOCA [um]");
   h1Y->SetTitle("RDrift - UDOCA [um]");
   if (h1Y->GetEntries() >= 500){
      h1Y->Fit(fDOCA, "", "");
   } else {
      std::cout << "Skipping global DOCA double-Gaussian fit: only "
                << h1Y->GetEntries() << " entries." << std::endl;
   }
   h1Y->Draw();
         h1Y->GetXaxis()->SetRangeUser(-1000, 1000);

   Double_t meanCore = fDOCA->GetParameter(1);
   Double_t sigmaCore = fDOCA->GetParameter(2);
   Double_t fractionTail = fDOCA->GetParameter(5);
   Double_t meanTail = fDOCA->GetParameter(3);
   Double_t sigmaTail = fDOCA->GetParameter(4);
   auto legendY = new TLegend(0.4,0.64, 0.95,0.9);
   legendY->SetTextSize(0.035);
   snprintf(title, sizeof(title), "#splitline{#sigma_{1}=%.2f, #mu_{1}=%.1f [um]}{#splitline{#sigma_{2}=%.2f, #mu_{2}=%.1f [um]}{A_{2}/A_{1}=%.2f}}", sigmaCore,  meanCore, sigmaTail, meanTail, fractionTail);
   legendY->AddEntry(h1Y,title,"l");
   legendY->SetFillStyle(0);
   legendY->Draw();

   int bin0 = h1Y->FindFirstBinAbove(h1Y->GetMaximum()/2);
   int bin1 = h1Y->FindLastBinAbove(h1Y->GetMaximum()/2);
   float FWHM = h1Y->GetBinCenter(bin1) - h1Y->GetBinCenter(bin0);
   std::cout << "FWHM: " << FWHM << std::endl;
   snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/DOCA_resolution_%s_%s.pdf", name, name2);

   cDOCAResolution->SaveAs(title);
   fDOCA->SetParameter(0,5000);
   
   TF1 *ftanh = new TF1("ftanh","[0]*tanh([3]*x-[1]) + [2]",-12, 12);
   ftanh->SetParameter(0,500);
   ftanh->SetParameter(1,4.);
   ftanh->SetParameter(2,500);
   ftanh->SetParameter(3,.5);

   TH1F *hdriftTime[nPanelHistograms];

   for (Int_t istation = 0; istation < nStations; istation++){
      snprintf(title, sizeof(title), "cdriftTimePanel_Station%02i", istation);
      TCanvas *cdriftTimePanel = new TCanvas(title, title, 0, 0, 1350, 1500);
      cdriftTimePanel->Divide(3, 4);

      for (Int_t ihist = 0; ihist < nPanels; ihist++){
         const Int_t globalHist = istation * nPanels + ihist;
         cdriftTimePanel->cd(ihist+1);
         snprintf(title, sizeof(title), "a1/hdriftTimePanel%i", globalHist);
         f1->GetObject(title, hdriftTime[globalHist]);
         if (!hdriftTime[globalHist]) continue;
         hdriftTime[globalHist]->SetLineColor(kBlack);
         hdriftTime[globalHist]->GetXaxis()->SetTitle("Drift Time [ns]");
         snprintf(title, sizeof(title), "S%i, Drift Pnl %i, MN %i", istation, ihist, panelMNID[istation][ihist]);
         hdriftTime[globalHist]->SetTitle(title);
         hdriftTime[globalHist]->Fit(ftanh,"", "", -8, 10);
         hdriftTime[globalHist]->Draw("same");
         
         Double_t Adrift = ftanh->GetParameter(0);
         Double_t meandrift = ftanh->GetParameter(1);
         Double_t vertOffset = ftanh->GetParameter(2);
         Double_t slope = ftanh->GetParameter(3);
         cout << vertOffset << endl;
         Double_t halfWay = (std::atanh((Adrift-vertOffset)/Adrift) + meandrift)/slope;
         std::cout << "halfway: " << halfWay << std::endl;
         auto legendDrifTtime = new TLegend(0.55,0.64, 0.95,0.9);
         legendDrifTtime->SetTextSize(0.045);
         snprintf(title, sizeof(title), "#splitline{A=%.2f, #mu=%.1f [ns]}{Mid Pt: %.1f [ns]}", Adrift,  meandrift, halfWay);
         legendDrifTtime->AddEntry(hdriftTime[globalHist],title,"l");
         legendDrifTtime->SetFillStyle(0);
         legendDrifTtime->Draw("same");

      }
      cdriftTimePanel->cd(0);
      snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/driftTimePanel_Station%02i_%s_%s.pdf", istation, name, name2);
      cdriftTimePanel->SaveAs(title);
   }

   // Panel charge / energy-deposition histograms.  These are filled in the
   // analyzer with tshs._edep * 1000 and use the same global panel indexing
   // convention as the other 216 per-panel histograms.
   // No charge fit is performed here: use the histogram mean directly.
   TH1F *hEDepPanel[nPanelHistograms] = {};
   Double_t meanCharge[nPanelHistograms];
   Double_t chargeEntries[nPanelHistograms];
   for (Int_t ibin = 0; ibin < nPanelHistograms; ibin++){
      meanCharge[ibin] = NAN;
      chargeEntries[ibin] = 0.0;
   }

   for (Int_t istation = 0; istation < nStations; istation++){
      snprintf(title, sizeof(title), "cChargePanel_Station%02i", istation);
      TCanvas *cChargePanel = new TCanvas(title, title, 0, 0, 1350, 1500);
      cChargePanel->Divide(3, 4);

      for (Int_t ihist = 0; ihist < nPanels; ihist++){
         const Int_t globalHist = istation * nPanels + ihist;
         cChargePanel->cd(ihist + 1);

         snprintf(title, sizeof(title), "a1/hEDepPanel%i", globalHist);
         f1->GetObject(title, hEDepPanel[globalHist]);
         if (!hEDepPanel[globalHist]) continue;

         hEDepPanel[globalHist]->SetLineColor(kBlack);
         hEDepPanel[globalHist]->GetXaxis()->SetTitle("tshs._edep #times 1000");
         hEDepPanel[globalHist]->GetYaxis()->SetTitle("Entries");

         meanCharge[globalHist] = hEDepPanel[globalHist]->GetMean();
         chargeEntries[globalHist] = hEDepPanel[globalHist]->GetEntries();

         snprintf(title, sizeof(title),
                  "S%i, Charge Pnl %i, MN %i",
                  istation, ihist, panelMNID[istation][ihist]);
         hEDepPanel[globalHist]->SetTitle(title);
         hEDepPanel[globalHist]->Draw("hist");

         auto legendCharge = new TLegend(0.55, 0.76, 0.95, 0.90);
         legendCharge->SetTextSize(0.045);
         snprintf(title, sizeof(title),
                  "#splitline{Mean charge = %.3f}{Entries = %.0f}",
                  meanCharge[globalHist], chargeEntries[globalHist]);
         legendCharge->AddEntry(hEDepPanel[globalHist], title, "l");
         legendCharge->SetFillStyle(0);
         legendCharge->Draw("same");
      }

      snprintf(title, sizeof(title),
               "/pnfs/mu2e/scratch/users/dpalo/pdfs/ChargePanel_Station%02i_%s_%s.pdf",
               istation, name, name2);
      cChargePanel->SaveAs(title);
   }

   // One CSV row per tracker panel using the histogram mean directly.
   snprintf(title, sizeof(title), "%s/PanelCharge_%s_%s.csv", csvFolder, name, name2);
   std::ofstream chargeFile(title);
   chargeFile << "StationID,PanelID,MNID,AverageCharge,Entries" << std::endl;
   for (Int_t panelID = 0; panelID < nPanelHistograms; panelID++){
      const Int_t stationID = panelID / nPanels;
      const Int_t panelInStationID = panelID % nPanels;
      chargeFile << stationID << ","
                 << panelInStationID << ","
                 << panelMNID[stationID][panelInStationID] << ","
                 << meanCharge[panelID] << ","
                 << chargeEntries[panelID] << std::endl;
   }
   chargeFile.close();

   

   // ------------------------------------------------------------------
   // Per-panel global-Z sensitivity distributions.
   // hZSensitivityPanelN is filled in TrackAnalyzer before the Z-sensitivity
   // cut, so these plots show the full |(wireDir x trackDir)_unit,Z| range.
   // No fit is applied.
   // ------------------------------------------------------------------
   TH1F *hZSensitivityPanel[nPanelHistograms] = {};

   for (Int_t istation = 0; istation < nStations; ++istation){
      snprintf(title, sizeof(title), "cZSensitivityPanel_Station%02i", istation);
      TCanvas *cZSensitivityPanel = new TCanvas(title, title, 0, 0, 1350, 1500);
      cZSensitivityPanel->Divide(3, 4);

      for (Int_t ipanel = 0; ipanel < nPanels; ++ipanel){
         const Int_t globalHist = istation * nPanels + ipanel;
         cZSensitivityPanel->cd(ipanel + 1);

         snprintf(title, sizeof(title), "a1/hZSensitivityPanel%i", globalHist);
         f1->GetObject(title, hZSensitivityPanel[globalHist]);
         if (!hZSensitivityPanel[globalHist]) continue;

         hZSensitivityPanel[globalHist]->SetLineColor(kBlack);
         hZSensitivityPanel[globalHist]->GetXaxis()->SetTitle("Global Z Sensitivity");
         hZSensitivityPanel[globalHist]->GetYaxis()->SetTitle("Entries");
         hZSensitivityPanel[globalHist]->GetXaxis()->SetRangeUser(0., 1.);
         snprintf(title, sizeof(title), "S%i, Z Sensitivity Pnl %i, MN %i",
                  istation, ipanel, panelMNID[istation][ipanel]);
         hZSensitivityPanel[globalHist]->SetTitle(title);
         hZSensitivityPanel[globalHist]->Draw("hist");
      }

      snprintf(title, sizeof(title),
               "/pnfs/mu2e/scratch/users/dpalo/pdfs/ZSensitivityPanel_Station%02i_%s_%s.pdf",
               istation, name, name2);
      cZSensitivityPanel->SaveAs(title);
   }


   // ------------------------------------------------------------------
   // Station-level GLOBAL-Z residual diagnostics added in TrackAnalyzer.
   //
   // hZErrorVsXStationN : <Delta Z_global> vs global track X
   // hZErrorVsYStationN : <Delta Z_global> vs global track Y
   // hZErrorXYStationN  : <Delta Z_global> in coarse global X-Y bins
   //                      (TProfile2D bin content = average Delta Z).
   // ------------------------------------------------------------------
   TProfile *hZErrorVsXStation[nStations] = {};
   TProfile *hZErrorVsYStation[nStations] = {};
   TProfile2D *hZErrorXYStation[nStations] = {};

   // Linear-fit outputs for <Delta Z_global> vs global Y.
   // Defaults remain zero for stations with no profile/data.
   Double_t zYOffsetMM[nStations] = {};
   Double_t zYRotXRad[nStations] = {};

   for (Int_t istation = 0; istation < nStations; ++istation){
      // Z error versus global X.
      snprintf(title, sizeof(title), "a1/hZErrorVsXStation%i", istation);
      f1->GetObject(title, hZErrorVsXStation[istation]);
      if (hZErrorVsXStation[istation]){
         snprintf(title, sizeof(title), "cZErrorVsXStation%02i", istation);
         TCanvas *cZErrorVsXStation = new TCanvas(title, title, 0, 0, 1000, 800);

         hZErrorVsXStation[istation]->SetLineColor(kBlack);
         hZErrorVsXStation[istation]->SetMarkerStyle(20);
         hZErrorVsXStation[istation]->SetMarkerSize(0.7);
         hZErrorVsXStation[istation]->GetXaxis()->SetTitle("Global X [mm]");
         hZErrorVsXStation[istation]->GetYaxis()->SetTitle("<#Delta Z_{global}> [#mu m]");
         hZErrorVsXStation[istation]->GetXaxis()->SetRangeUser(-1000., 1000.);
         hZErrorVsXStation[istation]->GetYaxis()->SetRangeUser(-1000., 1000.);
         snprintf(title, sizeof(title), "Station %i, Global Z Error vs Global X", istation);
         hZErrorVsXStation[istation]->SetTitle(title);
         hZErrorVsXStation[istation]->Draw("E1");

         snprintf(title, sizeof(title),
                  "/pnfs/mu2e/scratch/users/dpalo/pdfs/ZErrorVsX_Station%02i_%s_%s.pdf",
                  istation, name, name2);
         cZErrorVsXStation->SaveAs(title);
      }

      // Z error versus global Y.
      snprintf(title, sizeof(title), "a1/hZErrorVsYStation%i", istation);
      f1->GetObject(title, hZErrorVsYStation[istation]);
      if (hZErrorVsYStation[istation]){
         snprintf(title, sizeof(title), "cZErrorVsYStation%02i", istation);
         TCanvas *cZErrorVsYStation = new TCanvas(title, title, 0, 0, 1000, 800);

         hZErrorVsYStation[istation]->SetLineColor(kBlack);
         hZErrorVsYStation[istation]->SetMarkerStyle(20);
         hZErrorVsYStation[istation]->SetMarkerSize(0.7);
         hZErrorVsYStation[istation]->GetXaxis()->SetTitle("Global Y [mm]");
         hZErrorVsYStation[istation]->GetYaxis()->SetTitle("<#Delta Z_{global}> [#mu m]");
         hZErrorVsYStation[istation]->GetXaxis()->SetRangeUser(-1000., 1000.);
         hZErrorVsYStation[istation]->GetYaxis()->SetRangeUser(-1000., 1000.);
         snprintf(title, sizeof(title), "Station %i, Global Z Error vs Global Y", istation);
         hZErrorVsYStation[istation]->SetTitle(title);
         hZErrorVsYStation[istation]->Draw("E1");

         // Fit DeltaZ_global [um] = intercept + slope * Y_global [mm].
         // Convert intercept um -> mm and slope um/mm -> rad (small-angle).
         if (hZErrorVsYStation[istation]->GetEntries() > 0){
            TF1 *fZvsY = new TF1(Form("fZvsY_Station%i", istation), "[0]+[1]*x", -1000., 1000.);
            hZErrorVsYStation[istation]->Fit(fZvsY, "", "", -1000., 1000.);
            zYOffsetMM[istation] = fZvsY->GetParameter(0) / 1000.0;
            zYRotXRad[istation] = fZvsY->GetParameter(1) / 1000.0;
         }

         snprintf(title, sizeof(title),
                  "/pnfs/mu2e/scratch/users/dpalo/pdfs/ZErrorVsY_Station%02i_%s_%s.pdf",
                  istation, name, name2);
         cZErrorVsYStation->SaveAs(title);
      }

      // Coarse 20 x 20 X-Y map.  TProfile2D stores mean global-Z residual
      // in each X-Y bin, so COLZ directly displays the requested average error.
      snprintf(title, sizeof(title), "a1/hZErrorXYStation%i", istation);
      f1->GetObject(title, hZErrorXYStation[istation]);
      if (hZErrorXYStation[istation]){
         snprintf(title, sizeof(title), "cZErrorXYStation%02i", istation);
         TCanvas *cZErrorXYStation = new TCanvas(title, title, 0, 0, 1000, 900);

         hZErrorXYStation[istation]->GetXaxis()->SetTitle("Global X [mm]");
         hZErrorXYStation[istation]->GetYaxis()->SetTitle("Global Y [mm]");
         hZErrorXYStation[istation]->GetZaxis()->SetTitle("<#Delta Z_{global}> [#mu m]");
         hZErrorXYStation[istation]->GetXaxis()->SetRangeUser(-1000., 1000.);
         hZErrorXYStation[istation]->GetYaxis()->SetRangeUser(-1000., 1000.);
         hZErrorXYStation[istation]->GetZaxis()->SetRangeUser(-1000., 1000.);
         snprintf(title, sizeof(title), "Station %i, Average Global Z Error vs Global X-Y", istation);
         hZErrorXYStation[istation]->SetTitle(title);
         hZErrorXYStation[istation]->Draw("COLZ");

         snprintf(title, sizeof(title),
                  "/pnfs/mu2e/scratch/users/dpalo/pdfs/ZErrorXY_Station%02i_%s_%s.pdf",
                  istation, name, name2);
         cZErrorXYStation->SaveAs(title);
      }
   }

   // One row per station.  Missing/no-data stations remain all zeros.
   snprintf(title, sizeof(title), "%s/ZErrorVsYFit_%s_%s.csv", csvFolder, name, name2);
   std::ofstream zYFitFile(title);
   zYFitFile << "Station,dZ_mm,rotX_rad" << std::endl;
   for (Int_t istation = 0; istation < nStations; ++istation){
      zYFitFile << istation << ", "
                << zYOffsetMM[istation] << ", "
                << zYRotXRad[istation] << std::endl;
   }
   zYFitFile.close();

// Per-panel DOCA residual diagnostics. The analyzer stores the mean squared
   // residual in TProfiles; convert each bin to sqrt(<residual^2>) here.
   TProfile *hDocaErrorUDOCAPanel[nPanelHistograms] = {};
   TProfile *hDocaErrorSquaredUDOCAPanel[nPanelHistograms] = {};
   TH1F *hUDOCAEntriesPanel[nPanelHistograms] = {};
   TProfile *hDocaErrorSquaredUPosPanel[nPanelHistograms] = {};
   TH1F *hUPosEntriesPanel[nPanelHistograms] = {};
   TH1D *hDocaRMSUDOCAPanel[nPanelHistograms] = {};
   TH1D *hDocaRMSUPosPanel[nPanelHistograms] = {};

   for (Int_t istation = 0; istation < nStations; istation++){
      snprintf(title, sizeof(title), "cDocaErrorUDOCAPanel_Station%02i", istation);
      TCanvas *cDocaErrorUDOCA = new TCanvas(title, title, 0, 0, 1350, 1500);
      cDocaErrorUDOCA->Divide(3, 4);

      snprintf(title, sizeof(title), "cDocaRMSUDOCAPanel_Station%02i", istation);
      TCanvas *cDocaRMSUDOCA = new TCanvas(title, title, 0, 0, 1350, 1500);
      cDocaRMSUDOCA->Divide(3, 4);

      snprintf(title, sizeof(title), "cUDOCAEntriesPanel_Station%02i", istation);
      TCanvas *cUDOCAEntries = new TCanvas(title, title, 0, 0, 1350, 1500);
      cUDOCAEntries->Divide(3, 4);

      snprintf(title, sizeof(title), "cDocaRMSUPosPanel_Station%02i", istation);
      TCanvas *cDocaRMSUPos = new TCanvas(title, title, 0, 0, 1350, 1500);
      cDocaRMSUPos->Divide(3, 4);

      snprintf(title, sizeof(title), "cUPosEntriesPanel_Station%02i", istation);
      TCanvas *cUPosEntries = new TCanvas(title, title, 0, 0, 1350, 1500);
      cUPosEntries->Divide(3, 4);

      for (Int_t ipanel = 0; ipanel < nPanels; ipanel++){
         const Int_t globalHist = istation * nPanels + ipanel;

         snprintf(title, sizeof(title), "a1/hDocaErrorUDOCAPanel%i", globalHist);
         f1->GetObject(title, hDocaErrorUDOCAPanel[globalHist]);
         if (hDocaErrorUDOCAPanel[globalHist]){
            cDocaErrorUDOCA->cd(ipanel + 1);
            hDocaErrorUDOCAPanel[globalHist]->SetLineColor(kBlack);
            hDocaErrorUDOCAPanel[globalHist]->SetMarkerStyle(20);
            hDocaErrorUDOCAPanel[globalHist]->SetMarkerSize(0.55);
            hDocaErrorUDOCAPanel[globalHist]->GetXaxis()->SetTitle("Track UDOCA [mm]");
            hDocaErrorUDOCAPanel[globalHist]->GetYaxis()->SetTitle("<#DeltaDOCA> [#mu m]");
            hDocaErrorUDOCAPanel[globalHist]->GetYaxis()->SetRangeUser(-500., 500.);
            snprintf(title, sizeof(title), "S%i, Signed DOCA Pnl %i, MN %i",
                     istation, ipanel, panelMNID[istation][ipanel]);
            hDocaErrorUDOCAPanel[globalHist]->SetTitle(title);
            hDocaErrorUDOCAPanel[globalHist]->Draw("E1");
         }

         snprintf(title, sizeof(title), "a1/hDocaErrorSquaredUDOCAPanel%i", globalHist);
         f1->GetObject(title, hDocaErrorSquaredUDOCAPanel[globalHist]);
         if (hDocaErrorSquaredUDOCAPanel[globalHist]){
            snprintf(title, sizeof(title), "hDocaRMSUDOCAPanel%i", globalHist);
            hDocaRMSUDOCAPanel[globalHist] = hDocaErrorSquaredUDOCAPanel[globalHist]->ProjectionX(title);
            for (Int_t ibin = 1; ibin <= hDocaRMSUDOCAPanel[globalHist]->GetNbinsX(); ibin++){
               const Double_t meanSquare = hDocaErrorSquaredUDOCAPanel[globalHist]->GetBinContent(ibin);
               const Double_t meanSquareError = hDocaErrorSquaredUDOCAPanel[globalHist]->GetBinError(ibin);
               if (meanSquare > 0.0){
                  const Double_t rms = std::sqrt(meanSquare);
                  hDocaRMSUDOCAPanel[globalHist]->SetBinContent(ibin, rms);
                  hDocaRMSUDOCAPanel[globalHist]->SetBinError(ibin, meanSquareError / (2.0 * rms));
               } else {
                  hDocaRMSUDOCAPanel[globalHist]->SetBinContent(ibin, 0.0);
                  hDocaRMSUDOCAPanel[globalHist]->SetBinError(ibin, 0.0);
               }
            }
            cDocaRMSUDOCA->cd(ipanel + 1);
            hDocaRMSUDOCAPanel[globalHist]->SetLineColor(kBlack);
            hDocaRMSUDOCAPanel[globalHist]->SetMarkerStyle(20);
            hDocaRMSUDOCAPanel[globalHist]->SetMarkerSize(0.55);
            hDocaRMSUDOCAPanel[globalHist]->GetXaxis()->SetTitle("Track UDOCA [mm]");
            hDocaRMSUDOCAPanel[globalHist]->GetYaxis()->SetTitle("#sqrt{<#DeltaDOCA^{2}>} [#mu m]");
            hDocaRMSUDOCAPanel[globalHist]->GetYaxis()->SetRangeUser(0., 600.);
            snprintf(title, sizeof(title), "S%i, DOCA RMS Pnl %i, MN %i",
                     istation, ipanel, panelMNID[istation][ipanel]);
            hDocaRMSUDOCAPanel[globalHist]->SetTitle(title);
            hDocaRMSUDOCAPanel[globalHist]->Draw("E1");
         }

         snprintf(title, sizeof(title), "a1/hUDOCAEntriesPanel%i", globalHist);
         f1->GetObject(title, hUDOCAEntriesPanel[globalHist]);
         if (hUDOCAEntriesPanel[globalHist]){
            cUDOCAEntries->cd(ipanel + 1);
            hUDOCAEntriesPanel[globalHist]->SetLineColor(kBlack);
            hUDOCAEntriesPanel[globalHist]->GetXaxis()->SetTitle("Track UDOCA [mm]");
            hUDOCAEntriesPanel[globalHist]->GetYaxis()->SetTitle("Entries");
            snprintf(title, sizeof(title), "S%i, UDOCA Entries Pnl %i, MN %i",
                     istation, ipanel, panelMNID[istation][ipanel]);
            hUDOCAEntriesPanel[globalHist]->SetTitle(title);
            hUDOCAEntriesPanel[globalHist]->Draw("hist");
         }

         snprintf(title, sizeof(title), "a1/hDocaErrorSquaredUPosPanel%i", globalHist);
         f1->GetObject(title, hDocaErrorSquaredUPosPanel[globalHist]);
         if (hDocaErrorSquaredUPosPanel[globalHist]){
            snprintf(title, sizeof(title), "hDocaRMSUPosPanel%i", globalHist);
            hDocaRMSUPosPanel[globalHist] = hDocaErrorSquaredUPosPanel[globalHist]->ProjectionX(title);
            for (Int_t ibin = 1; ibin <= hDocaRMSUPosPanel[globalHist]->GetNbinsX(); ibin++){
               const Double_t meanSquare = hDocaErrorSquaredUPosPanel[globalHist]->GetBinContent(ibin);
               const Double_t meanSquareError = hDocaErrorSquaredUPosPanel[globalHist]->GetBinError(ibin);
               if (meanSquare > 0.0){
                  const Double_t rms = std::sqrt(meanSquare);
                  hDocaRMSUPosPanel[globalHist]->SetBinContent(ibin, rms);
                  hDocaRMSUPosPanel[globalHist]->SetBinError(ibin, meanSquareError / (2.0 * rms));
               } else {
                  hDocaRMSUPosPanel[globalHist]->SetBinContent(ibin, 0.0);
                  hDocaRMSUPosPanel[globalHist]->SetBinError(ibin, 0.0);
               }
            }
            cDocaRMSUPos->cd(ipanel + 1);
            hDocaRMSUPosPanel[globalHist]->SetLineColor(kBlack);
            hDocaRMSUPosPanel[globalHist]->SetMarkerStyle(20);
            hDocaRMSUPosPanel[globalHist]->SetMarkerSize(0.55);
            hDocaRMSUPosPanel[globalHist]->GetXaxis()->SetTitle("tshs._uupos [mm]");
            hDocaRMSUPosPanel[globalHist]->GetYaxis()->SetTitle("#sqrt{<#DeltaDOCA^{2}>} [#mu m]");
            hDocaRMSUPosPanel[globalHist]->GetYaxis()->SetRangeUser(0., 600.);
            snprintf(title, sizeof(title), "S%i, DOCA RMS vs U Pos Pnl %i, MN %i",
                     istation, ipanel, panelMNID[istation][ipanel]);
            hDocaRMSUPosPanel[globalHist]->SetTitle(title);
            hDocaRMSUPosPanel[globalHist]->Draw("E1");
         }

         snprintf(title, sizeof(title), "a1/hUPosEntriesPanel%i", globalHist);
         f1->GetObject(title, hUPosEntriesPanel[globalHist]);
         if (hUPosEntriesPanel[globalHist]){
            cUPosEntries->cd(ipanel + 1);
            hUPosEntriesPanel[globalHist]->SetLineColor(kBlack);
            hUPosEntriesPanel[globalHist]->GetXaxis()->SetTitle("tshs._uupos [mm]");
            hUPosEntriesPanel[globalHist]->GetYaxis()->SetTitle("Entries");
            snprintf(title, sizeof(title), "S%i, U Pos Entries Pnl %i, MN %i",
                     istation, ipanel, panelMNID[istation][ipanel]);
            hUPosEntriesPanel[globalHist]->SetTitle(title);
            hUPosEntriesPanel[globalHist]->Draw("hist");
         }
      }

      snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/DocaErrorUDOCA_Station%02i_%s_%s.pdf", istation, name, name2);
      cDocaErrorUDOCA->SaveAs(title);
      snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/DocaRMSUDOCA_Station%02i_%s_%s.pdf", istation, name, name2);
      cDocaRMSUDOCA->SaveAs(title);
      snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/UDOCAEntries_Station%02i_%s_%s.pdf", istation, name, name2);
      cUDOCAEntries->SaveAs(title);
      snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/DocaRMSUPos_Station%02i_%s_%s.pdf", istation, name, name2);
      cDocaRMSUPos->SaveAs(title);
      snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/UPosEntries_Station%02i_%s_%s.pdf", istation, name, name2);
      cUPosEntries->SaveAs(title);
   }

   TH1F *hDOCAErrorPanel[nPanelHistograms];
    Double_t A1[nPanelHistograms];
    Double_t Sig1[nPanelHistograms];
    Double_t Sig2[nPanelHistograms];
    Double_t Mu1[nPanelHistograms];
    Double_t Mu2[nPanelHistograms];
    for (Int_t ibin = 0; ibin < nPanelHistograms; ibin++){
       A1[ibin] = 0.0;
       Sig1[ibin] = 0.0;
       Sig2[ibin] = 0.0;
       Mu1[ibin] = 0.0;
       Mu2[ibin] = 0.0;
    }

    
    for (Int_t istation = 0; istation < nStations; istation++){
        snprintf(title, sizeof(title), "cDOCAErrorPanel_Station%02i", istation);
        TCanvas *cDOCAErrorPanel = new TCanvas(title, title, 0, 0, 1350, 1500);
        cDOCAErrorPanel->Divide(3, 4);
        
        for (Int_t ihist = 0; ihist < nPanels; ihist++){
            const Int_t globalHist = istation * nPanels + ihist;
            cDOCAErrorPanel->cd(ihist+1);
            snprintf(title, sizeof(title), "a1/hDOCAErrorPanel%i", globalHist);
            f1->GetObject(title, hDOCAErrorPanel[globalHist]);
            if (!hDOCAErrorPanel[globalHist]) continue;
            hDOCAErrorPanel[globalHist]->SetLineColor(kBlack);
            hDOCAErrorPanel[globalHist]->GetXaxis()->SetTitle("DOCA Error [um]");
            snprintf(title, sizeof(title), "S%i, DOCA Err Pnl %i, MN %i", istation, ihist, panelMNID[istation][ihist]);
            hDOCAErrorPanel[globalHist]->SetTitle(title);
            hDOCAErrorPanel[globalHist]->GetXaxis()->SetRangeUser(-1250, 1250);
            hDOCAErrorPanel[globalHist]->Draw();

            const Double_t entriesDOCA = hDOCAErrorPanel[globalHist]->GetEntries();
            if (entriesDOCA >= 500){
               hDOCAErrorPanel[globalHist]->Fit(fDOCA, "", "", -1500, 1500);

               Double_t meanCore2 = fDOCA->GetParameter(1);
               Double_t sigmaCore2 = fDOCA->GetParameter(2);
               Double_t fractionTail2 = fDOCA->GetParameter(5);
               Double_t meanTail2 = fDOCA->GetParameter(3);
               Double_t sigmaTail2 = fDOCA->GetParameter(4);

               const Bool_t docaOutlier =
                   (sigmaCore2 > 400.0) ||
                   (sigmaTail2 > 1000.0) ||
                   (std::abs(meanCore2) > 150.0) ||
                   (std::abs(meanTail2) > 300.0);

               if (docaOutlier){
                  snprintf(title, sizeof(title),
                           "S%i, DOCA Err Pnl %i, MN %i [OUTLIER]",
                           istation, ihist, panelMNID[istation][ihist]);
                  hDOCAErrorPanel[globalHist]->SetTitle(title);
               }

               A1[globalHist] = fDOCA->GetParameter(0);
               Sig1[globalHist] = sigmaCore2;
               Sig2[globalHist] = sigmaTail2;
               Mu1[globalHist] = meanCore2;
               Mu2[globalHist] = meanTail2;

               auto legendDOCAPanel = new TLegend(0.55,0.64, 0.95,0.9);
               legendDOCAPanel->SetTextSize(0.045);
               snprintf(title, sizeof(title),
                        "#splitline{#sigma_{1}=%.2f, #mu_{1}=%.1f [um]}{#splitline{#sigma_{2}=%.2f, #mu_{2}=%.1f [um]}{A_{2}/A_{1}=%.2f}}",
                        sigmaCore2, meanCore2, sigmaTail2, meanTail2, fractionTail2);
               legendDOCAPanel->AddEntry(hDOCAErrorPanel[globalHist], title, "l");
               legendDOCAPanel->SetFillStyle(0);
               legendDOCAPanel->Draw("same");
            } else {
               std::cout << "Skipping DOCA panel fit for global panel " << globalHist
                         << ": " << entriesDOCA << " entries." << std::endl;
            }
            
        }
        
        snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/docaerrorPanel_Station%02i_%s_%s.pdf", istation, name, name2);
        cDOCAErrorPanel->SaveAs(title);
    }
    snprintf(title, sizeof(title), "%s/DOCA_Histograms_%s_%s.csv", csvFolder, name, name2);
    std::ofstream myfile(title);
    for (Int_t ibin = 0; ibin < nPanelHistograms; ibin++){
       const Int_t istationCSV = ibin / nPanels;
       const Int_t ipanelCSV = ibin % nPanels;
       myfile << ibin << ", "
              << A1[ibin] << ", "
              << Sig1[ibin] << ", "
              << Sig2[ibin] << ", "
              << Mu1[ibin] << ", "
              << Mu2[ibin] << ", "
              << panelMNID[istationCSV][ipanelCSV] << std::endl;
    }
    
    
   TH1F *hULocal[nPanelHistograms];

   for (Int_t istation = 0; istation < nStations; istation++){
      snprintf(title, sizeof(title), "cULocalPanel_Station%02i", istation);
      TCanvas *cULocalPanel = new TCanvas(title, title, 0, 0, 1350, 1500);
      cULocalPanel->Divide(3, 4);

      for (Int_t ihist = 0; ihist < nPanels; ihist++){
         const Int_t globalHist = istation * nPanels + ihist;
         cULocalPanel->cd(ihist+1);
         snprintf(title, sizeof(title), "a1/hULocal%i", globalHist);
         f1->GetObject(title, hULocal[globalHist]);
         if (!hULocal[globalHist]) continue;
         hULocal[globalHist]->SetLineColor(kBlack);
         hULocal[globalHist]->GetXaxis()->SetTitle("ULocal[ mm]");
         snprintf(title, sizeof(title), "S%i, U Local Pnl %i, MN %i", istation, ihist, panelMNID[istation][ihist]);
         hULocal[globalHist]->SetTitle(title);
         hULocal[globalHist]->Draw("same");
      }

      cULocalPanel->cd(0);
      snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/ULocalPanel_Station%02i_%s_%s.pdf", istation, name, name2);
      cULocalPanel->SaveAs(title);
   }

   TProfile *hdVULocal[nPanelHistograms];

   for (Int_t istation = 0; istation < nStations; istation++){
      snprintf(title, sizeof(title), "cdVULocalPanel_Station%02i", istation);
      TCanvas *cdVULocalPanel = new TCanvas(title, title, 0, 0, 1350, 1500);
      cdVULocalPanel->Divide(3, 4);

      for (Int_t ihist = 0; ihist < nPanels; ihist++){
         const Int_t globalHist = istation * nPanels + ihist;
         cdVULocalPanel->cd(ihist+1);
         snprintf(title, sizeof(title), "a1/hdVULocal%i", globalHist);
         f1->GetObject(title, hdVULocal[globalHist]);
         if (!hdVULocal[globalHist]) continue;
         hdVULocal[globalHist]->SetLineColor(kBlack);
         hdVULocal[globalHist]->GetXaxis()->SetTitle("U Local [ mm]");
         snprintf(title, sizeof(title), "S%i, dV:U Pnl %i, MN %i", istation, ihist, panelMNID[istation][ihist]);
         hdVULocal[globalHist]->SetTitle(title);
      hdVULocal[globalHist]->GetYaxis()->SetRangeUser(-300, 300);
         hdVULocal[globalHist]->Draw("same");
      }

      cdVULocalPanel->cd(0);
      snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/dVULocalPanel_Station%02i_%s_%s.pdf", istation, name, name2);
      cdVULocalPanel->SaveAs(title);
   }

   // Local V residual for the full panel, fit with the same double-Gaussian
   // model used for DOCA/X/Y/Z. Histograms with <500 entries are not fit and
   // their fit parameters are written as zeros.
   TH1F *hVErrorPanel[nPanelHistograms] = {};
   Double_t VA1[nPanelHistograms];
   Double_t VSig1[nPanelHistograms];
   Double_t VSig2[nPanelHistograms];
   Double_t VMu1[nPanelHistograms];
   Double_t VMu2[nPanelHistograms];
   for (Int_t ibin = 0; ibin < nPanelHistograms; ibin++){
      VA1[ibin] = 0.0;
      VSig1[ibin] = 0.0;
      VSig2[ibin] = 0.0;
      VMu1[ibin] = 0.0;
      VMu2[ibin] = 0.0;
   }

   for (Int_t istation = 0; istation < nStations; istation++){
      snprintf(title, sizeof(title), "cVErrorPanelPanel_Station%02i", istation);
      TCanvas *cVErrorPanelPanel = new TCanvas(title, title, 0, 0, 1350, 1500);
      cVErrorPanelPanel->Divide(3, 4);

      for (Int_t ihist = 0; ihist < nPanels; ihist++){
         const Int_t globalHist = istation * nPanels + ihist;
         cVErrorPanelPanel->cd(ihist+1);
         snprintf(title, sizeof(title), "a1/hVErrorPanel%i", globalHist);
         f1->GetObject(title, hVErrorPanel[globalHist]);
         if (!hVErrorPanel[globalHist]) continue;

         hVErrorPanel[globalHist]->SetLineColor(kBlack);
         hVErrorPanel[globalHist]->GetXaxis()->SetTitle("V Error [um]");
         hVErrorPanel[globalHist]->GetXaxis()->SetRangeUser(-1250, 1250);
         snprintf(title, sizeof(title), "S%i, V Err Pnl %i, MN %i",
                  istation, ihist, panelMNID[istation][ihist]);
         hVErrorPanel[globalHist]->SetTitle(title);
         hVErrorPanel[globalHist]->Draw();

         const Double_t entriesV = hVErrorPanel[globalHist]->GetEntries();
         if (entriesV >= 500){
            hVErrorPanel[globalHist]->Fit(fDOCA, "", "", -1500, 1500);

            const Double_t meanCoreV = fDOCA->GetParameter(1);
            const Double_t sigmaCoreV = fDOCA->GetParameter(2);
            const Double_t fractionTailV = fDOCA->GetParameter(5);
            const Double_t meanTailV = fDOCA->GetParameter(3);
            const Double_t sigmaTailV = fDOCA->GetParameter(4);

            VA1[globalHist] = fDOCA->GetParameter(0);
            VSig1[globalHist] = sigmaCoreV;
            VSig2[globalHist] = sigmaTailV;
            VMu1[globalHist] = meanCoreV;
            VMu2[globalHist] = meanTailV;

            auto legendVPanel = new TLegend(0.55, 0.64, 0.95, 0.9);
            legendVPanel->SetTextSize(0.045);
            snprintf(title, sizeof(title),
                     "#splitline{#sigma_{1}=%.2f, #mu_{1}=%.1f [um]}{#splitline{#sigma_{2}=%.2f, #mu_{2}=%.1f [um]}{A_{2}/A_{1}=%.2f}}",
                     sigmaCoreV, meanCoreV, sigmaTailV, meanTailV, fractionTailV);
            legendVPanel->AddEntry(hVErrorPanel[globalHist], title, "l");
            legendVPanel->SetFillStyle(0);
            legendVPanel->Draw("same");
         } else {
            std::cout << "Skipping V-error panel fit for global panel " << globalHist
                      << ": " << entriesV << " entries." << std::endl;
         }
      }

      snprintf(title, sizeof(title),
               "/pnfs/mu2e/scratch/users/dpalo/pdfs/VErrorPanelPanel_Station%02i_%s_%s.pdf",
               istation, name, name2);
      cVErrorPanelPanel->SaveAs(title);
   }

   snprintf(title, sizeof(title), "%s/VErrorPanel_%s_%s.csv", csvFolder, name, name2);
   std::ofstream vErrorFile(title);
   for (Int_t ibin = 0; ibin < nPanelHistograms; ibin++){
      const Int_t istationCSV = ibin / nPanels;
      const Int_t ipanelCSV = ibin % nPanels;
      vErrorFile << ibin << ", "
                 << VA1[ibin] << ", "
                 << VSig1[ibin] << ", "
                 << VSig2[ibin] << ", "
                 << VMu1[ibin] << ", "
                 << VMu2[ibin] << ", "
                 << panelMNID[istationCSV][ipanelCSV] << std::endl;
   }

   // Global X residual for the full panel, fit with the same double-Gaussian
   // model used for the DOCA-error panel histograms.
   TH1F *hXErrorPanel[nPanelHistograms] = {};
   Double_t XA1[nPanelHistograms];
   Double_t XSig1[nPanelHistograms];
   Double_t XSig2[nPanelHistograms];
   Double_t XMu1[nPanelHistograms];
   Double_t XMu2[nPanelHistograms];
   for (Int_t ibin = 0; ibin < nPanelHistograms; ibin++){
      XA1[ibin] = 0.0;
      XSig1[ibin] = 0.0;
      XSig2[ibin] = 0.0;
      XMu1[ibin] = 0.0;
      XMu2[ibin] = 0.0;
   }

   for (Int_t istation = 0; istation < nStations; istation++){
      snprintf(title, sizeof(title), "cXErrorPanel_Station%02i", istation);
      TCanvas *cXErrorPanel = new TCanvas(title, title, 0, 0, 1350, 1500);
      cXErrorPanel->Divide(3, 4);

      for (Int_t ihist = 0; ihist < nPanels; ihist++){
         const Int_t globalHist = istation * nPanels + ihist;
         cXErrorPanel->cd(ihist+1);
         snprintf(title, sizeof(title), "a1/hXErrorPanel%i", globalHist);
         f1->GetObject(title, hXErrorPanel[globalHist]);
         if (!hXErrorPanel[globalHist]) continue;

         hXErrorPanel[globalHist]->SetLineColor(kBlack);
         hXErrorPanel[globalHist]->GetXaxis()->SetTitle("Global X Error [um]");
         snprintf(title, sizeof(title), "S%i, Gl X Err Pnl %i, MN %i", istation, ihist, panelMNID[istation][ihist]);
         hXErrorPanel[globalHist]->SetTitle(title);
         hXErrorPanel[globalHist]->GetXaxis()->SetRangeUser(-1250, 1250);
         hXErrorPanel[globalHist]->Draw("same");

         const Double_t entriesX = hXErrorPanel[globalHist]->GetEntries();
         if (entriesX >= 500){
         hXErrorPanel[globalHist]->Fit(fDOCA, "", "", -1500, 1500);
   
            const Double_t meanCoreX = fDOCA->GetParameter(1);
            const Double_t sigmaCoreX = fDOCA->GetParameter(2);
            const Double_t fractionTailX = fDOCA->GetParameter(5);
            const Double_t meanTailX = fDOCA->GetParameter(3);
            const Double_t sigmaTailX = fDOCA->GetParameter(4);
   
            XA1[globalHist] = fDOCA->GetParameter(0);
            XSig1[globalHist] = sigmaCoreX;
            XSig2[globalHist] = sigmaTailX;
            XMu1[globalHist] = meanCoreX;
            XMu2[globalHist] = meanTailX;
   
            auto legendXPanel = new TLegend(0.55, 0.64, 0.95, 0.9);
            legendXPanel->SetTextSize(0.045);
            snprintf(title, sizeof(title),
                     "#splitline{#sigma_{1}=%.2f, #mu_{1}=%.1f [um]}{#splitline{#sigma_{2}=%.2f, #mu_{2}=%.1f [um]}{A_{2}/A_{1}=%.2f}}",
                     sigmaCoreX, meanCoreX, sigmaTailX, meanTailX, fractionTailX);
            legendXPanel->AddEntry(hXErrorPanel[globalHist], title, "l");
            legendXPanel->SetFillStyle(0);
            legendXPanel->Draw("same");
         } else {
            std::cout << "Skipping X-error panel fit for global panel " << globalHist
                      << ": " << entriesX << " entries." << std::endl;
         }
      }

      snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/XErrorPanel_Station%02i_%s_%s.pdf", istation, name, name2);
      cXErrorPanel->SaveAs(title);
   }

   // Keep the same six-column format as DOCA_Histograms: global panel index,
   // A1, sigma1, sigma2, mu1, mu2.
   snprintf(title, sizeof(title), "%s/XErrorPanel_%s_%s.csv", csvFolder, name, name2);
   std::ofstream xErrorFile(title);
   for (Int_t ibin = 0; ibin < nPanelHistograms; ibin++){
      const Int_t istationCSV = ibin / nPanels;

      const Int_t ipanelCSV = ibin % nPanels;

      xErrorFile << ibin << ", "

                 << XA1[ibin] << ", "

                 << XSig1[ibin] << ", "

                 << XSig2[ibin] << ", "

                 << XMu1[ibin] << ", "

                 << XMu2[ibin] << ", "

                 << panelMNID[istationCSV][ipanelCSV] << std::endl;
   }

   // Global Y residual for the full panel, fit with the same double-Gaussian
   // model used for XErrorPanel and DOCA-error panel histograms.
   TH1F *hYErrorPanel[nPanelHistograms] = {};
   Double_t YA1[nPanelHistograms];
   Double_t YSig1[nPanelHistograms];
   Double_t YSig2[nPanelHistograms];
   Double_t YMu1[nPanelHistograms];
   Double_t YMu2[nPanelHistograms];
   for (Int_t ibin = 0; ibin < nPanelHistograms; ibin++){
      YA1[ibin] = 0.0;
      YSig1[ibin] = 0.0;
      YSig2[ibin] = 0.0;
      YMu1[ibin] = 0.0;
      YMu2[ibin] = 0.0;
   }

   for (Int_t istation = 0; istation < nStations; istation++){
      snprintf(title, sizeof(title), "cYErrorPanel_Station%02i", istation);
      TCanvas *cYErrorPanel = new TCanvas(title, title, 0, 0, 1350, 1500);
      cYErrorPanel->Divide(3, 4);

      for (Int_t ihist = 0; ihist < nPanels; ihist++){
         const Int_t globalHist = istation * nPanels + ihist;
         cYErrorPanel->cd(ihist+1);
         snprintf(title, sizeof(title), "a1/hYErrorPanel%i", globalHist);
         f1->GetObject(title, hYErrorPanel[globalHist]);
         if (!hYErrorPanel[globalHist]) continue;

         hYErrorPanel[globalHist]->SetLineColor(kBlack);
         hYErrorPanel[globalHist]->GetXaxis()->SetTitle("Global Y Error [um]");
         snprintf(title, sizeof(title), "S%i, Gl Y Err Pnl %i, MN %i", istation, ihist, panelMNID[istation][ihist]);
         hYErrorPanel[globalHist]->SetTitle(title);
         hYErrorPanel[globalHist]->GetXaxis()->SetRangeUser(-1250, 1250);
         hYErrorPanel[globalHist]->Draw("same");

         const Double_t entriesY = hYErrorPanel[globalHist]->GetEntries();
         if (entriesY >= 500){
         hYErrorPanel[globalHist]->Fit(fDOCA, "", "", -1500, 1500);
   
            const Double_t meanCoreY = fDOCA->GetParameter(1);
            const Double_t sigmaCoreY = fDOCA->GetParameter(2);
            const Double_t fractionTailY = fDOCA->GetParameter(5);
            const Double_t meanTailY = fDOCA->GetParameter(3);
            const Double_t sigmaTailY = fDOCA->GetParameter(4);
   
            YA1[globalHist] = fDOCA->GetParameter(0);
            YSig1[globalHist] = sigmaCoreY;
            YSig2[globalHist] = sigmaTailY;
            YMu1[globalHist] = meanCoreY;
            YMu2[globalHist] = meanTailY;
   
            auto legendYPanel = new TLegend(0.55, 0.64, 0.95, 0.9);
            legendYPanel->SetTextSize(0.045);
            snprintf(title, sizeof(title),
                     "#splitline{#sigma_{1}=%.2f, #mu_{1}=%.1f [um]}{#splitline{#sigma_{2}=%.2f, #mu_{2}=%.1f [um]}{A_{2}/A_{1}=%.2f}}",
                     sigmaCoreY, meanCoreY, sigmaTailY, meanTailY, fractionTailY);
            legendYPanel->AddEntry(hYErrorPanel[globalHist], title, "l");
            legendYPanel->SetFillStyle(0);
            legendYPanel->Draw("same");
         } else {
            std::cout << "Skipping Y-error panel fit for global panel " << globalHist
                      << ": " << entriesY << " entries." << std::endl;
         }
      }

      snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/YErrorPanel_Station%02i_%s_%s.pdf", istation, name, name2);
      cYErrorPanel->SaveAs(title);
   }

   // Same six-column format as XErrorPanel: global panel index,
   // A1, sigma1, sigma2, mu1, mu2.
   snprintf(title, sizeof(title), "%s/YErrorPanel_%s_%s.csv", csvFolder, name, name2);
   std::ofstream yErrorFile(title);
   for (Int_t ibin = 0; ibin < nPanelHistograms; ibin++){
      const Int_t istationCSV = ibin / nPanels;

      const Int_t ipanelCSV = ibin % nPanels;

      yErrorFile << ibin << ", "

                 << YA1[ibin] << ", "

                 << YSig1[ibin] << ", "

                 << YSig2[ibin] << ", "

                 << YMu1[ibin] << ", "

                 << YMu2[ibin] << ", "

                 << panelMNID[istationCSV][ipanelCSV] << std::endl;
   }

   // Global Z residual for the full panel, fit with the same double-Gaussian
   // model used for XErrorPanel, YErrorPanel, and DOCA-error panel histograms.
   TH1F *hZErrorPanel[nPanelHistograms] = {};
   Double_t ZA1[nPanelHistograms];
   Double_t ZSig1[nPanelHistograms];
   Double_t ZSig2[nPanelHistograms];
   Double_t ZMu1[nPanelHistograms];
   Double_t ZMu2[nPanelHistograms];
   for (Int_t ibin = 0; ibin < nPanelHistograms; ibin++){
      ZA1[ibin] = 0.0;
      ZSig1[ibin] = 0.0;
      ZSig2[ibin] = 0.0;
      ZMu1[ibin] = 0.0;
      ZMu2[ibin] = 0.0;
   }

   for (Int_t istation = 0; istation < nStations; istation++){
      snprintf(title, sizeof(title), "cZErrorPanel_Station%02i", istation);
      TCanvas *cZErrorPanel = new TCanvas(title, title, 0, 0, 1350, 1500);
      cZErrorPanel->Divide(3, 4);

      for (Int_t ihist = 0; ihist < nPanels; ihist++){
         const Int_t globalHist = istation * nPanels + ihist;
         cZErrorPanel->cd(ihist+1);
         snprintf(title, sizeof(title), "a1/hZErrorPanel%i", globalHist);
         f1->GetObject(title, hZErrorPanel[globalHist]);
         if (!hZErrorPanel[globalHist]) continue;

         hZErrorPanel[globalHist]->SetLineColor(kBlack);
         hZErrorPanel[globalHist]->GetXaxis()->SetTitle("Global Z Error [um]");
         snprintf(title, sizeof(title), "S%i, Gl Z Err Pnl %i, MN %i", istation, ihist, panelMNID[istation][ihist]);
         hZErrorPanel[globalHist]->SetTitle(title);
         hZErrorPanel[globalHist]->GetXaxis()->SetRangeUser(-1250, 1250);
         hZErrorPanel[globalHist]->Draw("same");

         const Double_t entriesZ = hZErrorPanel[globalHist]->GetEntries();
         if (entriesZ >= 500){
         hZErrorPanel[globalHist]->Fit(fDOCA, "", "", -1500, 1500);
   
            const Double_t meanCoreZ = fDOCA->GetParameter(1);
            const Double_t sigmaCoreZ = fDOCA->GetParameter(2);
            const Double_t fractionTailZ = fDOCA->GetParameter(5);
            const Double_t meanTailZ = fDOCA->GetParameter(3);
            const Double_t sigmaTailZ = fDOCA->GetParameter(4);
   
            ZA1[globalHist] = fDOCA->GetParameter(0);
            ZSig1[globalHist] = sigmaCoreZ;
            ZSig2[globalHist] = sigmaTailZ;
            ZMu1[globalHist] = meanCoreZ;
            ZMu2[globalHist] = meanTailZ;
   
            auto legendZPanel = new TLegend(0.55, 0.64, 0.95, 0.9);
            legendZPanel->SetTextSize(0.045);
            snprintf(title, sizeof(title),
                     "#splitline{#sigma_{1}=%.2f, #mu_{1}=%.1f [um]}{#splitline{#sigma_{2}=%.2f, #mu_{2}=%.1f [um]}{A_{2}/A_{1}=%.2f}}",
                     sigmaCoreZ, meanCoreZ, sigmaTailZ, meanTailZ, fractionTailZ);
            legendZPanel->AddEntry(hZErrorPanel[globalHist], title, "l");
            legendZPanel->SetFillStyle(0);
            legendZPanel->Draw("same");
         } else {
            std::cout << "Skipping Z-error panel fit for global panel " << globalHist
                      << ": " << entriesZ << " entries." << std::endl;
         }
      }

      snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/ZErrorPanel_Station%02i_%s_%s.pdf", istation, name, name2);
      cZErrorPanel->SaveAs(title);
   }

   // Same six-column format as XErrorPanel/YErrorPanel:
   // global panel index, A1, sigma1, sigma2, mu1, mu2.
   snprintf(title, sizeof(title), "%s/ZErrorPanel_%s_%s.csv", csvFolder, name, name2);
   std::ofstream zErrorFile(title);
   for (Int_t ibin = 0; ibin < nPanelHistograms; ibin++){
      const Int_t istationCSV = ibin / nPanels;

      const Int_t ipanelCSV = ibin % nPanels;

      zErrorFile << ibin << ", "

                 << ZA1[ibin] << ", "

                 << ZSig1[ibin] << ", "

                 << ZSig2[ibin] << ", "

                 << ZMu1[ibin] << ", "

                 << ZMu2[ibin] << ", "

                 << panelMNID[istationCSV][ipanelCSV] << std::endl;
   }

   // Residual versus straw U position. X/Y are global residuals; V is panel-local.
   struct PanelProfileSpec {
      const char* rootStem;
      const char* canvasStem;
      const char* yTitle;
      const char* plotLabel;
      const char* pdfStem;
   };

   const PanelProfileSpec panelProfiles[] = {
      {"hXErrorZPanelP", "cXErrorZPanelP", "Global X Error [um]", "Global X Error vs U", "XErrorZPanelP"},
      {"hYErrorZPanelP", "cYErrorZPanelP", "Global Y Error [um]", "Global Y Error vs U", "YErrorZPanelP"},
      {"hVErrorZPanelP", "cVErrorZPanelP", "Local V Error [um]",  "Local V Error vs U",  "VErrorZPanelP"}
   };

   for (const auto& spec : panelProfiles){
      TProfile *profile[nPanelHistograms] = {};
      for (Int_t istation = 0; istation < nStations; istation++){
         snprintf(title, sizeof(title), "%s_Station%02i", spec.canvasStem, istation);
         TCanvas *canvas = new TCanvas(title, title, 0, 0, 1350, 1500);
         canvas->Divide(3, 4);

         for (Int_t ihist = 0; ihist < nPanels; ihist++){
            const Int_t globalHist = istation * nPanels + ihist;
            canvas->cd(ihist+1);
            snprintf(title, sizeof(title), "a1/%s%i", spec.rootStem, globalHist);
            f1->GetObject(title, profile[globalHist]);
            if (!profile[globalHist]) continue;
            profile[globalHist]->SetLineColor(kBlack);
            profile[globalHist]->GetXaxis()->SetTitle("Straw U Position [mm]");
            profile[globalHist]->GetYaxis()->SetTitle(spec.yTitle);
            snprintf(title, sizeof(title), "S%i, %s Pnl %i, MN %i", istation, spec.plotLabel, ihist, panelMNID[istation][ihist]);
            profile[globalHist]->SetTitle(title);
	    profile[globalHist]->GetYaxis()->SetRangeUser(-300, 300);
            profile[globalHist]->Draw();
         }

         snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/%s_Station%02i_%s_%s.pdf", spec.pdfStem, istation, name, name2);
         canvas->SaveAs(title);
      }
   }

   TProfile *_hDocaErrorRDrift;
   TProfile *_hDocaErrorUDOCA;
   TProfile *_hDocaErrorCDrift;

   TCanvas *DocaErrorProfile = new TCanvas("DocaErrorProfile", "DocaErrorProfile", 0, 0, 1350, 1500);
   DocaErrorProfile->Divide(1);

   snprintf(title, sizeof(title), "a1/_hDocaErrorRDrift");
   f1->GetObject(title, _hDocaErrorRDrift);
   snprintf(title, sizeof(title), "Doca Residual Profiles");

   _hDocaErrorRDrift->SetLineColor(kBlack);
   _hDocaErrorRDrift->GetXaxis()->SetTitle("DOCA [mm]");
   _hDocaErrorRDrift->SetTitle(title);

   _hDocaErrorRDrift->GetYaxis()->SetRangeUser(-200, 200);

   _hDocaErrorRDrift->Draw();
   
   snprintf(title, sizeof(title), "a1/_hDocaErrorCDrift");
   f1->GetObject(title, _hDocaErrorCDrift);
   _hDocaErrorCDrift->SetLineColor(kRed);
   _hDocaErrorCDrift->GetXaxis()->SetTitle("DOCA [mm]");
   _hDocaErrorCDrift->SetTitle(title);
   _hDocaErrorCDrift->Draw("SAME");
   
   snprintf(title, sizeof(title), "a1/_hDocaErrorUDOCA");
   f1->GetObject(title, _hDocaErrorUDOCA);
   _hDocaErrorUDOCA->SetLineColor(kGreen);
   _hDocaErrorUDOCA->GetXaxis()->SetTitle("DOCA [mm]");
   _hDocaErrorUDOCA->SetTitle(title);
   _hDocaErrorUDOCA->Draw("SAME");
   

      auto legendDOCAErrorP = new TLegend(0.4,0.64, 0.95,0.9);
      legendDOCAErrorP->SetTextSize(0.035);
      snprintf(title, sizeof(title), "RDrift - UDOCA: RDrift");
      legendDOCAErrorP->AddEntry(_hDocaErrorRDrift,title,"l");
      snprintf(title, sizeof(title), "CDrift - UDOCA: CDrift");
      legendDOCAErrorP->AddEntry(_hDocaErrorCDrift,title,"l");
      snprintf(title, sizeof(title), "RDrift - UDOCA: UDOCA");
      legendDOCAErrorP->AddEntry(_hDocaErrorUDOCA,title,"l");

      legendDOCAErrorP->SetFillStyle(0);
      legendDOCAErrorP->Draw("same");

   DocaErrorProfile->cd(1);
   snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/DocaErrorP_%s_%s.pdf", name, name2);
   DocaErrorProfile->SaveAs(title);

   // Global track/hit diagnostics that are useful for checking the selection.
   struct Hist1DSpec {
      const char* rootName;
      const char* axisTitle;
      const char* plotTitle;
   };

   const Hist1DSpec global1D[] = {
      {"a1/hNTracks",          "Number of Tracks",                 "Tracks per Event"},
      {"a1/hnDOF",             "Track Fit NDOF",                   "Track Fit NDOF"},
      {"a1/ht0",               "Track Time [ns]",                  "Track Time"},
      {"a1/hp",                "Momentum [MeV/c]",                 "Track Momentum"},
      {"a1/hpErr",             "Momentum Error [MeV/c]",           "Track Momentum Error"},
      {"a1/_hNHits",           "Number of Hits",                   "Track Hit Count"},
      {"a1/_hNHitsQuality",    "Number of Quality Hits",           "Quality Hit Count"},
      {"a1/_hNPanels",         "Number of Panels",                 "Panels Intersected"},
      {"a1/_hNPanelsQuality",  "Number of Quality Panels",         "Quality Panels Intersected"},
      {"a1/_hchi2_dof",        "#chi^{2}/NDOF",                    "Track #chi^{2}/NDOF"},
      {"a1/_hqualityScore",    "Track Quality Score",              "Track Quality Score"},
      {"a1/_hmomX",            "p_{x}/p",                          "Track Momentum X Direction"},
      {"a1/_hmomY",            "p_{y}/p",                          "Track Momentum Y Direction"},
      {"a1/_hmomZ",            "p_{z}/p",                          "Track Momentum Z Direction"},
      {"a1/_hFullSpan",        "Full Hit Span [mm]",               "Full Hit Span"},
      {"a1/hXYSpan",           "Maximum XY Span [mm]",             "Quality-Hit XY Span"},
      {"a1/hZSpan",            "Z Span [mm]",                      "Quality-Hit Z Span"},
      {"a1/hdriftTime",        "Drift Time [ns]",                  "Drift Time"},
      {"a1/_hUDOCA",           "Track UDOCA [mm]",                 "Track UDOCA"},
      {"a1/_hRDrift",          "RDrift [mm]",                      "RDrift"},
      {"a1/_hCDrift",          "CDrift [mm]",                      "CDrift"},
      {"a1/hDOCAError",        "DOCA Error [um]",                  "DOCA Error"},
      {"a1/_hXPOCAError",      "Global X POCA Error [um]",         "Global X POCA Error"},
      {"a1/_hYPOCAError",      "Global Y POCA Error [um]",         "Global Y POCA Error"},
      {"a1/_hZPOCAError",      "Z POCA Error [um]",                "Z POCA Error"},
      {"a1/_hXPOCAErrorRatio", "X POCA Error / DOCA Error",        "X POCA / DOCA Error Ratio"},
      {"a1/_hYPOCAErrorRatio", "Y POCA Error / DOCA Error",        "Y POCA / DOCA Error Ratio"},
      {"a1/_hZPOCAErrorRatio", "Z POCA Error / DOCA Error",        "Z POCA / DOCA Error Ratio"},
      {"a1/_hTrackSigma", "Track Sigma [um]",        "Track Sigma [um"}


   };

   const Int_t nGlobal1D = sizeof(global1D)/sizeof(global1D[0]);
   const Int_t nPerCanvas = 12;
   for (Int_t first=0; first<nGlobal1D; first+=nPerCanvas){
      const Int_t page = first/nPerCanvas;
      snprintf(title, sizeof(title), "cGlobalDiagnostics_%02i", page);
      TCanvas *cGlobalDiagnostics = new TCanvas(title, title, 0, 0, 1350, 1500);
      cGlobalDiagnostics->Divide(3, 4);

      for (Int_t ipad=0; ipad<nPerCanvas && first+ipad<nGlobal1D; ipad++){
         cGlobalDiagnostics->cd(ipad+1);
         TH1 *hist = nullptr;
         f1->GetObject(global1D[first+ipad].rootName, hist);
         if (!hist) continue;
         hist->SetLineColor(kBlack);
         hist->GetXaxis()->SetTitle(global1D[first+ipad].axisTitle);
         hist->SetTitle(global1D[first+ipad].plotTitle);
         hist->Draw();
      }

      snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/GlobalDiagnostics_%02i_%s_%s.pdf", page, name, name2);
      cGlobalDiagnostics->SaveAs(title);
   }

   // Panel efficiency: number of geometrically intersected panels, the subset
   // with >= 1 track hit, and their ratio.
   TH1 *hPanelIntersected = nullptr;
   TH1 *hPanelIntersectedWithHit = nullptr;
   f1->GetObject("a1/_hPanelIntersected", hPanelIntersected);
   f1->GetObject("a1/_hPanelIntersectedWithHit", hPanelIntersectedWithHit);

   if (hPanelIntersected && hPanelIntersectedWithHit){
      TCanvas *cPanelEfficiency = new TCanvas("cPanelEfficiency", "cPanelEfficiency", 0, 0, 1350, 1500);
      cPanelEfficiency->Divide(1, 3);

      cPanelEfficiency->cd(1);
      hPanelIntersected->SetLineColor(kBlack);
      hPanelIntersected->GetXaxis()->SetTitle("Panel ID");
      hPanelIntersected->GetYaxis()->SetTitle("Tracks");
      hPanelIntersected->SetTitle("Geometrically Intersected Panels");
      hPanelIntersected->Draw("hist");

      cPanelEfficiency->cd(2);
      hPanelIntersectedWithHit->SetLineColor(kBlack);
      hPanelIntersectedWithHit->GetXaxis()->SetTitle("Panel ID");
      hPanelIntersectedWithHit->GetYaxis()->SetTitle("Tracks");
      hPanelIntersectedWithHit->SetTitle("Intersected Panels with Track Hit");
      hPanelIntersectedWithHit->Draw("hist");

      // Do not clone either TH1I for the ratio: use a floating-point histogram
      // so efficiencies between 0 and 1 are retained.
      TH1F *hPanelEfficiency = new TH1F("hPanelEfficiency",
                                        "Panel Efficiency;Panel ID;Efficiency",
                                        216, -0.5, 215.5);
      hPanelEfficiency->SetDirectory(0);
      hPanelEfficiency->Divide(hPanelIntersectedWithHit, hPanelIntersected,
                               1.0, 1.0, "B");

      cPanelEfficiency->cd(3);
      hPanelEfficiency->SetLineColor(kBlack);
      hPanelEfficiency->SetMarkerColor(kBlack);
      hPanelEfficiency->SetMarkerStyle(20);
      hPanelEfficiency->SetMarkerSize(0.6);
      hPanelEfficiency->GetYaxis()->SetRangeUser(0.0, 1.05);
      hPanelEfficiency->Draw("E1");

      // Write panels with efficiency < 0.8 to CSV.
      // Panels with zero geometric intersections are skipped because their
      // efficiency is undefined rather than zero.
      snprintf(title, sizeof(title),
               "%s/EfficiencyPanels_%s_%s.csv",
               csvFolder, name, name2);
      std::ofstream lowEfficiencyFile(title);
      lowEfficiencyFile << "PanelID,StationID,PanelInStationID,MNID,Efficiency" << std::endl;

      for (Int_t panelID = 0; panelID < nPanelHistograms; panelID++){
         const Int_t bin = hPanelEfficiency->FindBin(panelID);
         const Double_t nIntersected = hPanelIntersected->GetBinContent(bin);

         if (nIntersected <= 0.0) continue;

         const Double_t efficiency = hPanelEfficiency->GetBinContent(bin);
         if (efficiency < 1.2){
            const Int_t stationID = panelID / nPanels;
            const Int_t panelInStationID = panelID % nPanels;
            const Int_t mnID = panelMNID[stationID][panelInStationID];

            lowEfficiencyFile << panelID << ","
                              << stationID << ","
                              << panelInStationID << ","
                              << mnID << ","
                              << efficiency << std::endl;
         }
      }
      lowEfficiencyFile.close();

      snprintf(title, sizeof(title),
               "/pnfs/mu2e/scratch/users/dpalo/pdfs/PanelEfficiency_%s_%s.pdf",
               name, name2);
      cPanelEfficiency->SaveAs(title);
   } else {
      std::cerr << "WARNING: Could not find a1/_hPanelIntersected and/or "
                << "a1/_hPanelIntersectedWithHit" << std::endl;
   }

   // Global 2D POCA diagnostics.
   TCanvas *cPOCAGlobal2D = new TCanvas("cPOCAGlobal2D", "cPOCAGlobal2D", 0, 0, 1350, 700);
   cPOCAGlobal2D->Divide(2, 1);

   TH2 *hXYPOCA = nullptr;
   cPOCAGlobal2D->cd(1);
   f1->GetObject("a1/hXYPOCA", hXYPOCA);
   if (hXYPOCA){
      hXYPOCA->GetXaxis()->SetTitle("Global X POCA [mm]");
      hXYPOCA->GetYaxis()->SetTitle("Global Y POCA [mm]");
      hXYPOCA->Draw("colz");
   }

   TH2 *hYZPOCA = nullptr;
   cPOCAGlobal2D->cd(2);
   f1->GetObject("a1/hYZPOCA", hYZPOCA);
   if (hYZPOCA){
      hYZPOCA->GetXaxis()->SetTitle("Global Y POCA [mm]");
      hYZPOCA->GetYaxis()->SetTitle("Global Z POCA [mm]");
      hYZPOCA->Draw("colz");
   }

   snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/POCAGlobal2D_%s_%s.pdf", name, name2);
   cPOCAGlobal2D->SaveAs(title);

   // Panel-level wire occupancy with bad-wire overlays.
   TH1F *hWire[nPanelHistograms];
   TH1F *hWireBad[nPanelHistograms];
   TH1F *hWireBadShort[nPanelHistograms];

   for (Int_t istation = 0; istation < nStations; istation++){
      snprintf(title, sizeof(title), "cWire_Station%02i", istation);
      TCanvas *cWire = new TCanvas(title, title, 0, 0, 1350, 1500);
      cWire->Divide(3, 4);

      for (Int_t ihist = 0; ihist < nPanels; ihist++){
         const Int_t globalHist = istation * nPanels + ihist;
         cWire->cd(ihist + 1);

         snprintf(title, sizeof(title), "a1/hWire%i", globalHist);
         f1->GetObject(title, hWire[globalHist]);
         if (!hWire[globalHist]) continue;

         hWire[globalHist]->SetLineColor(kBlack);
         hWire[globalHist]->GetXaxis()->SetTitle("Wire Number");
         //hWire[globalHist]->GetYaxis()->SetRangeUser(0, 100);
         snprintf(title, sizeof(title), "S%i, Wire Occ Pnl %i, MN %i", istation, ihist, panelMNID[istation][ihist]);
         hWire[globalHist]->SetTitle(title);
         hWire[globalHist]->Draw("hist");

         snprintf(title, sizeof(title), "hWireBad_Station%02i_Panel%02i", istation, ihist);
         hWireBad[globalHist] = (TH1F*)hWire[globalHist]->Clone(title);
         hWireBad[globalHist]->Reset("ICES");
         hWireBad[globalHist]->SetDirectory(0);
         hWireBad[globalHist]->SetLineColor(kRed);
         hWireBad[globalHist]->SetFillColorAlpha(kRed, 0.35);

         snprintf(title, sizeof(title), "hWireBadShort_Station%02i_Panel%02i", istation, ihist);
         hWireBadShort[globalHist] = (TH1F*)hWire[globalHist]->Clone(title);
         hWireBadShort[globalHist]->Reset("ICES");
         hWireBadShort[globalHist]->SetDirectory(0);
         hWireBadShort[globalHist]->SetLineColor(kBlue);
         hWireBadShort[globalHist]->SetFillColorAlpha(kBlue, 0.35);

         for (Int_t iwire = 0; iwire < nWiresPerPanel; iwire++){
            if (badWire[istation][ihist][iwire]){
               const Int_t ibin = hWireBad[globalHist]->FindBin(iwire - 1);
               if (ibin >= 1 && ibin <= hWireBad[globalHist]->GetNbinsX()){
                  hWireBad[globalHist]->SetBinContent(ibin, 50.0);
               }
            }

            if (badWireShort[istation][ihist][iwire]){
               const Int_t ibin = hWireBadShort[globalHist]->FindBin(iwire - 1);
               if (ibin >= 1 && ibin <= hWireBadShort[globalHist]->GetNbinsX()){
                  hWireBadShort[globalHist]->SetBinContent(ibin, 50.0);
               }
            }
         }

         hWireBad[globalHist]->Draw("same");
         hWireBadShort[globalHist]->Draw("same");
      }

      snprintf(title, sizeof(title),
               "/pnfs/mu2e/scratch/users/dpalo/pdfs/WireOccupancy_Station%02i_%s_%s.pdf",
               istation, name, name2);
      cWire->SaveAs(title);
   }

   TH1F *hZPOCAErrorPanel[nPanelHistograms];
   fDOCA->SetParLimits(5,0.05, 0.3);

   for (Int_t istation = 0; istation < nStations; istation++){
      snprintf(title, sizeof(title), "cZPOCAErrorPanel_Station%02i", istation);
      TCanvas *cZPOCAErrorPanel = new TCanvas(title, title, 0, 0, 1350, 1500);
      cZPOCAErrorPanel->Divide(3, 4);

      for (Int_t ihist = 0; ihist < nPanels; ihist++){
         const Int_t globalHist = istation * nPanels + ihist;
         cZPOCAErrorPanel->cd(ihist+1);
         snprintf(title, sizeof(title), "a1/_hZPOCAErrorPanel%i", globalHist);
         f1->GetObject(title, hZPOCAErrorPanel[globalHist]);
         if (!hZPOCAErrorPanel[globalHist]) continue;
         hZPOCAErrorPanel[globalHist]->SetLineColor(kBlack);
         hZPOCAErrorPanel[globalHist]->GetXaxis()->SetTitle("ZPOCA Error [um]");
         snprintf(title, sizeof(title), "S%i, ZPOCA Err Pnl %i, MN %i", istation, ihist, panelMNID[istation][ihist]);
         hZPOCAErrorPanel[globalHist]->SetTitle(title);
         hZPOCAErrorPanel[globalHist]->Draw();
         const Double_t entriesZPOCA = hZPOCAErrorPanel[globalHist]->GetEntries();
         if (entriesZPOCA >= 500){
         hZPOCAErrorPanel[globalHist]->Fit(fDOCA,"", "", -1500, 1500);
   Double_t meanCore2 = fDOCA->GetParameter(1);
            Double_t sigmaCore2 = fDOCA->GetParameter(2);
            Double_t fractionTail2 = fDOCA->GetParameter(5);
            Double_t meanTail2 = fDOCA->GetParameter(3);
            Double_t sigmaTail2 = fDOCA->GetParameter(4);
            auto legendZPOCAPanel = new TLegend(0.55,0.64, 0.95,0.9);
            legendZPOCAPanel->SetTextSize(0.045);
            snprintf(title, sizeof(title), "#splitline{#sigma_{1}=%.2f, #mu_{1}=%.1f [um]}{#splitline{#sigma_{2}=%.2f, #mu_{2}=%.1f [um]}{A_{2}/A_{1}=%.2f}}", sigmaCore2,  meanCore2, sigmaTail2, meanTail2, fractionTail2);
            legendZPOCAPanel->AddEntry(hZPOCAErrorPanel[globalHist],title,"l");
            legendZPOCAPanel->SetFillStyle(0);
            legendZPOCAPanel->Draw("same");
         } else {
            std::cout << "Skipping ZPOCA panel fit for global panel " << globalHist
                      << ": " << entriesZPOCA << " entries." << std::endl;
         }
      }
      cZPOCAErrorPanel->cd(0);
      snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/ZPOCAerrorPanel_Station%02i_%s_%s.pdf", istation, name, name2);
      cZPOCAErrorPanel->SaveAs(title);
   }

   
    // X POCA residual per panel.
   TH1F *hXPOCAErrorPanel[nPanelHistograms];
   fDOCA->SetParLimits(5,0.05, 0.3);

   for (Int_t istation = 0; istation < nStations; istation++){
      snprintf(title, sizeof(title), "cXPOCAErrorPanel_Station%02i", istation);
      TCanvas *cXPOCAErrorPanel = new TCanvas(title, title, 0, 0, 1350, 1500);
      cXPOCAErrorPanel->Divide(3, 4);

      for (Int_t ihist = 0; ihist < nPanels; ihist++){
         const Int_t globalHist = istation * nPanels + ihist;
         cXPOCAErrorPanel->cd(ihist+1);
         snprintf(title, sizeof(title), "a1/_hXPOCAErrorPanel%i", globalHist);
         f1->GetObject(title, hXPOCAErrorPanel[globalHist]);
         if (!hXPOCAErrorPanel[globalHist]) continue;
         hXPOCAErrorPanel[globalHist]->SetLineColor(kBlack);
         hXPOCAErrorPanel[globalHist]->GetXaxis()->SetTitle("X POCA Error [um]");
         snprintf(title, sizeof(title), "S%i, XPOCA Err Pnl %i, MN %i", istation, ihist, panelMNID[istation][ihist]);
         hXPOCAErrorPanel[globalHist]->SetTitle(title);
         hXPOCAErrorPanel[globalHist]->Draw();
         const Double_t entriesXPOCA = hXPOCAErrorPanel[globalHist]->GetEntries();
         if (entriesXPOCA >= 500){
         hXPOCAErrorPanel[globalHist]->Fit(fDOCA,"", "", -1500, 1500);
   Double_t meanCore2 = fDOCA->GetParameter(1);
            Double_t sigmaCore2 = fDOCA->GetParameter(2);
            Double_t fractionTail2 = fDOCA->GetParameter(5);
            Double_t meanTail2 = fDOCA->GetParameter(3);
            Double_t sigmaTail2 = fDOCA->GetParameter(4);
            auto legendXPOCAPanel = new TLegend(0.55,0.64, 0.95,0.9);
            legendXPOCAPanel->SetTextSize(0.045);
            snprintf(title, sizeof(title), "#splitline{#sigma_{1}=%.2f, #mu_{1}=%.1f [um]}{#splitline{#sigma_{2}=%.2f, #mu_{2}=%.1f [um]}{A_{2}/A_{1}=%.2f}}", sigmaCore2, meanCore2, sigmaTail2, meanTail2, fractionTail2);
            legendXPOCAPanel->AddEntry(hXPOCAErrorPanel[globalHist],title,"l");
            legendXPOCAPanel->SetFillStyle(0);
            legendXPOCAPanel->Draw("same");
         } else {
            std::cout << "Skipping XPOCA panel fit for global panel " << globalHist
                      << ": " << entriesXPOCA << " entries." << std::endl;
         }
      }

      snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/XPOCAerrorPanel_Station%02i_%s_%s.pdf", istation, name, name2);
      cXPOCAErrorPanel->SaveAs(title);
   }

   TH1F *hYPOCAErrorPanel[nPanelHistograms];
   fDOCA->SetParLimits(5,0.05, 0.3);

   for (Int_t istation = 0; istation < nStations; istation++){
      snprintf(title, sizeof(title), "cYPOCAErrorPanel_Station%02i", istation);
      TCanvas *cYPOCAErrorPanel = new TCanvas(title, title, 0, 0, 1350, 1500);
      cYPOCAErrorPanel->Divide(3, 4);

      for (Int_t ihist = 0; ihist < nPanels; ihist++){
         const Int_t globalHist = istation * nPanels + ihist;
         cYPOCAErrorPanel->cd(ihist+1);
         snprintf(title, sizeof(title), "a1/_hYPOCAErrorPanel%i", globalHist);
         f1->GetObject(title, hYPOCAErrorPanel[globalHist]);
         if (!hYPOCAErrorPanel[globalHist]) continue;
         hYPOCAErrorPanel[globalHist]->SetLineColor(kBlack);
         hYPOCAErrorPanel[globalHist]->GetXaxis()->SetTitle("YPOCA Error [um]");
         snprintf(title, sizeof(title), "S%i, YPOCA Err Pnl %i, MN %i", istation, ihist, panelMNID[istation][ihist]);
         hYPOCAErrorPanel[globalHist]->SetTitle(title);
         hYPOCAErrorPanel[globalHist]->Draw();
         const Double_t entriesYPOCA = hYPOCAErrorPanel[globalHist]->GetEntries();
         if (entriesYPOCA >= 500){
         hYPOCAErrorPanel[globalHist]->Fit(fDOCA,"", "", -1500, 1500);
   Double_t meanCore2 = fDOCA->GetParameter(1);
            Double_t sigmaCore2 = fDOCA->GetParameter(2);
            Double_t fractionTail2 = fDOCA->GetParameter(5);
            Double_t meanTail2 = fDOCA->GetParameter(3);
            Double_t sigmaTail2 = fDOCA->GetParameter(4);
            auto legendYPOCAPanel = new TLegend(0.55,0.64, 0.95,0.9);
            legendYPOCAPanel->SetTextSize(0.045);
            snprintf(title, sizeof(title), "#splitline{#sigma_{1}=%.2f, #mu_{1}=%.1f [um]}{#splitline{#sigma_{2}=%.2f, #mu_{2}=%.1f [um]}{A_{2}/A_{1}=%.2f}}", sigmaCore2,  meanCore2, sigmaTail2, meanTail2, fractionTail2);
            legendYPOCAPanel->AddEntry(hYPOCAErrorPanel[globalHist],title,"l");
            legendYPOCAPanel->SetFillStyle(0);
            legendYPOCAPanel->Draw("same");
         } else {
            std::cout << "Skipping YPOCA panel fit for global panel " << globalHist
                      << ": " << entriesYPOCA << " entries." << std::endl;
         }
      }
      cYPOCAErrorPanel->cd(0);
      snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/YPOCAerrorPanel_Station%02i_%s_%s.pdf", istation, name, name2);
      cYPOCAErrorPanel->SaveAs(title);
   }

   

   
   TProfile *hZPOCAErrorWirePanelEven[nPanelHistograms];
   TProfile *hZPOCAErrorWirePanelOdd[nPanelHistograms];

   for (Int_t istation = 0; istation < nStations; istation++){
      snprintf(title, sizeof(title), "cZPOCAErrorWirePanel_Station%02i", istation);
      TCanvas *cZPOCAErrorWirePanel = new TCanvas(title, title, 0, 0, 1350, 1500);
      cZPOCAErrorWirePanel->Divide(3, 4);

      for (Int_t ihist = 0; ihist < nPanels; ihist++){
         const Int_t globalHist = istation * nPanels + ihist;
         cZPOCAErrorWirePanel->cd(ihist+1);
         snprintf(title, sizeof(title), "a1/hZPOCAErrorWireEven%i", globalHist);
         f1->GetObject(title, hZPOCAErrorWirePanelEven[globalHist]);
         if (!hZPOCAErrorWirePanelEven[globalHist]) continue;
         hZPOCAErrorWirePanelEven[globalHist]->SetLineColor(kBlack);
         hZPOCAErrorWirePanelEven[globalHist]->GetXaxis()->SetTitle("Wire Number");
         hZPOCAErrorWirePanelEven[globalHist]->GetYaxis()->SetTitle("Z POCA Error [um]");
         hZPOCAErrorWirePanelEven[globalHist]->GetYaxis()->SetRangeUser(-150, 150);
         snprintf(title, sizeof(title), "S%i, ZPOCA/Wire Pnl %i, MN %i", istation, ihist, panelMNID[istation][ihist]);
         hZPOCAErrorWirePanelEven[globalHist]->SetTitle(title);
         hZPOCAErrorWirePanelEven[globalHist]->Draw("same");

         snprintf(title, sizeof(title), "a1/hZPOCAErrorWireOdd%i", globalHist);
         f1->GetObject(title, hZPOCAErrorWirePanelOdd[globalHist]);
         if (!hZPOCAErrorWirePanelOdd[globalHist]) continue;
         hZPOCAErrorWirePanelOdd[globalHist]->SetLineColor(kRed);
         hZPOCAErrorWirePanelOdd[globalHist]->GetXaxis()->SetTitle("Wire Number");
         hZPOCAErrorWirePanelOdd[globalHist]->GetYaxis()->SetTitle("Z POCA Error [um]");
         hZPOCAErrorWirePanelOdd[globalHist]->GetYaxis()->SetRangeUser(-150, 150);
         hZPOCAErrorWirePanelOdd[globalHist]->SetTitle(title);
         hZPOCAErrorWirePanelOdd[globalHist]->Draw("same");
      }
      cZPOCAErrorWirePanel->cd(0);
      snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/ZPOCAErrorWirePanel_Station%02i_%s_%s.pdf", istation, name, name2);
      cZPOCAErrorWirePanel->SaveAs(title);
   }
   

   TProfile *hVPOCAErrorWirePanelEven[nPanelHistograms];
   TProfile *hVPOCAErrorWirePanelOdd[nPanelHistograms];

   for (Int_t istation = 0; istation < nStations; istation++){
      snprintf(title, sizeof(title), "cVPOCAErrorWirePanel_Station%02i", istation);
      TCanvas *cVPOCAErrorWirePanel = new TCanvas(title, title, 0, 0, 1350, 1500);
      cVPOCAErrorWirePanel->Divide(3, 4);

      for (Int_t ihist = 0; ihist < nPanels; ihist++){
         const Int_t globalHist = istation * nPanels + ihist;
         cVPOCAErrorWirePanel->cd(ihist+1);
         snprintf(title, sizeof(title), "a1/hVPOCAErrorWireEven%i", globalHist);
         f1->GetObject(title, hVPOCAErrorWirePanelEven[globalHist]);
         if (!hVPOCAErrorWirePanelEven[globalHist]) continue;
         hVPOCAErrorWirePanelEven[globalHist]->SetLineColor(kBlack);
         hVPOCAErrorWirePanelEven[globalHist]->GetXaxis()->SetTitle("Wire Number");
         hVPOCAErrorWirePanelEven[globalHist]->GetYaxis()->SetTitle("V POCA Error [um]");
         hVPOCAErrorWirePanelEven[globalHist]->GetYaxis()->SetRangeUser(-400, 400);
         snprintf(title, sizeof(title), "S%i, VPOCA/Wire Pnl %i, MN %i", istation, ihist, panelMNID[istation][ihist]);
         hVPOCAErrorWirePanelEven[globalHist]->SetTitle(title);
         hVPOCAErrorWirePanelEven[globalHist]->Draw("same");

         snprintf(title, sizeof(title), "a1/hVPOCAErrorWireOdd%i", globalHist);
         f1->GetObject(title, hVPOCAErrorWirePanelOdd[globalHist]);
         if (!hVPOCAErrorWirePanelOdd[globalHist]) continue;
         hVPOCAErrorWirePanelOdd[globalHist]->SetLineColor(kRed);
         hVPOCAErrorWirePanelOdd[globalHist]->GetXaxis()->SetTitle("Wire Number");
         hVPOCAErrorWirePanelOdd[globalHist]->GetYaxis()->SetTitle("V POCA Error [um]");
         hVPOCAErrorWirePanelOdd[globalHist]->GetYaxis()->SetRangeUser(-400, 1400);
         hVPOCAErrorWirePanelOdd[globalHist]->SetTitle(title);
         hVPOCAErrorWirePanelOdd[globalHist]->Draw("same");
      }
      cVPOCAErrorWirePanel->cd(0);
      snprintf(title, sizeof(title), "/pnfs/mu2e/scratch/users/dpalo/pdfs/VPOCAErrorWirePanel_Station%02i_%s_%s.pdf", istation, name, name2);
      cVPOCAErrorWirePanel->SaveAs(title);
   }

   // The following histograms are one-per-wire (~20k total) and are expensive to draw.
   if (drawWireHistograms){
      TCanvas *VErrorZWireCanvas[nStations][nPanels][nSectorsPerPanel];
      TCanvas *WErrorZWireCanvas[nStations][nPanels][nSectorsPerPanel];

      TH2F *hVErrorZWire2D[nWireHistograms] = {};
      TH2F *hWErrorZWire2D[nWireHistograms] = {};
      TProfile *hVErrorZWireP[nWireHistograms] = {};
      TProfile *hWErrorZWireP[nWireHistograms] = {};

      float intercepts[nWireHistograms] = {};
      float slopes[nWireHistograms] = {};
      float interceptsError[nWireHistograms] = {};
      float slopesError[nWireHistograms] = {};

      for (Int_t istation = 0; istation < nStations; istation++){
         for (Int_t ipanel = 0; ipanel < nPanels; ipanel++){
            for (Int_t isector = 0; isector < nSectorsPerPanel; isector++){
               snprintf(title, sizeof(title),
                        "VErrorZWireCanvas_Station%02i_Panel%02i_Sector%02i",
                        istation, ipanel, isector);
               VErrorZWireCanvas[istation][ipanel][isector] =
                   new TCanvas(title, title, 0, 0, 1000, 1000);
               VErrorZWireCanvas[istation][ipanel][isector]->Divide(3, 4);

               snprintf(title, sizeof(title),
                        "WErrorZWireCanvas_Station%02i_Panel%02i_Sector%02i",
                        istation, ipanel, isector);
               WErrorZWireCanvas[istation][ipanel][isector] =
                   new TCanvas(title, title, 0, 0, 1000, 1000);
               WErrorZWireCanvas[istation][ipanel][isector]->Divide(3, 4);
            }
         }
      }

      for (Int_t istation = 0; istation < nStations; istation++){
         for (Int_t ipanel = 0; ipanel < nPanels; ipanel++){
            for (Int_t isector = 0; isector < nSectorsPerPanel; isector++){
               for (Int_t iwire = 0; iwire < nWiresPerSector; iwire++){
                  const Int_t wireIndex = istation*nPanels*nWiresPerPanel
                                        + ipanel*nWiresPerPanel
                                        + isector*nWiresPerSector
                                        + iwire;

                  VErrorZWireCanvas[istation][ipanel][isector]->cd(iwire+1);

                  snprintf(title, sizeof(title), "a1/hVErrorZWire2D%i", wireIndex);
                  f1->GetObject(title, hVErrorZWire2D[wireIndex]);
                  snprintf(title, sizeof(title), "a1/hVErrorZWireP%i", wireIndex);
                  f1->GetObject(title, hVErrorZWireP[wireIndex]);
                  if (!hVErrorZWire2D[wireIndex] || !hVErrorZWireP[wireIndex]) continue;

                  hVErrorZWire2D[wireIndex]->GetYaxis()->SetRangeUser(-150, 150);
                  hVErrorZWire2D[wireIndex]->Draw("colz");
                  hVErrorZWireP[wireIndex]->SetStats(0);
                  hVErrorZWireP[wireIndex]->Draw("same");

                  snprintf(title, sizeof(title), "fVWire_%i", wireIndex);
                  TF1 *fVWire = new TF1(title, "[0] + [1]*x", -250, 250);
                  fVWire->SetParameters(0, 0);
                  hVErrorZWireP[wireIndex]->Fit(fVWire, "Q", "", -250, 250);
                  fVWire->Draw("same");
               }

               snprintf(title, sizeof(title),
                        "pdfs/VErrorZWire_Station%02i_Panel%02i_Sector%02i_%s_%s.pdf",
                        istation, ipanel, isector, name, name2);
               VErrorZWireCanvas[istation][ipanel][isector]->SaveAs(title);
            }
         }
      }

      TProfile::Approximate();
      for (Int_t istation = 0; istation < nStations; istation++){
         for (Int_t ipanel = 0; ipanel < nPanels; ipanel++){
            for (Int_t isector = 0; isector < nSectorsPerPanel; isector++){
               for (Int_t iwire = 0; iwire < nWiresPerSector; iwire++){
                  const Int_t wireIndex = istation*nPanels*nWiresPerPanel
                                        + ipanel*nWiresPerPanel
                                        + isector*nWiresPerSector
                                        + iwire;

                  WErrorZWireCanvas[istation][ipanel][isector]->cd(iwire+1);

                  snprintf(title, sizeof(title), "a1/hWErrorZWire2D%i", wireIndex);
                  f1->GetObject(title, hWErrorZWire2D[wireIndex]);
                  snprintf(title, sizeof(title), "a1/hWErrorZWireP%i", wireIndex);
                  f1->GetObject(title, hWErrorZWireP[wireIndex]);
                  if (!hWErrorZWire2D[wireIndex] || !hWErrorZWireP[wireIndex]) continue;

                  hWErrorZWireP[wireIndex]->Approximate();
                  hWErrorZWire2D[wireIndex]->GetYaxis()->SetRangeUser(-300, 300);
                  hWErrorZWireP[wireIndex]->GetYaxis()->SetRangeUser(-300, 300);
                  hWErrorZWire2D[wireIndex]->Draw("colz");
                  hWErrorZWireP[wireIndex]->SetStats(0);
                  hWErrorZWireP[wireIndex]->Draw("same");

                  snprintf(title, sizeof(title), "fWWire_%i", wireIndex);
                  TF1 *fWWire = new TF1(title, "[0] + [1]*x", -500, 500);
                  fWWire->SetParameters(0, 0);
                  hWErrorZWireP[wireIndex]->Fit(fWWire, "Q", "", -500, 500);
                  fWWire->Draw("same");

                  if (hWErrorZWireP[wireIndex]->GetEntries() > 400){
                     intercepts[wireIndex] = fWWire->GetParameter(0);
                     slopes[wireIndex] = fWWire->GetParameter(1);
                     interceptsError[wireIndex] = fWWire->GetParError(0);
                     slopesError[wireIndex] = fWWire->GetParError(1);
                  }
               }

               snprintf(title, sizeof(title),
                        "pdfs/WErrorZWire_Station%02i_Panel%02i_Sector%02i_%s_%s.pdf",
                        istation, ipanel, isector, name, name2);
               WErrorZWireCanvas[istation][ipanel][isector]->SaveAs(title);
            }
         }
      }

      std::ofstream myfile(csvName);
      for (Int_t ibin = 0; ibin < nWireHistograms; ibin++){
         myfile << ibin << ", "
                << intercepts[ibin] << ", "
                << slopes[ibin] << ", "
                << interceptsError[ibin] << ", "
                << slopesError[ibin] << std::endl;
      }
   }

   

}
