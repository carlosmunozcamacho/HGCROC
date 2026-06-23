#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>

#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TLine.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

void analyze_LED(int runNumber = 691,
                 const char *inputDir = ".",
                 const char *outputDir = "analysis_output",
                 int integralStartSample = 2,
                 int integralEndSample = 8,
                 int integralBins = 240,
                 int totBins = 120)
{
  const int nChannels = 144;
  const int nSamples = 20;
  const int channels[16] = {
      36, 37, 38, 39,
      40, 41, 42, 43,
      45, 46, 47, 48,
      49, 50, 51, 52};

  const TString inputFileName = TString::Format("%s/Run%d.root", inputDir, runNumber);
  TFile *inputFile = TFile::Open(inputFileName, "READ");
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Could not open " << inputFileName << std::endl;
    return;
  }

  TTree *tree = nullptr;
  inputFile->GetObject("events", tree);
  if (!tree) {
    std::cerr << "Could not find TTree named 'events' in " << inputFileName << std::endl;
    inputFile->ls();
    inputFile->Close();
    return;
  }

  const TString runOutputDir = TString::Format("%s/Run%d", outputDir, runNumber);
  gSystem->mkdir(runOutputDir, kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  TH1::AddDirectory(kFALSE);

  UInt_t adc[nChannels][nSamples] = {};
  UInt_t tot[nChannels][nSamples] = {};
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("adc", 1);
  tree->SetBranchStatus("tot", 1);
  tree->SetBranchAddress("adc", adc);
  tree->SetBranchAddress("tot", tot);

  integralStartSample = std::max(0, std::min(integralStartSample, nSamples - 1));
  integralEndSample = std::max(0, std::min(integralEndSample, nSamples - 1));
  if (integralStartSample > integralEndSample) {
    std::swap(integralStartSample, integralEndSample);
  }
  integralBins = std::max(20, integralBins);
  totBins = std::max(20, totBins);

  std::vector<std::vector<double>> integralValues(16);
  for (auto &values : integralValues) {
    values.reserve(tree->GetEntries());
  }
  std::vector<std::vector<double>> totIntegralValues(16);
  for (auto &values : totIntegralValues) {
    values.reserve(tree->GetEntries());
  }

  double globalMaximum = 0.0;
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
    tree->GetEntry(entry);
    for (int channelIndex = 0; channelIndex < 16; ++channelIndex) {
      const int channel = channels[channelIndex];
      double integral = 0.0;
      double totIntegral = 0.0;
      for (int sample = 0; sample < nSamples; ++sample) {
        globalMaximum = TMath::Max(globalMaximum, static_cast<double>(adc[channel][sample]));
        if (sample >= integralStartSample && sample <= integralEndSample) {
          integral += adc[channel][sample];
        }
        totIntegral += tot[channel][sample];
      }
      integralValues[channelIndex].push_back(integral);
      totIntegralValues[channelIndex].push_back(totIntegral);
    }
  }
  const double yMaximum = 1.05 * globalMaximum;

  TCanvas *cWaveforms = new TCanvas("cWaveforms", "Selected channel waveforms", 1600, 1200);
  cWaveforms->Divide(4, 4, 0.001, 0.001);

  std::vector<TH2D *> waveformHistograms;
  waveformHistograms.reserve(16);
  std::vector<TLine *> integrationWindowLines;
  integrationWindowLines.reserve(32);

  const double integrationWindowLowEdge = integralStartSample - 0.5;
  const double integrationWindowHighEdge = integralEndSample + 0.5;

  for (int i = 0; i < 16; ++i) {
    const int channel = channels[i];
    cWaveforms->cd(i + 1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->SetRightMargin(0.14);
    gPad->SetLeftMargin(0.10);
    gPad->SetBottomMargin(0.10);
    gPad->SetTopMargin(0.12);

    TH2D *histogram = new TH2D(TString::Format("hWaveform_ch%d", channel),
                               TString::Format("Run %d  Channel %d", runNumber, channel),
                               nSamples, -0.5, nSamples - 0.5,
                               120, 0.0, yMaximum);
    histogram->SetDirectory(nullptr);
    waveformHistograms.push_back(histogram);

    histogram->SetTitle(TString::Format("Run %d  Channel %d", runNumber, channel));
    histogram->GetYaxis()->SetRangeUser(0.0, yMaximum);
    histogram->GetXaxis()->SetTitle("Sample");
    histogram->GetYaxis()->SetTitle("ADC");
    for (int sample = 0; sample < nSamples; ++sample) {
      histogram->GetXaxis()->SetBinLabel(sample + 1, TString::Format("%d", sample));
    }
    histogram->GetXaxis()->LabelsOption("h");
    histogram->GetXaxis()->SetTitleSize(0.045);
    histogram->GetYaxis()->SetTitleSize(0.045);
    histogram->GetXaxis()->SetLabelSize(0.040);
    histogram->GetYaxis()->SetLabelSize(0.040);
    histogram->SetTitleSize(0.070, "t");

    for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
      tree->GetEntry(entry);
      for (int sample = 0; sample < nSamples; ++sample) {
        histogram->Fill(sample, adc[channel][sample]);
      }
    }
    histogram->Draw("colz");

    TLine *lowLine = new TLine(integrationWindowLowEdge, 0.0, integrationWindowLowEdge, yMaximum);
    TLine *highLine = new TLine(integrationWindowHighEdge, 0.0, integrationWindowHighEdge, yMaximum);
    for (TLine *line : {lowLine, highLine}) {
      line->SetLineColor(kRed + 1);
      line->SetLineStyle(2);
      line->SetLineWidth(2);
      line->Draw("same");
      integrationWindowLines.push_back(line);
    }
  }

  cWaveforms->Update();

  const TString pdfName = TString::Format("%s/Run%d_C1_selected_channel_waveforms.pdf",
                                          runOutputDir.Data(), runNumber);
  const TString pngName = TString::Format("%s/Run%d_C1_selected_channel_waveforms.png",
                                          runOutputDir.Data(), runNumber);
  cWaveforms->SaveAs(pdfName);
  cWaveforms->SaveAs(pngName);

  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(1110);

  TCanvas *cIntegrals = new TCanvas("cIntegrals", "Selected channel waveform integrals", 1600, 1200);
  cIntegrals->Divide(4, 4, 0.001, 0.001);

  std::vector<TH1D *> integralHistograms;
  integralHistograms.reserve(16);
  std::vector<TF1 *> gaussianFits;
  gaussianFits.reserve(16);
  std::vector<double> rawFitMeans(16, 0.0);

  double globalIntegralMinimum = integralValues[0][0];
  double globalIntegralMaximum = integralValues[0][0];
  for (const auto &values : integralValues) {
    const auto minmax = std::minmax_element(values.begin(), values.end());
    globalIntegralMinimum = std::min(globalIntegralMinimum, *minmax.first);
    globalIntegralMaximum = std::max(globalIntegralMaximum, *minmax.second);
  }
  const double integralPadding = std::max(10.0, 0.05 * (globalIntegralMaximum - globalIntegralMinimum));
  const double integralXMin = globalIntegralMinimum - integralPadding;
  const double integralXMax = globalIntegralMaximum + integralPadding;
  const double integralBinWidth = (integralXMax - integralXMin) / integralBins;

  for (int i = 0; i < 16; ++i) {
    const int channel = channels[i];
    TH1D *histogram = new TH1D(TString::Format("hIntegral_ch%d", channel),
                               TString::Format("Run %d  Channel %d  Integral samples %d-%d",
                                               runNumber, channel, integralStartSample, integralEndSample),
                               integralBins, integralXMin, integralXMax);
    histogram->SetDirectory(nullptr);
    histogram->SetLineColor(kBlue + 1);
    histogram->SetLineWidth(2);
    histogram->GetXaxis()->SetTitle("Integrated ADC");
    histogram->GetYaxis()->SetTitle("Entries");
    histogram->GetXaxis()->SetTitleSize(0.045);
    histogram->GetYaxis()->SetTitleSize(0.045);
    histogram->GetXaxis()->SetLabelSize(0.040);
    histogram->GetYaxis()->SetLabelSize(0.040);
    histogram->SetTitleSize(0.060, "t");

    for (double value : integralValues[i]) {
      histogram->Fill(value);
    }
    integralHistograms.push_back(histogram);
  }

  double globalIntegralYMaximum = 0.0;
  for (TH1D *histogram : integralHistograms) {
    globalIntegralYMaximum = std::max(globalIntegralYMaximum, histogram->GetMaximum());
  }
  const double integralYMaximum = 1.15 * globalIntegralYMaximum;

  for (int i = 0; i < 16; ++i) {
    const int channel = channels[i];
    TH1D *histogram = integralHistograms[i];

    cIntegrals->cd(i + 1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.11);
    gPad->SetTopMargin(0.12);

    histogram->SetMinimum(0.0);
    histogram->SetMaximum(integralYMaximum);

    const int peakBin = histogram->GetMaximumBin();
    const double peakX = histogram->GetBinCenter(peakBin);
    const double rms = histogram->GetRMS();
    const double fitMin = std::max(histogram->GetXaxis()->GetXmin(), peakX - 2.0 * rms);
    const double fitMax = std::min(histogram->GetXaxis()->GetXmax(), peakX + 2.0 * rms);

    TF1 *fit = new TF1(TString::Format("fGaussian_ch%d", channel), "gaus", fitMin, fitMax);
    fit->SetLineColor(kRed + 1);
    fit->SetLineWidth(2);
    gaussianFits.push_back(fit);

    histogram->Fit(fit, "RQ");
    rawFitMeans[i] = fit->GetParameter(1);
    histogram->Draw();
    fit->Draw("same");
    gPad->Update();
  }

  cIntegrals->Update();

  const TString integralPdfName = TString::Format("%s/Run%d_C2_integrals_samples_%d_%d.pdf",
                                                 runOutputDir.Data(), runNumber,
                                                 integralStartSample, integralEndSample);
  const TString integralPngName = TString::Format("%s/Run%d_C2_integrals_samples_%d_%d.png",
                                                 runOutputDir.Data(), runNumber,
                                                 integralStartSample, integralEndSample);
  cIntegrals->SaveAs(integralPdfName);
  cIntegrals->SaveAs(integralPngName);

  const double averageRawPeak =
      std::accumulate(rawFitMeans.begin(), rawFitMeans.end(), 0.0) / rawFitMeans.size();

  std::vector<std::vector<double>> calibratedIntegralValues(16);
  for (int i = 0; i < 16; ++i) {
    const double calibrationFactor = (rawFitMeans[i] != 0.0) ? averageRawPeak / rawFitMeans[i] : 1.0;
    calibratedIntegralValues[i].reserve(integralValues[i].size());
    for (double value : integralValues[i]) {
      calibratedIntegralValues[i].push_back(value * calibrationFactor);
    }
  }

  TCanvas *cCalibratedIntegrals =
      new TCanvas("cCalibratedIntegrals", "Calibrated selected channel waveform integrals", 1600, 1200);
  cCalibratedIntegrals->Divide(4, 4, 0.001, 0.001);

  double globalCalibratedIntegralMinimum = calibratedIntegralValues[0][0];
  double globalCalibratedIntegralMaximum = calibratedIntegralValues[0][0];
  for (const auto &values : calibratedIntegralValues) {
    const auto minmax = std::minmax_element(values.begin(), values.end());
    globalCalibratedIntegralMinimum = std::min(globalCalibratedIntegralMinimum, *minmax.first);
    globalCalibratedIntegralMaximum = std::max(globalCalibratedIntegralMaximum, *minmax.second);
  }
  const double calibratedIntegralPadding =
      std::max(10.0, 0.05 * (globalCalibratedIntegralMaximum - globalCalibratedIntegralMinimum));
  const double calibratedIntegralXMin = globalCalibratedIntegralMinimum - calibratedIntegralPadding;
  const double calibratedIntegralXMax = globalCalibratedIntegralMaximum + calibratedIntegralPadding;
  const int calibratedIntegralBins =
      std::max(20, TMath::CeilNint((calibratedIntegralXMax - calibratedIntegralXMin) / integralBinWidth));

  std::vector<TH1D *> calibratedIntegralHistograms;
  calibratedIntegralHistograms.reserve(16);
  std::vector<TF1 *> calibratedGaussianFits;
  calibratedGaussianFits.reserve(16);
  std::vector<double> calibratedFitMeans(16, 0.0);
  std::vector<double> calibratedFitMeanErrors(16, 0.0);
  std::vector<double> calibratedFitSigmas(16, 0.0);
  std::vector<double> calibratedFitSigmaErrors(16, 0.0);

  for (int i = 0; i < 16; ++i) {
    const int channel = channels[i];
    TH1D *histogram = new TH1D(TString::Format("hCalibratedIntegral_ch%d", channel),
                               TString::Format("Run %d  Channel %d  Calibrated integral samples %d-%d",
                                               runNumber, channel, integralStartSample, integralEndSample),
                               calibratedIntegralBins, calibratedIntegralXMin, calibratedIntegralXMax);
    histogram->SetDirectory(nullptr);
    histogram->SetLineColor(kBlue + 1);
    histogram->SetLineWidth(2);
    histogram->GetXaxis()->SetTitle("Calibrated integrated ADC");
    histogram->GetYaxis()->SetTitle("Entries");
    histogram->GetXaxis()->SetTitleSize(0.045);
    histogram->GetYaxis()->SetTitleSize(0.045);
    histogram->GetXaxis()->SetLabelSize(0.040);
    histogram->GetYaxis()->SetLabelSize(0.040);
    histogram->SetTitleSize(0.055, "t");

    for (double value : calibratedIntegralValues[i]) {
      histogram->Fill(value);
    }
    calibratedIntegralHistograms.push_back(histogram);
  }

  double globalCalibratedIntegralYMaximum = 0.0;
  for (TH1D *histogram : calibratedIntegralHistograms) {
    globalCalibratedIntegralYMaximum = std::max(globalCalibratedIntegralYMaximum, histogram->GetMaximum());
  }
  const double calibratedIntegralYMaximum = 1.15 * globalCalibratedIntegralYMaximum;

  for (int i = 0; i < 16; ++i) {
    const int channel = channels[i];
    TH1D *histogram = calibratedIntegralHistograms[i];

    cCalibratedIntegrals->cd(i + 1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.11);
    gPad->SetTopMargin(0.12);

    histogram->SetMinimum(0.0);
    histogram->SetMaximum(calibratedIntegralYMaximum);

    const int peakBin = histogram->GetMaximumBin();
    const double peakX = histogram->GetBinCenter(peakBin);
    const double rms = histogram->GetRMS();
    const double fitMin = std::max(histogram->GetXaxis()->GetXmin(), peakX - 2.0 * rms);
    const double fitMax = std::min(histogram->GetXaxis()->GetXmax(), peakX + 2.0 * rms);

    TF1 *fit = new TF1(TString::Format("fCalibratedGaussian_ch%d", channel), "gaus", fitMin, fitMax);
    fit->SetLineColor(kRed + 1);
    fit->SetLineWidth(2);
    calibratedGaussianFits.push_back(fit);

    histogram->Fit(fit, "RQ");
    calibratedFitMeans[i] = fit->GetParameter(1);
    calibratedFitMeanErrors[i] = fit->GetParError(1);
    calibratedFitSigmas[i] = fit->GetParameter(2);
    calibratedFitSigmaErrors[i] = fit->GetParError(2);
    histogram->Draw();
    fit->Draw("same");
    gPad->Update();
  }

  cCalibratedIntegrals->Update();

  const TString calibratedIntegralPdfName =
      TString::Format("%s/Run%d_C3_calibrated_integrals_samples_%d_%d.pdf",
                      runOutputDir.Data(), runNumber, integralStartSample, integralEndSample);
  const TString calibratedIntegralPngName =
      TString::Format("%s/Run%d_C3_calibrated_integrals_samples_%d_%d.png",
                      runOutputDir.Data(), runNumber, integralStartSample, integralEndSample);
  cCalibratedIntegrals->SaveAs(calibratedIntegralPdfName);
  cCalibratedIntegrals->SaveAs(calibratedIntegralPngName);

  TCanvas *cResolutionSummary =
      new TCanvas("cResolutionSummary", "Calibrated resolution summary", 1400, 1200);
  cResolutionSummary->Divide(1, 2);

  std::vector<double> graphX(16, 0.0);
  std::vector<double> graphY(16, 0.0);
  std::vector<double> graphXErrors(16, 0.0);
  std::vector<double> graphYErrors(16, 0.0);
  double maxResolutionPercent = 0.0;

  for (int i = 0; i < 16; ++i) {
    graphX[i] = channels[i];
    if (calibratedFitMeans[i] != 0.0) {
      const double resolution = calibratedFitSigmas[i] / calibratedFitMeans[i];
      const double resolutionError =
          TMath::Sqrt(TMath::Power(calibratedFitSigmaErrors[i] / calibratedFitMeans[i], 2) +
                      TMath::Power(calibratedFitSigmas[i] * calibratedFitMeanErrors[i] /
                                       (calibratedFitMeans[i] * calibratedFitMeans[i]),
                                   2));
      graphY[i] = 100.0 * resolution;
      graphYErrors[i] = 100.0 * resolutionError;
      maxResolutionPercent = std::max(maxResolutionPercent, graphY[i] + graphYErrors[i]);
    }
  }

  cResolutionSummary->cd(1);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gPad->SetLeftMargin(0.10);
  gPad->SetRightMargin(0.03);
  gPad->SetBottomMargin(0.12);
  gPad->SetTopMargin(0.10);

  TGraphErrors *resolutionGraph =
      new TGraphErrors(16, graphX.data(), graphY.data(), graphXErrors.data(), graphYErrors.data());
  resolutionGraph->SetTitle(TString::Format("Run %d calibrated channel resolution;Channel;#sigma / #mu [%%]",
                                            runNumber));
  resolutionGraph->SetMarkerStyle(20);
  resolutionGraph->SetMarkerSize(1.2);
  resolutionGraph->SetMarkerColor(kBlue + 1);
  resolutionGraph->SetLineColor(kBlue + 1);
  resolutionGraph->SetLineWidth(2);
  resolutionGraph->Draw("AP");
  resolutionGraph->GetXaxis()->SetLimits(35.0, 53.0);
  resolutionGraph->SetMinimum(0.0);
  resolutionGraph->SetMaximum(1.20 * maxResolutionPercent);
  resolutionGraph->GetXaxis()->SetNdivisions(518, kFALSE);
  resolutionGraph->GetXaxis()->SetTitleSize(0.045);
  resolutionGraph->GetYaxis()->SetTitleSize(0.045);
  resolutionGraph->GetXaxis()->SetLabelSize(0.040);
  resolutionGraph->GetYaxis()->SetLabelSize(0.040);

  cResolutionSummary->cd(2);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gPad->SetLeftMargin(0.10);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.12);
  gPad->SetTopMargin(0.10);

  TH1D *combinedCalibratedHistogram =
      new TH1D("hCombinedCalibratedIntegral",
               TString::Format("Run %d  Combined calibrated integrals, channels 36-43 and 45-52;"
                               "Calibrated integrated ADC;Entries",
                               runNumber),
               calibratedIntegralBins, calibratedIntegralXMin, calibratedIntegralXMax);
  combinedCalibratedHistogram->SetDirectory(nullptr);
  combinedCalibratedHistogram->SetLineColor(kBlue + 1);
  combinedCalibratedHistogram->SetLineWidth(2);
  for (const auto &values : calibratedIntegralValues) {
    for (double value : values) {
      combinedCalibratedHistogram->Fill(value);
    }
  }

  const int combinedPeakBin = combinedCalibratedHistogram->GetMaximumBin();
  const double combinedPeakX = combinedCalibratedHistogram->GetBinCenter(combinedPeakBin);
  const double combinedRms = combinedCalibratedHistogram->GetRMS();
  const double combinedFitMin =
      std::max(combinedCalibratedHistogram->GetXaxis()->GetXmin(), combinedPeakX - 2.0 * combinedRms);
  const double combinedFitMax =
      std::min(combinedCalibratedHistogram->GetXaxis()->GetXmax(), combinedPeakX + 2.0 * combinedRms);

  TF1 *combinedFit = new TF1("fCombinedCalibratedGaussian", "gaus", combinedFitMin, combinedFitMax);
  combinedFit->SetLineColor(kRed + 1);
  combinedFit->SetLineWidth(2);
  combinedCalibratedHistogram->Fit(combinedFit, "RQ");
  combinedCalibratedHistogram->Draw();
  combinedFit->Draw("same");

  const double combinedMean = combinedFit->GetParameter(1);
  const double combinedSigma = combinedFit->GetParameter(2);
  const double combinedMeanError = combinedFit->GetParError(1);
  const double combinedSigmaError = combinedFit->GetParError(2);
  const double combinedResolution = (combinedMean != 0.0) ? combinedSigma / combinedMean : 0.0;
  const double combinedResolutionError =
      (combinedMean != 0.0)
          ? TMath::Sqrt(TMath::Power(combinedSigmaError / combinedMean, 2) +
                        TMath::Power(combinedSigma * combinedMeanError / (combinedMean * combinedMean), 2))
          : 0.0;

  TPaveText *combinedResolutionText = new TPaveText(0.16, 0.72, 0.46, 0.87, "NDC");
  combinedResolutionText->SetFillColor(kWhite);
  combinedResolutionText->SetBorderSize(1);
  combinedResolutionText->SetTextAlign(12);
  combinedResolutionText->AddText(TString::Format("#sigma/#mu = %.4f #pm %.4f",
                                                  combinedResolution, combinedResolutionError));
  combinedResolutionText->AddText(TString::Format("#sigma/#mu = %.2f #pm %.2f %%",
                                                  100.0 * combinedResolution,
                                                  100.0 * combinedResolutionError));
  combinedResolutionText->Draw("same");
  gPad->Update();

  cResolutionSummary->Update();

  const TString resolutionSummaryPdfName =
      TString::Format("%s/Run%d_C4_resolution_summary_samples_%d_%d.pdf",
                      runOutputDir.Data(), runNumber, integralStartSample, integralEndSample);
  const TString resolutionSummaryPngName =
      TString::Format("%s/Run%d_C4_resolution_summary_samples_%d_%d.png",
                      runOutputDir.Data(), runNumber, integralStartSample, integralEndSample);
  cResolutionSummary->SaveAs(resolutionSummaryPdfName);
  cResolutionSummary->SaveAs(resolutionSummaryPngName);

  TCanvas *cTotIntegrals = new TCanvas("cTotIntegrals", "Selected channel TOT integrals", 1600, 1200);
  cTotIntegrals->Divide(4, 4, 0.001, 0.001);

  double globalTotMinimum = totIntegralValues[0][0];
  double globalTotMaximum = totIntegralValues[0][0];
  for (const auto &values : totIntegralValues) {
    const auto minmax = std::minmax_element(values.begin(), values.end());
    globalTotMinimum = std::min(globalTotMinimum, *minmax.first);
    globalTotMaximum = std::max(globalTotMaximum, *minmax.second);
  }
  const double totPadding = std::max(1.0, 0.05 * (globalTotMaximum - globalTotMinimum));
  const double totXMin = globalTotMinimum - totPadding;
  const double totXMax = globalTotMaximum + totPadding;
  const double totBinWidth = (totXMax - totXMin) / totBins;

  std::vector<bool> totChannelActive(16, false);
  std::vector<int> activeTotChannelIndices;
  activeTotChannelIndices.reserve(16);
  for (int i = 0; i < 16; ++i) {
    const auto minmax = std::minmax_element(totIntegralValues[i].begin(), totIntegralValues[i].end());
    if (*minmax.second > 0.0) {
      totChannelActive[i] = true;
      activeTotChannelIndices.push_back(i);
    }
  }

  std::vector<TH1D *> totIntegralHistograms;
  totIntegralHistograms.reserve(16);
  std::vector<TF1 *> totGaussianFits;
  totGaussianFits.reserve(16);
  std::vector<double> totRawFitMeans(16, 0.0);

  for (int i = 0; i < 16; ++i) {
    const int channel = channels[i];
    TH1D *histogram = new TH1D(TString::Format("hTotIntegral_ch%d", channel),
                               TString::Format("Run %d  Channel %d  TOT sum, all 20 samples",
                                               runNumber, channel),
                               totBins, totXMin, totXMax);
    histogram->SetDirectory(nullptr);
    histogram->SetLineColor(kBlue + 1);
    histogram->SetLineWidth(2);
    histogram->GetXaxis()->SetTitle("Summed TOT");
    histogram->GetYaxis()->SetTitle("Entries");
    histogram->GetXaxis()->SetTitleSize(0.045);
    histogram->GetYaxis()->SetTitleSize(0.045);
    histogram->GetXaxis()->SetLabelSize(0.040);
    histogram->GetYaxis()->SetLabelSize(0.040);
    histogram->SetTitleSize(0.060, "t");

    for (double value : totIntegralValues[i]) {
      histogram->Fill(value);
    }
    totIntegralHistograms.push_back(histogram);
  }

  double globalTotYMaximum = 0.0;
  for (TH1D *histogram : totIntegralHistograms) {
    globalTotYMaximum = std::max(globalTotYMaximum, histogram->GetMaximum());
  }
  const double totYMaximum = 1.15 * globalTotYMaximum;

  for (int i = 0; i < 16; ++i) {
    const int channel = channels[i];
    TH1D *histogram = totIntegralHistograms[i];

    cTotIntegrals->cd(i + 1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.11);
    gPad->SetTopMargin(0.12);

    histogram->SetMinimum(0.0);
    histogram->SetMaximum(totYMaximum);

    histogram->Draw();
    if (totChannelActive[i]) {
      const int peakBin = histogram->GetMaximumBin();
      const double peakX = histogram->GetBinCenter(peakBin);
      const double rms = histogram->GetRMS();
      const double fitMin = std::max(histogram->GetXaxis()->GetXmin(), peakX - 2.0 * rms);
      const double fitMax = std::min(histogram->GetXaxis()->GetXmax(), peakX + 2.0 * rms);

      TF1 *fit = new TF1(TString::Format("fTotGaussian_ch%d", channel), "gaus", fitMin, fitMax);
      fit->SetLineColor(kRed + 1);
      fit->SetLineWidth(2);
      totGaussianFits.push_back(fit);

      histogram->Fit(fit, "RQ");
      totRawFitMeans[i] = fit->GetParameter(1);
      fit->Draw("same");
    } else {
      TPaveText *ignoredText = new TPaveText(0.22, 0.72, 0.78, 0.84, "NDC");
      ignoredText->SetFillColor(kWhite);
      ignoredText->SetBorderSize(1);
      ignoredText->SetTextAlign(22);
      ignoredText->AddText("TOT empty");
      ignoredText->AddText("ignored for calibration");
      ignoredText->Draw("same");
    }
    gPad->Update();
  }

  cTotIntegrals->Update();

  const TString totIntegralPdfName =
      TString::Format("%s/Run%d_C5_tot_integrals_all_samples.pdf", runOutputDir.Data(), runNumber);
  const TString totIntegralPngName =
      TString::Format("%s/Run%d_C5_tot_integrals_all_samples.png", runOutputDir.Data(), runNumber);
  cTotIntegrals->SaveAs(totIntegralPdfName);
  cTotIntegrals->SaveAs(totIntegralPngName);

  std::vector<int> validTotChannelIndices;
  validTotChannelIndices.reserve(activeTotChannelIndices.size());
  double averageTotRawPeak = 0.0;
  for (int channelIndex : activeTotChannelIndices) {
    if (totRawFitMeans[channelIndex] > 0.0) {
      averageTotRawPeak += totRawFitMeans[channelIndex];
      validTotChannelIndices.push_back(channelIndex);
    }
  }

  TString calibratedTotPdfName;
  TString calibratedTotPngName;
  TString totResolutionSummaryPdfName;
  TString totResolutionSummaryPngName;
  double combinedTotResolution = 0.0;
  double combinedTotResolutionError = 0.0;

  if (!validTotChannelIndices.empty()) {
    averageTotRawPeak /= validTotChannelIndices.size();

    std::vector<std::vector<double>> calibratedTotValues(16);
    for (int channelIndex : validTotChannelIndices) {
      const double calibrationFactor = averageTotRawPeak / totRawFitMeans[channelIndex];
      calibratedTotValues[channelIndex].reserve(totIntegralValues[channelIndex].size());
      for (double value : totIntegralValues[channelIndex]) {
        calibratedTotValues[channelIndex].push_back(value * calibrationFactor);
      }
    }

    bool firstCalibratedTotValue = true;
    double globalCalibratedTotMinimum = 0.0;
    double globalCalibratedTotMaximum = 0.0;
    for (int channelIndex : validTotChannelIndices) {
      const auto minmax = std::minmax_element(calibratedTotValues[channelIndex].begin(),
                                              calibratedTotValues[channelIndex].end());
      if (firstCalibratedTotValue) {
        globalCalibratedTotMinimum = *minmax.first;
        globalCalibratedTotMaximum = *minmax.second;
        firstCalibratedTotValue = false;
      } else {
        globalCalibratedTotMinimum = std::min(globalCalibratedTotMinimum, *minmax.first);
        globalCalibratedTotMaximum = std::max(globalCalibratedTotMaximum, *minmax.second);
      }
    }

    const double calibratedTotPadding =
        std::max(1.0, 0.05 * (globalCalibratedTotMaximum - globalCalibratedTotMinimum));
    const double calibratedTotXMin = globalCalibratedTotMinimum - calibratedTotPadding;
    const double calibratedTotXMax = globalCalibratedTotMaximum + calibratedTotPadding;
    const int calibratedTotBins =
        std::max(20, TMath::CeilNint((calibratedTotXMax - calibratedTotXMin) / totBinWidth));

    TCanvas *cCalibratedTotIntegrals =
        new TCanvas("cCalibratedTotIntegrals", "Calibrated selected channel TOT integrals", 1600, 1200);
    cCalibratedTotIntegrals->Divide(4, 4, 0.001, 0.001);

    std::vector<TH1D *> calibratedTotHistograms;
    calibratedTotHistograms.reserve(16);
    std::vector<TF1 *> calibratedTotFits;
    calibratedTotFits.reserve(validTotChannelIndices.size());
    std::vector<double> calibratedTotFitMeans(16, 0.0);
    std::vector<double> calibratedTotFitMeanErrors(16, 0.0);
    std::vector<double> calibratedTotFitSigmas(16, 0.0);
    std::vector<double> calibratedTotFitSigmaErrors(16, 0.0);

    for (int i = 0; i < 16; ++i) {
      const int channel = channels[i];
      TH1D *histogram = new TH1D(TString::Format("hCalibratedTotIntegral_ch%d", channel),
                                 TString::Format("Run %d  Channel %d  Calibrated TOT sum",
                                                 runNumber, channel),
                                 calibratedTotBins, calibratedTotXMin, calibratedTotXMax);
      histogram->SetDirectory(nullptr);
      histogram->SetLineColor(kBlue + 1);
      histogram->SetLineWidth(2);
      histogram->GetXaxis()->SetTitle("Calibrated summed TOT");
      histogram->GetYaxis()->SetTitle("Entries");
      histogram->GetXaxis()->SetTitleSize(0.045);
      histogram->GetYaxis()->SetTitleSize(0.045);
      histogram->GetXaxis()->SetLabelSize(0.040);
      histogram->GetYaxis()->SetLabelSize(0.040);
      histogram->SetTitleSize(0.055, "t");

      if (totChannelActive[i] && totRawFitMeans[i] > 0.0) {
        for (double value : calibratedTotValues[i]) {
          histogram->Fill(value);
        }
      }
      calibratedTotHistograms.push_back(histogram);
    }

    double globalCalibratedTotYMaximum = 0.0;
    for (TH1D *histogram : calibratedTotHistograms) {
      globalCalibratedTotYMaximum = std::max(globalCalibratedTotYMaximum, histogram->GetMaximum());
    }
    const double calibratedTotYMaximum = 1.15 * globalCalibratedTotYMaximum;

    for (int i = 0; i < 16; ++i) {
      const int channel = channels[i];
      TH1D *histogram = calibratedTotHistograms[i];

      cCalibratedTotIntegrals->cd(i + 1);
      gPad->SetGridx(1);
      gPad->SetGridy(1);
      gPad->SetRightMargin(0.05);
      gPad->SetLeftMargin(0.12);
      gPad->SetBottomMargin(0.11);
      gPad->SetTopMargin(0.12);

      histogram->SetMinimum(0.0);
      histogram->SetMaximum(calibratedTotYMaximum);
      histogram->Draw();

      if (totChannelActive[i] && totRawFitMeans[i] > 0.0) {
        const int peakBin = histogram->GetMaximumBin();
        const double peakX = histogram->GetBinCenter(peakBin);
        const double rms = histogram->GetRMS();
        const double fitMin = std::max(histogram->GetXaxis()->GetXmin(), peakX - 2.0 * rms);
        const double fitMax = std::min(histogram->GetXaxis()->GetXmax(), peakX + 2.0 * rms);

        TF1 *fit = new TF1(TString::Format("fCalibratedTotGaussian_ch%d", channel), "gaus", fitMin, fitMax);
        fit->SetLineColor(kRed + 1);
        fit->SetLineWidth(2);
        calibratedTotFits.push_back(fit);

        histogram->Fit(fit, "RQ");
        calibratedTotFitMeans[i] = fit->GetParameter(1);
        calibratedTotFitMeanErrors[i] = fit->GetParError(1);
        calibratedTotFitSigmas[i] = fit->GetParameter(2);
        calibratedTotFitSigmaErrors[i] = fit->GetParError(2);
        fit->Draw("same");
      } else {
        TPaveText *ignoredText = new TPaveText(0.22, 0.72, 0.78, 0.84, "NDC");
        ignoredText->SetFillColor(kWhite);
        ignoredText->SetBorderSize(1);
        ignoredText->SetTextAlign(22);
        ignoredText->AddText("TOT empty");
        ignoredText->AddText("ignored");
        ignoredText->Draw("same");
      }
      gPad->Update();
    }

    cCalibratedTotIntegrals->Update();

    calibratedTotPdfName =
        TString::Format("%s/Run%d_C6_calibrated_tot_integrals_all_samples.pdf",
                        runOutputDir.Data(), runNumber);
    calibratedTotPngName =
        TString::Format("%s/Run%d_C6_calibrated_tot_integrals_all_samples.png",
                        runOutputDir.Data(), runNumber);
    cCalibratedTotIntegrals->SaveAs(calibratedTotPdfName);
    cCalibratedTotIntegrals->SaveAs(calibratedTotPngName);

    TCanvas *cTotResolutionSummary =
        new TCanvas("cTotResolutionSummary", "Calibrated TOT resolution summary", 1400, 1200);
    cTotResolutionSummary->Divide(1, 2);

    std::vector<double> totGraphX;
    std::vector<double> totGraphY;
    std::vector<double> totGraphXErrors;
    std::vector<double> totGraphYErrors;
    totGraphX.reserve(validTotChannelIndices.size());
    totGraphY.reserve(validTotChannelIndices.size());
    totGraphXErrors.reserve(validTotChannelIndices.size());
    totGraphYErrors.reserve(validTotChannelIndices.size());
    double maxTotResolutionPercent = 0.0;

    for (int channelIndex : validTotChannelIndices) {
      if (calibratedTotFitMeans[channelIndex] == 0.0) {
        continue;
      }
      const double resolution = calibratedTotFitSigmas[channelIndex] / calibratedTotFitMeans[channelIndex];
      const double resolutionError =
          TMath::Sqrt(TMath::Power(calibratedTotFitSigmaErrors[channelIndex] /
                                       calibratedTotFitMeans[channelIndex],
                                   2) +
                      TMath::Power(calibratedTotFitSigmas[channelIndex] *
                                       calibratedTotFitMeanErrors[channelIndex] /
                                       (calibratedTotFitMeans[channelIndex] *
                                        calibratedTotFitMeans[channelIndex]),
                                   2));
      totGraphX.push_back(channels[channelIndex]);
      totGraphY.push_back(100.0 * resolution);
      totGraphXErrors.push_back(0.0);
      totGraphYErrors.push_back(100.0 * resolutionError);
      maxTotResolutionPercent = std::max(maxTotResolutionPercent, 100.0 * resolution + 100.0 * resolutionError);
    }

    cTotResolutionSummary->cd(1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.03);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.10);

    TGraphErrors *totResolutionGraph =
        new TGraphErrors(totGraphX.size(), totGraphX.data(), totGraphY.data(),
                         totGraphXErrors.data(), totGraphYErrors.data());
    totResolutionGraph->SetTitle(TString::Format("Run %d calibrated TOT channel resolution;Channel;#sigma / #mu [%%]",
                                                 runNumber));
    totResolutionGraph->SetMarkerStyle(20);
    totResolutionGraph->SetMarkerSize(1.2);
    totResolutionGraph->SetMarkerColor(kBlue + 1);
    totResolutionGraph->SetLineColor(kBlue + 1);
    totResolutionGraph->SetLineWidth(2);
    totResolutionGraph->Draw("AP");
    totResolutionGraph->GetXaxis()->SetLimits(35.0, 53.0);
    totResolutionGraph->SetMinimum(0.0);
    totResolutionGraph->SetMaximum(1.20 * maxTotResolutionPercent);
    totResolutionGraph->GetXaxis()->SetNdivisions(518, kFALSE);
    totResolutionGraph->GetXaxis()->SetTitleSize(0.045);
    totResolutionGraph->GetYaxis()->SetTitleSize(0.045);
    totResolutionGraph->GetXaxis()->SetLabelSize(0.040);
    totResolutionGraph->GetYaxis()->SetLabelSize(0.040);

    cTotResolutionSummary->cd(2);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.10);

    TH1D *combinedCalibratedTotHistogram =
        new TH1D("hCombinedCalibratedTotIntegral",
                 TString::Format("Run %d  Combined calibrated TOT sums, active selected channels;"
                                 "Calibrated summed TOT;Entries",
                                 runNumber),
                 calibratedTotBins, calibratedTotXMin, calibratedTotXMax);
    combinedCalibratedTotHistogram->SetDirectory(nullptr);
    combinedCalibratedTotHistogram->SetLineColor(kBlue + 1);
    combinedCalibratedTotHistogram->SetLineWidth(2);

    for (int channelIndex : validTotChannelIndices) {
      for (double value : calibratedTotValues[channelIndex]) {
        combinedCalibratedTotHistogram->Fill(value);
      }
    }

    const int combinedTotPeakBin = combinedCalibratedTotHistogram->GetMaximumBin();
    const double combinedTotPeakX = combinedCalibratedTotHistogram->GetBinCenter(combinedTotPeakBin);
    const double combinedTotRms = combinedCalibratedTotHistogram->GetRMS();
    const double combinedTotFitMin =
        std::max(combinedCalibratedTotHistogram->GetXaxis()->GetXmin(), combinedTotPeakX - 2.0 * combinedTotRms);
    const double combinedTotFitMax =
        std::min(combinedCalibratedTotHistogram->GetXaxis()->GetXmax(), combinedTotPeakX + 2.0 * combinedTotRms);

    TF1 *combinedTotFit = new TF1("fCombinedCalibratedTotGaussian", "gaus",
                                  combinedTotFitMin, combinedTotFitMax);
    combinedTotFit->SetLineColor(kRed + 1);
    combinedTotFit->SetLineWidth(2);
    combinedCalibratedTotHistogram->Fit(combinedTotFit, "RQ");
    combinedCalibratedTotHistogram->Draw();
    combinedTotFit->Draw("same");

    const double combinedTotMean = combinedTotFit->GetParameter(1);
    const double combinedTotSigma = combinedTotFit->GetParameter(2);
    const double combinedTotMeanError = combinedTotFit->GetParError(1);
    const double combinedTotSigmaError = combinedTotFit->GetParError(2);
    combinedTotResolution = (combinedTotMean != 0.0) ? combinedTotSigma / combinedTotMean : 0.0;
    combinedTotResolutionError =
        (combinedTotMean != 0.0)
            ? TMath::Sqrt(TMath::Power(combinedTotSigmaError / combinedTotMean, 2) +
                          TMath::Power(combinedTotSigma * combinedTotMeanError /
                                         (combinedTotMean * combinedTotMean),
                                       2))
            : 0.0;

    TPaveText *combinedTotResolutionText = new TPaveText(0.16, 0.72, 0.46, 0.87, "NDC");
    combinedTotResolutionText->SetFillColor(kWhite);
    combinedTotResolutionText->SetBorderSize(1);
    combinedTotResolutionText->SetTextAlign(12);
    combinedTotResolutionText->AddText(TString::Format("#sigma/#mu = %.4f #pm %.4f",
                                                       combinedTotResolution, combinedTotResolutionError));
    combinedTotResolutionText->AddText(TString::Format("#sigma/#mu = %.2f #pm %.2f %%",
                                                       100.0 * combinedTotResolution,
                                                       100.0 * combinedTotResolutionError));
    combinedTotResolutionText->Draw("same");
    gPad->Update();

    cTotResolutionSummary->Update();

    totResolutionSummaryPdfName =
        TString::Format("%s/Run%d_C7_tot_resolution_summary_all_samples.pdf",
                        runOutputDir.Data(), runNumber);
    totResolutionSummaryPngName =
        TString::Format("%s/Run%d_C7_tot_resolution_summary_all_samples.png",
                        runOutputDir.Data(), runNumber);
    cTotResolutionSummary->SaveAs(totResolutionSummaryPdfName);
    cTotResolutionSummary->SaveAs(totResolutionSummaryPngName);
  } else {
    std::cout << "No non-empty TOT channels found. Skipping TOT calibration and summary canvases." << std::endl;
  }

  std::cout << "Processed " << inputFileName << std::endl;
  std::cout << "Wrote " << pdfName << std::endl;
  std::cout << "Wrote " << pngName << std::endl;
  std::cout << "Wrote " << integralPdfName << std::endl;
  std::cout << "Wrote " << integralPngName << std::endl;
  std::cout << "Average raw peak used for calibration: " << averageRawPeak << std::endl;
  std::cout << "Wrote " << calibratedIntegralPdfName << std::endl;
  std::cout << "Wrote " << calibratedIntegralPngName << std::endl;
  std::cout << "Combined calibrated resolution: " << combinedResolution
            << " +/- " << combinedResolutionError << std::endl;
  std::cout << "Wrote " << resolutionSummaryPdfName << std::endl;
  std::cout << "Wrote " << resolutionSummaryPngName << std::endl;
  std::cout << "Wrote " << totIntegralPdfName << std::endl;
  std::cout << "Wrote " << totIntegralPngName << std::endl;
  if (calibratedTotPdfName.Length() > 0) {
    std::cout << "Average non-empty TOT peak used for calibration: " << averageTotRawPeak << std::endl;
    std::cout << "Wrote " << calibratedTotPdfName << std::endl;
    std::cout << "Wrote " << calibratedTotPngName << std::endl;
    std::cout << "Combined calibrated TOT resolution: " << combinedTotResolution
              << " +/- " << combinedTotResolutionError << std::endl;
    std::cout << "Wrote " << totResolutionSummaryPdfName << std::endl;
    std::cout << "Wrote " << totResolutionSummaryPngName << std::endl;
  }

  inputFile->Close();
}
