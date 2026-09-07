// Draw reviewer-ready logger1 temperature plots for the 2025 November beam period.
// Execute from this directory:
// /sw/packages/root/6.32.04/bin/root -l -q plot_logger1_temeparue.cc

#include <TBox.h>
#include <TCanvas.h>
#include <TDatime.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TTree.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGaxis.h>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <set>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

struct TemperatureRange {
  double min = std::numeric_limits<double>::max();
  double max = std::numeric_limits<double>::lowest();
  double sum = 0.;
  int count = 0;

  void Add(double value) {
    if (value < min) min = value;
    if (value > max) max = value;
    sum += value;
    ++count;
  }
  double Mean() const { return sum / count; }
};

struct PedestalPoint { double beamDay, peak, peakError, width, widthError; int run; };
struct RunTime { double on = 0.; double off = 0.; };

std::vector<PedestalPoint> FitPedestals(const char *runListFile, const char *commentFile,
                                        const char *rootDirectory, int minEntries,
                                        double fitHalfWidth, int firstDay)
{
  std::set<int> requestedRuns;
  std::ifstream runList(runListFile);
  std::string line;
  while (std::getline(runList, line)) {
    int run = -1;
    if (sscanf(line.c_str(), " %d:", &run) == 1 && run > 0) requestedRuns.insert(run);
  }

  std::map<int, RunTime> runTimes;
  std::ifstream comments(commentFile);
  while (std::getline(comments, line)) {
    int year, month, day, hour, minute, second, run;
    char state[32] = {};
    if (sscanf(line.c_str(), "%d %d/%d %d:%d:%d [RUN %d] %31s",
               &year, &month, &day, &hour, &minute, &second, &run, state) != 8) continue;
    if (!requestedRuns.count(run)) continue;
    const double stamp = TDatime(year, month, day, hour, minute, second).Convert();
    if (std::string(state) == "G_ON") {
      if (runTimes[run].on == 0.) runTimes[run].on = stamp;
    } else if (std::string(state) == "G_OFF") {
      runTimes[run].off = stamp;
    }
  }

  std::vector<PedestalPoint> result;
  const double beamStart = TDatime(2025, 11, firstDay, 0, 0, 0).Convert();
  for (const int run : requestedRuns) {
    const auto timeIt = runTimes.find(run);
    if (timeIt == runTimes.end() || timeIt->second.on == 0. || timeIt->second.off == 0.) continue;
    const TString fileName = Form("%s/run%05d_Hodoscope.root", rootDirectory, run);
    if (gSystem->AccessPathName(fileName)) continue;
    TFile file(fileName, "READ");
    if (file.IsZombie()) continue;
    TH1D *storedHistogram = dynamic_cast<TH1D *>(file.Get("BAC_ADC_seg4U"));
    if (!storedHistogram || storedHistogram->GetEntries() < minEntries) continue;
    TH1D histogram(*storedHistogram);
    histogram.SetName(Form("pedestal_run%d", run));
    histogram.SetDirectory(nullptr);
    const double seed = histogram.GetBinCenter(histogram.GetMaximumBin());
    TF1 fit(Form("pedestal_fit_run%d", run), "gaus", seed - fitHalfWidth, seed + fitHalfWidth);
    if (histogram.Fit(&fit, "QNR") != 0 || fit.GetParameter(2) <= 0.) continue;
    const double midpoint = 0.5 * (timeIt->second.on + timeIt->second.off);
    result.push_back({1. + (midpoint - beamStart) / 86400., fit.GetParameter(1),
                      fit.GetParError(1), fit.GetParameter(2), fit.GetParError(2), run});
  }
  std::sort(result.begin(), result.end(), [](const PedestalPoint &a, const PedestalPoint &b) {
    return a.beamDay < b.beamDay;
  });
  return result;
}

void AppendPedestalPages(TCanvas &canvas, const char *outputPdf,
                         const std::vector<double> &temperatureTime,
                         const std::vector<double> &temperature,
                         const std::vector<PedestalPoint> &pedestals,
                         int firstDay, int lastDay, double npeFactor, const char *rootDirectory, double fitHalfWidth, double maxPeakErrorNpe, double maxWidthErrorNpe)
{
  if (pedestals.empty()) return;
  const double windowStartDay = 7.8, windowEndDay = 8.8;
  std::vector<PedestalPoint> selected;
  for (const auto &point : pedestals) if (point.beamDay >= windowStartDay && point.beamDay <= windowEndDay && point.peakError * npeFactor <= maxPeakErrorNpe && point.widthError * npeFactor <= maxWidthErrorNpe) selected.push_back(point);
  if (selected.empty()) return;
  std::vector<double> x, delta, width, noXError;
  double meanPeak=0.; for(const auto &p:selected) meanPeak+=p.peak*npeFactor; meanPeak/=selected.size();
  for(const auto &p:selected){x.push_back((p.beamDay-windowStartDay)*24.);delta.push_back(p.peak*npeFactor-meanPeak);width.push_back(p.width*npeFactor);noXError.push_back(0.);}
  std::vector<double> temperatureHour; for(double t:temperatureTime) temperatureHour.push_back((t-windowStartDay)*24.);
  canvas.Clear(); canvas.SetGrid(); gPad->SetRightMargin(0.18);
  TGraph temperatureGraph(temperatureHour.size(),temperatureHour.data(),temperature.data());
  temperatureGraph.SetTitle("Temperature and BAC pedestal variation;Elapsed time [h];Temperature [^{#circ}C]");
  temperatureGraph.SetMarkerStyle(20);temperatureGraph.SetMarkerSize(0.95);temperatureGraph.SetMarkerColor(kBlack);temperatureGraph.Draw("AP");
  temperatureGraph.GetYaxis()->SetRangeUser(20.,25.); temperatureGraph.GetYaxis()->SetNdivisions(505); temperatureGraph.GetYaxis()->SetDecimals(kFALSE); temperatureGraph.GetXaxis()->SetLimits(0.,24.); temperatureGraph.GetXaxis()->SetLabelFont(132); temperatureGraph.GetYaxis()->SetLabelFont(132); temperatureGraph.GetXaxis()->SetTitleFont(132); temperatureGraph.GetYaxis()->SetTitleFont(132);
  const double scale=2.5/10.;std::vector<double> scaled,scaledWidth;for(size_t i=0;i<delta.size();++i){scaled.push_back(22.5+delta[i]*scale);scaledWidth.push_back(width[i]*scale);}
  TGraphErrors widthBand(x.size(),x.data(),scaled.data(),noXError.data(),scaledWidth.data()); widthBand.SetFillColorAlpha(kBlue-9,0.40); widthBand.SetLineColor(kBlue+1); widthBand.Draw("3 same"); TLegend widthLegend(0.18,0.78,0.42,0.86); widthLegend.SetBorderSize(0); widthLegend.SetFillColor(kWhite); widthLegend.SetFillStyle(1001); TBox widthLegendBox(0.,0.,1.,1.); widthLegendBox.SetFillColorAlpha(kBlue-9,0.40); widthLegendBox.SetLineColor(kWhite); widthLegend.AddEntry(&widthLegendBox,"Pedestal Width","f"); widthLegend.Draw(); TGraph pedestalGraph(x.size(),x.data(),scaled.data()); pedestalGraph.SetMarkerStyle(24); pedestalGraph.SetMarkerSize(0.95); pedestalGraph.SetMarkerColor(kBlue+1); pedestalGraph.Draw("P same");
  TGaxis rightAxis(24.,20.,24.,25.,-10.,10.,510,"+R"); rightAxis.SetTitle(""); rightAxis.SetLabelColor(kBlue+1); rightAxis.SetLineColor(kBlue+1); rightAxis.SetLabelOffset(0.050); rightAxis.SetLabelSize(0.045); rightAxis.SetLabelFont(132); rightAxis.Draw(); TLatex rightTitle; rightTitle.SetNDC(); rightTitle.SetTextFont(132); rightTitle.SetTextColor(kBlue+1); rightTitle.SetTextSize(0.060); rightTitle.SetTextAngle(270); rightTitle.DrawLatex(0.91,0.92,"Pedestal Mean [N_{p.e.}]");
  canvas.Print(outputPdf);

  // Gaussian-fit distributions: up to 16 completed runs per PDF page.
  for (size_t first = 0; first < selected.size(); first += 16) {
    canvas.Clear();
    canvas.Divide(4, 4);
    const size_t last = std::min(first + 16, selected.size());
    for (size_t index = first; index < last; ++index) {
      const PedestalPoint &point = selected[index];
      canvas.cd(index - first + 1);
      gPad->SetGrid();
      TFile file(Form("%s/run%05d_Hodoscope.root", rootDirectory, point.run), "READ");
      TH1D *source = dynamic_cast<TH1D *>(file.Get("BAC_ADC_seg4U"));
      if (!source) continue;
      TH1D *histogram = dynamic_cast<TH1D *>(source->Clone(Form("pedestal_panel_run%d", point.run)));
      histogram->SetDirectory(nullptr);
      histogram->SetTitle(Form("run%05d: #mu=%.3f, #sigma=%.3f NPE;ADC count;Entries", point.run, point.peak * npeFactor, point.width * npeFactor));
      histogram->SetLineColor(kBlack);
      histogram->GetXaxis()->SetRangeUser(std::max(0., point.peak - 5. * point.width), point.peak + 5. * point.width);
      histogram->Draw();
      TF1 *fit = new TF1(Form("pedestal_panel_fit_run%d", point.run), "gaus",
                          point.peak - fitHalfWidth, point.peak + fitHalfWidth);
      fit->SetParameters(histogram->GetMaximum(), point.peak, point.width);
      fit->SetLineColor(kRed + 1);
      histogram->Fit(fit, "QR");
    }
    canvas.Print(outputPdf);
  }
}

void plot_logger1_temeparue()
{
  // ===== User-configurable parameters =====
  const char *inputDirectory = "/gpfs/group/had/sks/E72/JPARC2025Nov/share/pg-monitor/csv/gl840";
  const char *outputPdf = "plot_logger1_temeparue.pdf";
  const char *pedestalPdf = "plot_logger1_pedestal_run02263-run02468.pdf";
  const std::string loggerName = "logger1.monitor.k18br";
  const int temperatureValueColumn = 8; // zero-based CSV index: ch04_value (second temperature)
  const int firstDay = 10;              // official experiment period, inclusive
  const int lastDay = 25;               // official experiment period, inclusive
  const int averageIntervalMinutes = 30;
  const double temperatureMin = 19.0; // degC
  const double temperatureMax = 23.0; // degC
  const char *runListFile = "/home/had/haein/work/e72/ana/e72/runmanager/runlist/hodo_735.yml";
  const char *commentFile = "/home/had/haein/raw_data/JPARC2025Nov/e72_2026apr/misc/comment.txt";
  const char *hodoscopeDirectory = "/home/had/haein/data/JPARC2025Nov_root/hodo/mom-735";
  const int minPedestalEntries = 1000;
  const double pedestalFitHalfWidth = 25.0; // ADC count
  const double npeFactor = 37.7331 / 398.085; // NPE per ADC count, beta_compare.cc
  const double maxPedestalPeakErrorNpe = 1.0;
  const double maxPedestalWidthErrorNpe = 1.0;
  // ========================================

  const std::string command = "gzip -cd " + std::string(inputDirectory)
                            + "/gl840_2025-11-*.csv.gz";
  FILE *pipe = popen(command.c_str(), "r");
  if (!pipe) {
    Error("plot_logger1_temeparue", "Cannot read November 2025 CSV files");
    return;
  }

  // Key: 10-minute-bin start in Unix seconds. Value: (temperature sum, sample count).
  const long long intervalSeconds = 60LL * averageIntervalMinutes;
  std::map<long long, std::pair<double, int>> bins;
  char buffer[8192];
  while (fgets(buffer, sizeof(buffer), pipe)) {
    std::vector<std::string> field;
    std::stringstream stream(buffer);
    std::string item;
    while (std::getline(stream, item, ',')) field.push_back(item);
    if (field.size() <= static_cast<size_t>(temperatureValueColumn) || field[1] != loggerName)
      continue; // skips headers and non-logger1 rows

    int year, month, day, hour, minute, second;
    if (sscanf(field[0].c_str(), "%d-%d-%d %d:%d:%d",
               &year, &month, &day, &hour, &minute, &second) != 6)
      continue;
    if (year != 2025 || month != 11 || day < firstDay || day > lastDay)
      continue;
    try {
      TDatime timestamp(year, month, day, hour, minute, second);
      const long long binStart = (timestamp.Convert() / intervalSeconds) * intervalSeconds;
      auto &bin = bins[binStart];
      bin.first += std::stod(field[temperatureValueColumn]);
      ++bin.second;
    } catch (...) {
      Warning("plot_logger1_temeparue", "Skipping malformed CSV value");
    }
  }
  pclose(pipe);

  if (bins.empty()) {
    Error("plot_logger1_temeparue", "No logger1 ch04 values were found");
    return;
  }

  std::vector<double> allTime, allBeamDay, allTemperature;
  std::map<int, std::vector<std::pair<double, double>>> dailyPoints;
  std::map<int, TemperatureRange> dailyRange;
  std::map<int, std::map<int, TemperatureRange>> hourlyRange;
  for (const auto &entry : bins) {
    const double binCenter = entry.first + intervalSeconds / 2.0;
    const double meanTemperature = entry.second.first / entry.second.second;
    TDatime dateTime(static_cast<UInt_t>(entry.first));
    const int day = dateTime.GetDay();
    const int hour = dateTime.GetHour();
    allTime.push_back(binCenter);
    const double beamStart = TDatime(2025, 11, firstDay, 0, 0, 0).Convert();
    allBeamDay.push_back(1. + (binCenter - beamStart) / (24. * 60. * 60.));
    allTemperature.push_back(meanTemperature);
    dailyPoints[day].emplace_back(binCenter, meanTemperature);
    dailyRange[day].Add(meanTemperature);
    hourlyRange[day][hour].Add(meanTemperature);
  }

  const std::vector<PedestalPoint> pedestals = FitPedestals(runListFile, commentFile, hodoscopeDirectory, minPedestalEntries, pedestalFitHalfWidth, firstDay);
  gStyle->SetOptStat(0);
  gStyle->SetTextFont(132); gStyle->SetLabelFont(132, "XYZ"); gStyle->SetTitleFont(132, "XYZ");
  gStyle->SetLegendFont(132); gROOT->ForceStyle();
  TCanvas canvas("canvas", "logger1 temperature", 1200, 800);
  canvas.Print((std::string(outputPdf) + "[").c_str());

  TLatex info;
  info.SetNDC();
  info.SetTextFont(132);
  info.SetTextSize(0.040);
  info.DrawLatex(0.10, 0.86, "Logger temperature monitor: official experiment period");
  info.SetTextSize(0.028);
  info.DrawLatex(0.10, 0.74, "Macro: plot_logger1_temeparue.cc");
  info.DrawLatex(0.10, 0.67, "Execution date: 2026-09-07");
  info.DrawLatex(0.10, 0.60, "Input dates: 2025-11-10 through 2025-11-25 (JST, UTC+09)");
  info.DrawLatex(0.10, 0.53, "Logger: logger1.monitor.k18br; quantity: ch04_value (second temperature) [degC]");
  info.DrawLatex(0.10, 0.46, Form("Points are %d-minute fixed-bin means.", averageIntervalMinutes));
  info.DrawLatex(0.10, 0.39, "Colored boxes show the min--max spread of the 10-minute means.");
  info.DrawLatex(0.10, 0.32, "Y-axis range on every graph: 20 to 25 degC");
  info.DrawLatex(0.10, 0.25, Form("Averaged points: %zu", allTime.size()));
  info.DrawLatex(0.10, 0.18, "Pedestal: Gaussian fit to stored BAC_ADC_seg4U (= bac_adc_u[4]) histogram.");
  info.DrawLatex(0.10, 0.11, "Runs: hodo_735.yml; only files with valid G_ON/G_OFF and ROOT output are used.");
  canvas.Print(outputPdf);

  auto setTimeAxis = [](TGraph &graph, const char *timeFormat) {
    graph.GetXaxis()->SetTimeDisplay(1);
    graph.GetXaxis()->SetTimeFormat(timeFormat);
    graph.GetXaxis()->SetTimeOffset(0, "local");
  };

  // Page 2: all averaged points versus elapsed beam-period time.
  std::vector<double> allElapsedDay = allBeamDay;
  for (double &elapsedDay : allElapsedDay) elapsedDay -= 1.;
  canvas.Clear();
  canvas.SetGrid();
  TGraph monthlyGraph(allElapsedDay.size(), allElapsedDay.data(), allTemperature.data());
  monthlyGraph.SetTitle("logger1 second temperature channel during beam period;Elapsed time [day];Temperature [^{#circ}C]");
  monthlyGraph.SetMarkerStyle(20);
  monthlyGraph.SetMarkerSize(0.62);
  monthlyGraph.SetMarkerColor(kBlack);
  monthlyGraph.Draw("AP");
  monthlyGraph.GetYaxis()->SetRangeUser(temperatureMin, temperatureMax);
  monthlyGraph.GetYaxis()->SetNdivisions(505);
  monthlyGraph.GetYaxis()->SetDecimals(kFALSE);
  monthlyGraph.GetXaxis()->SetLimits(-0.5, lastDay - firstDay + 0.5);
  monthlyGraph.GetXaxis()->SetNdivisions(510);
  canvas.Print(outputPdf);

  // Page 3: one colored daily min--max box and one mean point for each day.
  std::vector<double> dayTime, dayMean;
  for (int day = firstDay; day <= lastDay; ++day) {
    if (!dailyRange.count(day)) continue;
    dayTime.push_back(day - firstDay);
    dayMean.push_back(dailyRange[day].Mean());
  }
  canvas.Clear();
  canvas.SetGrid();
  TGraph dailyMeanGraph(dayTime.size(), dayTime.data(), dayMean.data());
  dailyMeanGraph.SetTitle("Daily temperature range during beam period;Elapsed time [day];Temperature [^{#circ}C]");
  dailyMeanGraph.SetMarkerStyle(20);
  dailyMeanGraph.SetMarkerSize(1.10);
  dailyMeanGraph.SetMarkerColor(kBlack);
  dailyMeanGraph.Draw("AP");
  dailyMeanGraph.GetYaxis()->SetRangeUser(temperatureMin, temperatureMax);
  dailyMeanGraph.GetYaxis()->SetNdivisions(505);
  dailyMeanGraph.GetYaxis()->SetDecimals(kFALSE);
  dailyMeanGraph.GetXaxis()->SetLimits(-0.5, lastDay - firstDay + 0.5);
  dailyMeanGraph.GetXaxis()->SetNdivisions(510);
  std::vector<TBox> dailyBoxes;
  dailyBoxes.reserve(lastDay - firstDay + 1);
  const double halfDayWidth = 0.42;
  for (int day = firstDay; day <= lastDay; ++day) {
    if (!dailyRange.count(day)) continue;
    const double center = day - firstDay;
    dailyBoxes.emplace_back(center - halfDayWidth, std::max(temperatureMin, dailyRange[day].min),
                            center + halfDayWidth, std::min(temperatureMax, dailyRange[day].max));
    dailyBoxes.back().SetFillColorAlpha(kAzure - 9, 0.55);
    dailyBoxes.back().SetLineColor(kWhite);
    dailyBoxes.back().Draw("same");
  }
  dailyMeanGraph.Draw("P same");
  if (!pedestals.empty()) {
    double minimum = pedestals.front().peak * npeFactor, maximum = minimum;
    for (const auto &point : pedestals) { minimum = std::min(minimum, point.peak * npeFactor); maximum = std::max(maximum, point.peak * npeFactor); }
    if (maximum <= minimum) maximum = minimum + 1.;
    std::vector<double> xPed, yPed, exPed, eyPed;
    for (const auto &point : pedestals) {
      xPed.push_back(point.beamDay - 1.); exPed.push_back(0.);
      yPed.push_back(temperatureMin + (point.peak * npeFactor - minimum) * (temperatureMax - temperatureMin) / (maximum - minimum));
      eyPed.push_back(point.peakError * npeFactor * (temperatureMax - temperatureMin) / (maximum - minimum));
    }
    TGraphErrors pedestalGraph(xPed.size(), xPed.data(), yPed.data(), exPed.data(), eyPed.data());
    pedestalGraph.SetMarkerStyle(21); pedestalGraph.SetMarkerSize(0.85); pedestalGraph.SetMarkerColor(kRed + 1); pedestalGraph.SetLineColor(kRed + 1);
    pedestalGraph.Draw("P same");
    TGaxis pedestalAxis(lastDay - firstDay + 1.5, temperatureMin, lastDay - firstDay + 1.5, temperatureMax, minimum, maximum, 510, "+L");
    pedestalAxis.SetTitle("Pedestal peak [NPE]"); pedestalAxis.SetLabelColor(kRed + 1); pedestalAxis.SetTitleColor(kRed + 1); pedestalAxis.Draw();
  }
  TBox dailyLegendBox(0., 0., 1., 1.);
  dailyLegendBox.SetFillColorAlpha(kAzure - 9, 0.55);
  dailyLegendBox.SetLineColor(kWhite);
  TLegend dailyLegend(0.29, 0.20, 0.62, 0.28);
  dailyLegend.SetBorderSize(0);
  dailyLegend.SetFillColor(kWhite);
  dailyLegend.SetFillStyle(1001);
  dailyLegend.SetTextSize(0.045);
  dailyLegend.AddEntry(&dailyLegendBox, "Daily temperature range (min--max)", "f");
  dailyLegend.Draw();
  canvas.Print(outputPdf);

  // Remaining pages: points plus hourly colored min--max boxes for each day.
  for (int day = firstDay; day <= lastDay; ++day) {
    const auto found = dailyPoints.find(day);
    if (found == dailyPoints.end()) continue;
    std::vector<double> time, temperature;
    for (const auto &point : found->second) {
      time.push_back(point.first);
      temperature.push_back(point.second);
    }
    canvas.Clear();
    canvas.SetGrid();
    TGraph graph(time.size(), time.data(), temperature.data());
    graph.SetTitle(Form("logger1 second temperature channel, 2025-11-%02d;Time (JST);Temperature [^{#circ}C]", day));
    graph.SetMarkerStyle(20);
    graph.SetMarkerSize(0.92);
    graph.SetMarkerColor(kBlack);
    graph.Draw("AP");
    graph.GetYaxis()->SetRangeUser(temperatureMin, temperatureMax); graph.GetYaxis()->SetNdivisions(505); graph.GetYaxis()->SetDecimals(kFALSE);
    setTimeAxis(graph, "%H:%M");

    std::vector<TBox> hourlyBoxes;
    hourlyBoxes.reserve(hourlyRange[day].size());
    for (const auto &hourEntry : hourlyRange[day]) {
      const int hour = hourEntry.first;
      const TemperatureRange &range = hourEntry.second;
      const double hourStart = TDatime(2025, 11, day, hour, 0, 0).Convert();
      hourlyBoxes.emplace_back(hourStart, range.min, hourStart + 3600., range.max);
      hourlyBoxes.back().SetFillColorAlpha(kOrange - 2, 0.30);
      hourlyBoxes.back().SetLineColor(kOrange + 7);
      hourlyBoxes.back().Draw("same");
    }
    graph.Draw("P same");
    TBox hourlyLegendBox(0., 0., 1., 1.);
    hourlyLegendBox.SetFillColorAlpha(kOrange - 2, 0.30);
    TLegend hourlyLegend(0.14, 0.73, 0.48, 0.87);
    hourlyLegend.SetBorderSize(0);
    hourlyLegend.SetFillStyle(0);
    hourlyLegend.AddEntry(&graph, "10-minute mean", "p");
    hourlyLegend.AddEntry(&hourlyLegendBox, "Hourly min--max range", "f");
    hourlyLegend.Draw();
    canvas.Print(outputPdf);
  }
  canvas.Print((std::string(outputPdf) + "]").c_str());
  TCanvas pedestalCanvas("pedestalCanvas", "BAC pedestal", 1200, 800);
  pedestalCanvas.Print((std::string(pedestalPdf) + "[").c_str());
  TLatex pedestalInfo; pedestalInfo.SetNDC(); pedestalInfo.SetTextFont(132); pedestalInfo.SetTextSize(0.040);
  pedestalInfo.DrawLatex(0.10, 0.84, "BAC pedestal stability during pedestal-run period");
  pedestalInfo.SetTextSize(0.028);
  pedestalInfo.DrawLatex(0.10, 0.72, "Macro: plot_logger1_temeparue.cc");
  pedestalInfo.DrawLatex(0.10, 0.65, "Input: BAC_ADC_seg4U (= bac_adc_u[4]), Gaussian fit");
  pedestalInfo.DrawLatex(0.10, 0.58, "Temperature: logger1 ch04, 10-minute means; pedestal peak/width in NPE");
  pedestalInfo.DrawLatex(0.10, 0.51, Form("Completed pedestal fits: %zu; panels: up to 16 runs/page", pedestals.size()));
  pedestalCanvas.Print(pedestalPdf);
  AppendPedestalPages(pedestalCanvas, pedestalPdf, allBeamDay, allTemperature, pedestals, firstDay, lastDay, npeFactor, hodoscopeDirectory, pedestalFitHalfWidth, maxPedestalPeakErrorNpe, maxPedestalWidthErrorNpe);
  pedestalCanvas.Print((std::string(pedestalPdf) + "]").c_str());
}
