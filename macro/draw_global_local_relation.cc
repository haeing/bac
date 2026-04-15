#include <TCanvas.h>
#include <TGraph.h>
#include <TArrow.h>
#include <TBox.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TMath.h>
#include <TVector3.h>
#include <vector>

// --------------------------------------------------
// helper: project 3D point to 2D screen coordinates
// using an orthographic projection with a chosen
// Geant4-like viewpoint (theta, phi)
// --------------------------------------------------
TVector3 ProjectToScreen(const TVector3& P,
                         const TVector3& eH,
                         const TVector3& eV)
{
  return TVector3(P.Dot(eH), P.Dot(eV), 0.0);
}

void draw_global_local_relation()
{
  //gStyle->SetOptStat(0);

  // ==================================================
  // Geometry / transform parameters from your code
  // ==================================================
  const double p  = 36.0;     // mm
  const double x0 = 0.0;      // mm
  const double y0 = 117.5;    // mm
  const double z0 = 84.0;     // mm

  // aerogel for reference (schematic box centered at global origin)
  const double aero_thick  = 37.0;   // mm
  const double aero_height = 115.0;  // mm

  // ==================================================
  // Viewer direction: use same as Geant4 command
  // /vis/viewer/set/viewpointThetaPhi 90 150
  // ==================================================
  const double theta = 90.0 * TMath::DegToRad();
  const double phi   = 150.0 * TMath::DegToRad();

  // View direction
  TVector3 viewDir(std::sin(theta)*std::cos(phi),
                   std::sin(theta)*std::sin(phi),
                   std::cos(theta));
  viewDir = viewDir.Unit();

  // Choose an "up" reference not parallel to viewDir
  TVector3 upRef(0,1,0);
  if (std::abs(viewDir.Dot(upRef)) > 0.99) upRef = TVector3(0,0,1);

  // Screen horizontal / vertical basis
  TVector3 eH = upRef.Cross(viewDir).Unit();   // to the right on screen
  TVector3 eV = viewDir.Cross(eH).Unit();      // upward on screen

  // ==================================================
  // Define local basis of reflector
  //
  // local parabola from your code:
  // y_loc = -x_loc^2 / (4p)
  //
  // define:
  //   v-axis = +x_loc   (transverse)
  //   u-axis = -y_loc   (symmetry axis)
  //   w-axis = +z_loc   (extrusion direction)
  // ==================================================
  TVector3 vLoc(1,0,0);
  TVector3 uLoc(0,-1,0);
  TVector3 wLoc(0,0,1);

  // Apply Geant4 rotation: first Y, then Z
  auto ApplyGeant4Rotation = [](TVector3 vec) {
    vec.RotateY(+90.0 * TMath::DegToRad());
    vec.RotateZ(+65.0 * TMath::DegToRad());
    return vec;
  };

  TVector3 vG = ApplyGeant4Rotation(vLoc).Unit();
  TVector3 uG = ApplyGeant4Rotation(uLoc).Unit();
  TVector3 wG = ApplyGeant4Rotation(wLoc).Unit();

  // Reflector local origin in global coordinates
  TVector3 R0(x0,y0,z0);

  // Global origin (aerogel center)
  TVector3 G0(0,0,0);

  // ==================================================
  // Canvas + frame in projected coordinates
  // ==================================================
  TCanvas* c = new TCanvas("c_proj","Projected local/global axes",900,900);
  c->SetMargin(0.10,0.05,0.10,0.05);

  // Project a few reference points to estimate frame range
  std::vector<TVector3> testPts;

  // aerogel corners
  testPts.push_back(TVector3(-aero_thick/2., -aero_height/2., 0));
  testPts.push_back(TVector3(+aero_thick/2., -aero_height/2., 0));
  testPts.push_back(TVector3(+aero_thick/2., +aero_height/2., 0));
  testPts.push_back(TVector3(-aero_thick/2., +aero_height/2., 0));

  // reflector origin and local-axis endpoints
  testPts.push_back(R0);
  testPts.push_back(R0 + 60*uG);
  testPts.push_back(R0 + 60*vG);

  double xmin=+1e9, xmax=-1e9, ymin=+1e9, ymax=-1e9;
  for (auto &P : testPts) {
    TVector3 q = ProjectToScreen(P, eH, eV);
    xmin = std::min(xmin, q.X());
    xmax = std::max(xmax, q.X());
    ymin = std::min(ymin, q.Y());
    ymax = std::max(ymax, q.Y());
  }

  double dx = xmax-xmin;
  double dy = ymax-ymin;
  xmin -= 0.25*dx; xmax += 0.25*dx;
  ymin -= 0.25*dy; ymax += 0.25*dy;

  TH2D* frame = new TH2D("frame","",100,xmin,xmax,100,ymin,ymax);
  frame->SetStats(0);
  frame->GetXaxis()->SetTitle("projected horizontal coordinate");
  frame->GetYaxis()->SetTitle("projected vertical coordinate");
  frame->Draw();

  TLatex latex;
  latex.SetTextFont(42);

  // ==================================================
  // Draw projected aerogel box
  // ==================================================
  TGraph* gAer = new TGraph(5);
  {
    TVector3 c1 = ProjectToScreen(TVector3(-aero_thick/2., -aero_height/2., 0), eH, eV);
    TVector3 c2 = ProjectToScreen(TVector3(+aero_thick/2., -aero_height/2., 0), eH, eV);
    TVector3 c3 = ProjectToScreen(TVector3(+aero_thick/2., +aero_height/2., 0), eH, eV);
    TVector3 c4 = ProjectToScreen(TVector3(-aero_thick/2., +aero_height/2., 0), eH, eV);
    gAer->SetPoint(0,c1.X(),c1.Y());
    gAer->SetPoint(1,c2.X(),c2.Y());
    gAer->SetPoint(2,c3.X(),c3.Y());
    gAer->SetPoint(3,c4.X(),c4.Y());
    gAer->SetPoint(4,c1.X(),c1.Y());
  }
  gAer->SetLineColor(kBlack);
  gAer->SetLineWidth(2);
  gAer->Draw("L SAME");

  // ==================================================
  // Draw projected global axes at aerogel center
  // ==================================================
  TVector3 qG0 = ProjectToScreen(G0, eH, eV);
  TVector3 qGY = ProjectToScreen(G0 + 40*TVector3(0,1,0), eH, eV);
  TVector3 qGZ = ProjectToScreen(G0 + 40*TVector3(0,0,1), eH, eV);

  TArrow* aY = new TArrow(qG0.X(), qG0.Y(), qGY.X(), qGY.Y(), 0.02, "|>");
  aY->SetLineColor(kGreen+2);
  aY->SetFillColor(kGreen+2);
  aY->SetLineWidth(3);
  aY->Draw();

  TArrow* aZ = new TArrow(qG0.X(), qG0.Y(), qGZ.X(), qGZ.Y(), 0.02, "|>");
  aZ->SetLineColor(kBlue+1);
  aZ->SetFillColor(kBlue+1);
  aZ->SetLineWidth(3);
  aZ->Draw();

  latex.SetTextSize(0.035);
  latex.SetTextColor(kGreen+2);
  latex.DrawLatex(qGY.X()+2, qGY.Y()+2, "y");

  latex.SetTextColor(kBlue+1);
  latex.DrawLatex(qGZ.X()+2, qGZ.Y()+2, "z");

  latex.SetTextColor(kBlack);
  latex.DrawLatex(qG0.X()+2, qG0.Y()-5, "global origin");
  latex.DrawLatex(qG0.X()+2, qG0.Y()-10, "(aerogel center)");

  // ==================================================
  // Draw projected local axes at reflector reference point
  // ==================================================
  TVector3 qR0 = ProjectToScreen(R0, eH, eV);
  TVector3 qRU = ProjectToScreen(R0 + 45*uG, eH, eV);
  TVector3 qRV = ProjectToScreen(R0 + 45*vG, eH, eV);

  TArrow* aU = new TArrow(qR0.X(), qR0.Y(), qRU.X(), qRU.Y(), 0.02, "|>");
  aU->SetLineColor(kRed+2);
  aU->SetFillColor(kRed+2);
  aU->SetLineWidth(3);
  aU->Draw();

  TArrow* aV = new TArrow(qR0.X(), qR0.Y(), qRV.X(), qRV.Y(), 0.02, "|>");
  aV->SetLineColor(kMagenta+1);
  aV->SetFillColor(kMagenta+1);
  aV->SetLineWidth(3);
  aV->Draw();

  latex.SetTextColor(kRed+2);
  latex.DrawLatex(qRU.X()+2, qRU.Y()+2, "u");

  latex.SetTextColor(kMagenta+1);
  latex.DrawLatex(qRV.X()+2, qRV.Y()+2, "v");

  // local origin marker
  TGraph* gLoc = new TGraph(1);
  gLoc->SetPoint(0, qR0.X(), qR0.Y());
  gLoc->SetMarkerStyle(20);
  gLoc->SetMarkerSize(1.2);
  gLoc->SetMarkerColor(kBlack);
  gLoc->Draw("P SAME");

  latex.SetTextColor(kBlack);
  latex.DrawLatex(qR0.X()+2, qR0.Y()+6, "local origin");

  // ==================================================
  // Draw projected parabola
  //
  // local point:
  //   (x_loc, y_loc, z_loc) = (v, -v^2/(4p), 0)
  // ==================================================
  TGraph* gPar = new TGraph();
  const int N = 400;
  for (int i=0; i<N; ++i) {
    double v = -70.0 + 140.0 * i / (N-1.0);
    double xLoc = v;
    double yLoc = -v*v/(4.0*p);
    double zLoc = 0.0;

    TVector3 P(xLoc,yLoc,zLoc);
    P = ApplyGeant4Rotation(P);
    P += R0;

    TVector3 q = ProjectToScreen(P, eH, eV);
    gPar->SetPoint(i, q.X(), q.Y());
  }
  gPar->SetLineColor(kOrange+7);
  gPar->SetLineWidth(3);
  gPar->Draw("L SAME");

  // labels
  latex.SetTextSize(0.030);
  latex.SetTextColor(kBlack);
  latex.DrawLatex(xmin+5, ymax-10, "Projected global/local relation");
  latex.DrawLatex(xmin+5, ymax-18, "Reflector placed with rotateY(90^{#circ}) then rotateZ(65^{#circ})");

  TLegend* leg = new TLegend(0.56,0.72,0.90,0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(gPar, "projected parabolic profile", "l");
  leg->AddEntry(aU, "local axis u", "l");
  leg->AddEntry(aV, "local axis v", "l");
  leg->Draw();

  //c->SaveAs("projected_local_global_axes.pdf");
  //c->SaveAs("projected_local_global_axes.png");
}
