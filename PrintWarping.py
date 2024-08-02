
#!/usr/bin/python

from ROOT import *
from ROOT import PlotUtils
try:
  import array as np
except:
  exit()
gROOT.SetBatch() #Don't render histograms to a window.  Also gets filled areas correct.

TH1.AddDirectory(False)
#variables = ["pzmu_vs_ptmu", "tpi_vs_thetapi_deg", "pmu_vs_tpi", "tpi_vs_ptmu", "enu_vs_tpi"]
variables = ["tpi_vs_ptmu"]
#variables = ["pzmu_vs_ptmu", "ptmu_vs_pzmu", "tpi_vs_pmu", "pmu_vs_tpi", "tpi_vs_thetapi_deg", "thetapi_deg_vs_tpi", "ptmu_vs_tpi", "tpi_vs_ptmu", "enu_vs_tpi", "tpi_vs_enu"]
date = "20240619"
warp = "WARP3"
plist="ALL"
corrfac="8.8"
scale="4.689984"
for var in variables:
  mcFile = TFile.Open("/minerva/data/users/granados/WarpingStudies/Warping/2DWarping/Warping_2D_{PLIST}_{DATE}_stat_scale{SCALE}_{WARP}_{VAR}_corfac{CORFAC}_exclude.root".format(PLIST=plist,DATE=date,WARP=warp,VAR=var,SCALE=scale,CORFAC=corrfac))
#  mcFile = TFile.Open("/minerva/data/users/granados/WarpingStudies/Warping/2DWarping/Warping_2D_{PLIST}_{DATE}_stat_scale{SCALE}_{WARP}_{VAR}.root".format(PLIST=plist,DATE=date,WARP=warp,VAR=var,SCALE=scale,CORFAC=corrfac))
#  mcFile = TFile.Open("/minerva/data/users/granados/WarpingStudies/Warping/2DWarping/Warping_2D_{PLIST}_{DATE}_{WARP}_{VAR}.root".format(PLIST=plist,DATE=date,WARP=warp,VAR=var,SCALE=scale,CORFAC=corrfac))

  lineWidth = 3

  chi2Iter = mcFile.Get("Chi2_Iteration_Dists/h_chi2_modelData_trueData_iter_chi2") 
  AverageChi2Iter = mcFile.Get("Chi2_Iteration_Dists/m_avg_chi2_modelData_trueData_iter_chi2") 
  truncatedChi2Iter = mcFile.Get("Chi2_Iteration_Dists/m_avg_chi2_modelData_trueData_iter_chi2_truncated") 
  MedianChi2Iter = mcFile.Get("Chi2_Iteration_Dists/h_median_chi2_modelData_trueData_iter_chi2")
  NomEffDenTrue = mcFile.Get("Input_Hists/h_mc_truth")

  nxbins = NomEffDenTrue.GetXaxis().GetNbins()
  nybins = NomEffDenTrue.GetYaxis().GetNbins()
  ndf = double(nybins) * double(nxbins)

  h_ndf = MedianChi2Iter.Clone()
  nbins_ndf = h_ndf.GetXaxis().GetNbins()
  for i in range(nbins_ndf):
    h_ndf.SetBinContent(i+1, double(ndf))

  titleName = ""
  if var == "ptmu_vs_tpi":
    titleName = "p^{t}_{#mu} vs T_{#pi}"
  if var == "pzmu_vs_ptmu":
    titleName = "p^{z}_{#mu} vs p^{t}_{#mu}"
  if var == "tpi_vs_thetapi_deg":
    titleName = "T_{#pi} vs #theta_{#pi}"
  if var == "tpi_vs_pmu":
    titleName = "T_{#pi} vs p_{#mu}"
  if var == "tpi_vs_ptmu":
    titleName = "T_{#pi} vs p^{t}_{#mu}"
  if var == "ptmu_vs_pzmu":
    titleName = "p^{T}_{#mu} vs p^{z}_{#mu}"
  if var == "thetapi_deg_vs_tpi":
    titleName = "#theta_{#pi} vs T_{#pi}"
  if var == "pmu_vs_tpi":
    titleName = "p_{#mu} vs T_{#pi}"
  if var == "enu_vs_tpi":
    titleName = "E_{#nu} vs T_{#pi}"
  if var == "tpi_vs_enu":
    titleName = "T_{#pi} vs E_{#nu}"

  c1 = TCanvas("Warping studies for")
  Title = TPaveText (8., 6000, 120., 10500.)
  Title.Clear()
  Title.AddText(titleName + " " + warp)
  Title.SetShadowColor(0)
  Title.SetLineColor(0)
  Title.SetFillColor(0)
  legend = TLegend(0.7,0.75,0.9,0.9)

  AverageChi2Iter.SetLineWidth(lineWidth)
  AverageChi2Iter.SetLineColor(kRed-4)

  truncatedChi2Iter.SetLineWidth(lineWidth)
  truncatedChi2Iter.SetLineColor(kMagenta-7)

  MedianChi2Iter.SetLineWidth(lineWidth)
  MedianChi2Iter.SetLineColor(kCyan-7)

  h_ndf.SetLineWidth(lineWidth)
  h_ndf.SetLineStyle(10)
  h_ndf.SetLineColor(kOrange-3)

#  c1.UseCurrentStyle()
  chi2Iter.GetYaxis().SetRangeUser(0.,100.)
  chi2Iter.GetXaxis().CenterTitle()
  chi2Iter.GetXaxis().SetTitleOffset(1.3)
  chi2Iter.GetXaxis().SetTitleSize(0.04)
  
  chi2Iter.GetYaxis().CenterTitle()
  chi2Iter.GetYaxis().SetTitleOffset(1.3)
  chi2Iter.GetYaxis().SetTitleSize(0.04)

  gStyle.SetPalette(kAvocado);
  chi2Iter.Draw("colz")
  gStyle.SetOptStat(0)
#  c1.UseCurrentStyle()
  c1.Update()

  AverageChi2Iter.Draw("SAME HIST")
  MedianChi2Iter.Draw("SAME HIST")
  truncatedChi2Iter.Draw("SAME HIST")
  h_ndf.Draw("SAME HIST")
  c1.SetLogx()
  
#  c1.SetLogy()
  #c1.BuildLegend(0.7, 0.6, 0.9, 0.9)
  legend.AddEntry(AverageChi2Iter, "Average", "l")
  legend.AddEntry(truncatedChi2Iter, "Truncated", "l")
  legend.AddEntry(MedianChi2Iter, "Median", "l")
  legend.AddEntry(h_ndf, "ndf", "l")
  legend.Draw()
  Title.Draw()
  c1.Print("WarpingPlots/Warping_{DATE}_{PLIST}_{VAR}_{WARP}_corfac{CORFAC}_exclude.png".format(DATE=date,VAR=var,WARP=warp,PLIST=plist,CORFAC=corrfac))
#  c1.Print("WarpingPlots/Warping_{DATE}_{PLIST}_{VAR}_{WARP}.png".format(DATE=date,VAR=var,WARP=warp,PLIST=plist,CORFAC=corrfac))
  c1.Clear()

