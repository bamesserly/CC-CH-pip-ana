#!/usr/bin/python

from ROOT import *
from ROOT import PlotUtils

try:
    import array as np
except:
    exit()
gROOT.SetBatch()  # Don't render histograms to a window.  Also gets filled areas correct.

TH1.AddDirectory(False)
variables = [
    "mixtpi",
    "mixthetapi_deg",
    "mtpi",
    "enu",
    "pmu",
    "ptmu",
    "pzmu",
    "q2",
    "thetamu_deg",
    "wexp",
]
# variables = ["mixtpi"]
date = "20230523"
warp = "WARP2"

for var in variables:
    mcFile = TFile.Open(
        "WarpingStudies/WexpModBinning/Warping_{WARP}_{VAR}.root".format(
            WARP=warp, VAR=var
        )
    )

    lineWidth = 3

    chi2Iter = mcFile.Get("Chi2_Iteration_Dists/h_chi2_modelData_trueData_iter_chi2")
    AverageChi2Iter = mcFile.Get(
        "Chi2_Iteration_Dists/m_avg_chi2_modelData_trueData_iter_chi2"
    )
    truncatedChi2Iter = mcFile.Get(
        "Chi2_Iteration_Dists/m_avg_chi2_modelData_trueData_iter_chi2_truncated"
    )
    MedianChi2Iter = mcFile.Get(
        "Chi2_Iteration_Dists/h_median_chi2_modelData_trueData_iter_chi2"
    )
    NomEffDenTrue = mcFile.Get("Input_Hists/h_mc_truth")

    nxbins = NomEffDenTrue.GetXaxis().GetNbins()
    nybins = NomEffDenTrue.GetYaxis().GetNbins()
    ndf = double(nybins) * double(nxbins)

    h_ndf = MedianChi2Iter.Clone()
    nbins_ndf = h_ndf.GetXaxis().GetNbins()
    for i in range(nbins_ndf):
        h_ndf.SetBinContent(i + 1, double(ndf))

    titleName = ""
    if var == "mixtpi":
        titleName = "T_{#pi}"
    if var == "pzmu":
        titleName = "p^{z}_{#mu}"
    if var == "ptmu":
        titleName = "p^{t}_{#mu}"
    if var == "thetapi_deg":
        titleName = "#theta_{#pi}"
    if var == "pmu":
        titleName = "p_{#mu}"
    if var == "q2":
        titleName = "q^{2}"
    if var == "thetamu_deg":
        titleName = "#theta_{#mu}"
    if var == "wexp":
        titleName = "W_{exp}"
    if var == "mixthetapi_deg":
        titleName = "#theta_{#pi}"
    if var == "enu":
        titleName = "E_{#nu}"
    c1 = TCanvas("Warping studies for")
    Title = TPaveText(8.0, 6000, 120.0, 10500.0)
    Title.Clear()
    Title.AddText(titleName + " " + warp)
    Title.SetShadowColor(0)
    Title.SetLineColor(0)
    Title.SetFillColor(0)
    legend = TLegend(0.7, 0.75, 0.9, 0.9)

    AverageChi2Iter.SetLineWidth(lineWidth)
    AverageChi2Iter.SetLineColor(kRed - 4)

    truncatedChi2Iter.SetLineWidth(lineWidth)
    truncatedChi2Iter.SetLineColor(kMagenta - 7)

    MedianChi2Iter.SetLineWidth(lineWidth)
    MedianChi2Iter.SetLineColor(kCyan - 7)

    h_ndf.SetLineWidth(lineWidth)
    h_ndf.SetLineStyle(10)
    h_ndf.SetLineColor(kOrange - 3)

    #  c1.UseCurrentStyle()

    chi2Iter.GetXaxis().CenterTitle()
    chi2Iter.GetXaxis().SetTitleOffset(1.3)
    chi2Iter.GetXaxis().SetTitleSize(0.04)

    chi2Iter.GetYaxis().CenterTitle()
    chi2Iter.GetYaxis().SetTitleOffset(1.3)
    chi2Iter.GetYaxis().SetTitleSize(0.04)

    gStyle.SetPalette(kAvocado)
    chi2Iter.Draw("colz")
    gStyle.SetOptStat(0)
    #  c1.UseCurrentStyle()
    c1.Update()

    AverageChi2Iter.Draw("SAME HIST")
    MedianChi2Iter.Draw("SAME HIST")
    truncatedChi2Iter.Draw("SAME HIST")
    h_ndf.Draw("SAME HIST")
    c1.SetLogx()
    c1.SetLogy()
    # c1.BuildLegend(0.7, 0.6, 0.9, 0.9)
    legend.AddEntry(AverageChi2Iter, "Average", "l")
    legend.AddEntry(truncatedChi2Iter, "Truncated", "l")
    legend.AddEntry(MedianChi2Iter, "Median", "l")
    legend.AddEntry(h_ndf, "ndf", "l")
    legend.Draw()
    Title.Draw()
    c1.Print(
        "WarpingStudies/WexpModBinning/Warping_{VAR}_{WARP}.png".format(
            VAR=var, WARP=warp
        )
    )
    c1.Clear()
