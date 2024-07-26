

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
    "enu",
    "pmu",
    "ptmu",
    "pzmu",
    "q2",
    "thetamu_deg",
    "wexp",
]
# variables = ["mixtpi"]
date = "20240424"
warp = "WARP7"
plist = "ALL"
lineWidth = 3
fprev = TFile.Open(
        "DataXSecInputs_20240424_me1A_mixed_NewEstimatorptmucut_noSys_p4_NOMINAL.root")

fval = TFile.Open(
        "DataXSecInputs_20240506_me1A_mixed_nosys_validationp6.root")


for var in variables:
    xsecdata = fprev.Get("cross_section_{V}".format(V=var))
    xsecdataval = fval.Get("cross_section_{V}".format(V=var))
    xsecmc = fprev.Get("mc_cross_section_{V}".format(V=var))
    xsecmcval = fval.Get("mc_cross_section_{V}".format(V=var))

    nxbins = xsecdata.GetXaxis().GetNbins()
    ratio = xsecdataval.Clone()
    ratio1 = xsecmcval.Clone()
    oneline = xsecdata.Clone()
    for i in range(nxbins):
        oneline.SetBinContent(i + 1, 1)

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
    
    c2 = TCanvas("validation data")

    Title1 = TPaveText(8.0, 6000, 120.0, 10500.0)
    Title1.Clear()
    Title1.AddText("Ratio Nominal vs " + warp + " " + titleName)
    Title1.SetShadowColor(0)
    Title1.SetLineColor(0)
    Title1.SetFillColor(0)
    legend1 = TLegend(0.7, 0.75, 0.9, 0.9)

    ratio.Divide(ratio, xsecdataval) 
    ratio.SetLineWidth(lineWidth)
    ratio.SetLineColor(kRed - 4)
    c2.UseCurrentStyle()

    oneline.SetLineWidth(lineWidth)
    oneline.SetLineStyle(2)
    oneline.SetLineColor(kGray)

    oneline.GetXaxis().CenterTitle()
    oneline.GetXaxis().SetTitleOffset(1.3)
    oneline.GetXaxis().SetTitleSize(0.04)

    oneline.GetYaxis().CenterTitle()
    oneline.GetYaxis().SetTitleOffset(1.3)
    oneline.GetYaxis().SetTitleSize(0.04)

    oneline.Draw("HIST")
    gStyle.SetOptStat(0)
    #c2.Update()
    #ratio.SetMinimum(0.9);
    #ratio.GetYaxis().SetRangeUser(0.9,1.1);
    ratio.Draw("SAME HIST")
    #c1.SetLogx()
    #c1.SetLogy()
    # c1.BuildLegend(0.7, 0.6, 0.9, 0.9)
    legend1.AddEntry(ratio, "prepp6/P4", "l")
    legend1.Draw()
    #Title1.Draw()
    c2.Print(
        "Validation/validation_data_{VAR}_.png".format(
            VAR=var)
    )
    c3 = TCanvas("validation mc")

    Title = TPaveText(8.0, 6000, 120.0, 10500.0)
    Title.Clear()
    Title.AddText("Ratio Nominal vs " + warp + " " + titleName)
    Title.SetShadowColor(0)
    Title.SetLineColor(0)
    Title.SetFillColor(0)
    legend = TLegend(0.7, 0.75, 0.9, 0.9)

    ratio1.Divide(ratio1, xsecmc) 
    ratio1.SetLineWidth(lineWidth)
    ratio1.SetLineColor(kRed - 4)
    c3.UseCurrentStyle()

    oneline.SetLineWidth(lineWidth)
    oneline.SetLineStyle(2)
    oneline.SetLineColor(kGray)

    oneline.GetXaxis().CenterTitle()
    oneline.GetXaxis().SetTitleOffset(1.3)
    oneline.GetXaxis().SetTitleSize(0.04)

    oneline.GetYaxis().CenterTitle()
    oneline.GetYaxis().SetTitleOffset(1.3)
    oneline.GetYaxis().SetTitleSize(0.04)

    oneline.Draw("HIST")
    gStyle.SetOptStat(0)
    #c2.Update()

    #ratio1.GetYaxis().SetRangeUser(0.9, 1.1);
    ratio1.Draw("SAME HIST")
    #c1.SetLogx()
    #c1.SetLogy()
    # c1.BuildLegend(0.7, 0.6, 0.9, 0.9)
    legend.AddEntry(ratio1, "prep6/P4", "l")
    legend.Draw()
    #Title1.Draw()
    c3.Print(
        "Validation/validation_mc_{VAR}_.png".format(
            VAR=var
        )
    )

    c2.Clear()
    c3.Clear()


