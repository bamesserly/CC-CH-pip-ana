# Fit a line to tpi vs range and plot it.
# Script from Mehreen, modernized and tweaked by Ben
#
# The latest and greatest input file to the best of my knowledge is:
# /exp/minerva/data/users/sultana/March2024_Analysis_Code/ \
#   July182024_sl7_newstudies_2Ddistcut_AvailECut12_sysOn/ \
#     summedandscaled/scaled_MC.root
from ROOT import (
    PlotUtils,
    TH1,
    TFile,
    TCanvas,
    TGraph,
    TGraphErrors,
    TF1,
    kFALSE,
    TMath,
)
import os
import sys
import array as ar

TH1.AddDirectory(kFALSE)

mcfile = TFile.Open(sys.argv[1])

FIT_FORM = "[0]*x + [1]*sqrt(x)"

# Get input 2D hist

# N.B. the x-axis LABEL of this plot is WRONG: it's tpi not range
range_tpi_hist = mcfile.Get("Pion_Range_vs_Pion_KE_signal")
range_tpi_hist.SetDirectory(0)

# re-use the hist name
range_tpi_hist = range_tpi_hist.GetCVHistoWithStatError()
n_tpi_bins = range_tpi_hist.GetNbinsX() + 1
n_range_bins = range_tpi_hist.GetNbinsY() + 1
print(
    "There are ", n_tpi_bins, " NXBins (PionKE) and ", n_range_bins, " NYBins (Range)"
)

# ==============================================================================
# Attempt to get tpi median for each bin of range using two different methods.
#
# (1) ProjectionX + GetQuantiles. Allows the median to fall inbetween some
# fraction of the way through a bin.
#
# (2) direct access of bin contents + TMath. "Rounds" the median to a tpi bin
# center.
#
# M went with the second method, though there was a bug in her implementation
# of the first. I don't have any particular objection to using the second
# method.
# ==============================================================================

# ==============================================================================
# (1) GetProjectionX + GetQuantiles. Median can fall part way through a bin.
#
# We don't currently use this method.
#
# Fixed a bug with how ProjectionX was being called.
#
# Median in detail:
# - Find the bin b inside of which 50% of the integral lies above and below.
#   In other words:
#   let x =0.5*(integral from 0 to N)
#   find first b such that (integral from 0 to b) >= x
# - Then, linearly interpolate:
#   p = (x - hist.GetBinLowEdge(b)) / hist.GetBinContent(b)
# - The median lies p% of the way through bin b.
#   median = hist.GetBinWidth(b) * p + hist.GetBinLowEdge(b)
# ==============================================================================
for range_bin in range(1, n_range_bins):
    # range bin center
    range_val = range_tpi_hist.GetYaxis().GetBinCenter(range_bin)

    # tpi hist for a given range bin
    tpihist_bini = range_tpi_hist.ProjectionX(
        f"tpi_slice{range_bin}", range_bin, range_bin, "e"
    )

    # compute median
    median_definition = ar.array("d", [0.5])  # median definition
    median_tpi_array = ar.array("d", [0.0])  # where to store result
    z = tpihist_bini.GetQuantiles(1, median_tpi_array, median_definition)
    median_tpi = median_tpi_array[0]

    print(f"range bin center {range_val}, median tpi {median_tpi}")

# ==============================================================================
# (2) Print the median Tpi for each bin of range getting bin contents manually
# and using TMath. The median rounds to a single bin center value.
# Median method: https://root.cern.ch/root/html524/TMath.html#TMath:Median
# Save to file.
# ===============================================================================
range_bin_centers = ar.array("d", [0] * n_range_bins)
tpi_medians = ar.array("d", [0] * n_range_bins)
tpierrs = ar.array("d", [0] * n_range_bins)
rangeerrs = ar.array("d", [0] * n_range_bins)

for range_bin in range(1, n_range_bins):
    range_bin_center = range_tpi_hist.GetYaxis().GetBinCenter(range_bin)

    tpi_bin_centers = ar.array("d", [0] * n_tpi_bins)
    z = ar.array("d", [0] * n_tpi_bins)
    # get array of 2D hist contents
    for tpi_bin in range(1, n_tpi_bins):
        tpi_bin_centers[tpi_bin] = range_tpi_hist.GetXaxis().GetBinCenter(tpi_bin)
        z[tpi_bin] = range_tpi_hist.GetBinContent(tpi_bin, range_bin)
    sumz = sum(z)
    median_tpi = TMath.Median(n_tpi_bins, tpi_bin_centers, z)
    meanz = TMath.RMS(n_tpi_bins, tpi_bin_centers, z)
    tpierrs[range_bin] = meanz / TMath.Sqrt(sumz) if sumz > 0 else 999

    range_bin_centers[range_bin] = range_bin_center
    tpi_medians[range_bin] = median_tpi

    print(f"range {range_bin_center}, median tpi {median_tpi}")

# ==============================================================================
# Plot the medians vs range, fit, and save to a file
# ==============================================================================
outfile = TFile.Open("tpi_from_range_fit.root", "RECREATE")

# medians vs range plot
range_bin_centers = range_bin_centers[1:-1]
tpi_medians = tpi_medians[1:-1]
g = TGraph(len(range_bin_centers), range_bin_centers, tpi_medians)
g.SetName("TpiEst")
g.GetXaxis().SetTitle("Pion Range (mm)")
g.GetYaxis().SetTitle("Pion KE (MeV)")
g.Draw("AC*")
g.Write()

# medians vs range with errors plot
ge = TGraphErrors(
    len(range_bin_centers), range_bin_centers, tpi_medians, rangeerrs, tpierrs
)
ge.SetName("TpiEst_Errs")
ge.GetXaxis().SetTitle("Pion Range (mm)")
ge.GetYaxis().SetTitle("Pion KE (MeV)")
ge.Draw("ALP")
ge.Write()

# do the fit
overall = TCanvas("Fit")
func = TF1("func", FIT_FORM, 0.0, 350.0)
# func2 = TF1('func2', 'gaus',350., 3000.)
# func3 = TF1('func3', '[0]*x + [1]*sqrt(x)+gaus(3)', 0.0, 3000.)
# func = TF1('func', 'landau', 0.0, 3000.)
gefit = ge.Clone()
gefit.SetName("TpiEst_Errs_Fit")
fit = gefit.Fit("func", "S")
# gefit.Fit('func2','R+')

# par = ar.array('d',[]*2)

# func.GetParameters(par[1])
# func2.GetParameters(par[3])

# func3.SetParameters(par)
# fit = gefit.Fit(func3, 'R+')

rangepoints = ar.array("d", [0] * 2400)
tpipoints = ar.array("d", [0] * 2400)
for q in range(0, 2400):
    rangepoints[q] = float(q)
    tpipoints[q] = func.Eval(q)

# draw the fit and save it to file
gefit.Draw("ALP")
overall.Print("fitResult_1.root")
overall.Print("fitResult_1.png")
chiperndf = fit.Chi2() / fit.Ndf()

# Chi 2s
print("Chi Squared: ", fit.Chi2(), " NdF: ", fit.Ndf())
print("Chi Squared Per Degree of Freedom: ", chiperndf)
gefit.Write()
# fitline = gefit.GetCurve()
# fitline.Draw()

# Also TF1, sure
curve = TGraph(len(rangepoints), rangepoints, tpipoints)
# fit.Scan(curve, 0.0, 2400.)
curve.SetName("FitCurve")
curve.GetXaxis().SetTitle("Pion Range (mm)")
curve.GetYaxis().SetTitle("Pion KE (MeV)")
curve.Draw("AC*")
curve.Write()

# ==============================================================================
# finally, save the input hists to the output file as well
# ==============================================================================
for key in mcfile.GetListOfKeys():
    obj = key.ReadObj()

    if "Pion_Range_vs_Pion" in obj.GetName():
        obj.Write()

    if "POTUsed" in obj.GetName():
        obj.Write()

outfile.Close()
mcfile.Close()
