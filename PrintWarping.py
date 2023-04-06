
#!/usr/bin/python

from ROOT import *
from ROOT import PlotUtils
try:
  import array as np
except:
  exit()
gROOT.SetBatch() #Don't render histograms to a window.  Also gets filled areas correct.

TH1.AddDirectory(False)
mcFile = TFile.Open("/minerva/data/users/granados/WarpingStudies/Warping/2DWarping/AaronBinning/Warping_2DWARP1_pzmu_vs_ptmu.root")

lineWidth = 2

chi2Iter = mcFile.Get("Chi2_Iteration_Dists/h_chi2_modelData_trueData_iter_chi2") 
AverageChi2Iter = mcFile.Get("Chi2_Iteration_Dists/m_avg_chi2_modelData_trueData_iter_chi2") 
truncatedChi2Iter = mcFile.Get("Chi2_Iteration_Dists/m_avg_chi2_modelData_trueData_iter_chi2_truncated") 
MedianChi2Iter = mcFile.Get("Chi2_Iteration_Dists/h_median_chi2_modelData_trueData_iter_chi2")

c1 = TCanvas("Warping studies for")

legend = TLegend(0.75,0.8,0.9,0.9)

AverageChi2Iter.SetLineWidth(lineWidth)
AverageChi2Iter.SetLineColor(kRed-4)

truncatedChi2Iter.SetLineWidth(lineWidth)
truncatedChi2Iter.SetLineColor(kMagenta-7)

MedianChi2Iter.SetLineWidth(lineWidth)
MedianChi2Iter.SetLineColor(kCyan-7)

gStyle.SetPalette(kAvocado);
chi2Iter.Draw("colz")
c1.UseCurrentStyle()

c1.Update()
gStyle.SetOptStat(0)

AverageChi2Iter.Draw("SAME HIST")
MedianChi2Iter.Draw("SAME HIST")
truncatedChi2Iter.Draw("SAME HIST")
c1.SetLogx()
c1.SetLogy()
#c1.BuildLegend(0.7, 0.6, 0.9, 0.9)
legend.AddEntry(AverageChi2Iter, "Average", "l")
legend.AddEntry(truncatedChi2Iter, "Truncated", "l")
legend.AddEntry(MedianChi2Iter, "Median", "l")
legend.Draw()
c1.Print("WarpingTest.png")
c1.Clear()

