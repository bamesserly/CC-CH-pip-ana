
#!/usr/bin/python

from ROOT import *
from ROOT import PlotUtils
from array import array
try:
  import array as np
except:
  exit()
gROOT.SetBatch() #Don't render histograms to a window.  Also gets filled areas correct.

TH1.AddDirectory(False)
lowq2MParamFiles = TFile.Open("../opt/etc/MParamFiles/data/Reweight/lowQ2pi_weights.root")
MyLowq2 = TFile.Open("../NewLowQ2/lowQ2pi_weights.root")

g_lowq2MP = lowq2MParamFiles.Get("MENU1PI_weights")
g_MyLowq2 = MyLowq2.Get("MENU1PI_weights")

nPMP = g_lowq2MP.GetN()
nPMy = g_MyLowq2.GetN()

print (nPMP)
print (nPMy)

rval = array ('d')
xval = array ('d')

xMp = double(0)
yMp = double(0)
xMy = double(0)
yMy = double(0)

for p in range(nPMP):
  g_lowq2MP.GetPoint(p, xMp, yMp)
  g_MyLowq2.GetPoint(p, xMy, yMy)
  rval.append( yMp / yMy) 
  xval.append(xMp)

for i in range(10):
  print (rval[i])
  print (xval[i])

c1 = TCanvas("Warping studies for")
ratio = TGraph(nPMP, xval, rval) 
c1.cd()
ratio.GetYaxis().SetRangeUser(0.9,1.1)
ratio.GetYaxis().SetTitle("ratio")
ratio.GetXaxis().SetTitle("q^2 GeV")
ratio.Draw()
c1.Update()
c1.Print("RatioLowq2.png")






