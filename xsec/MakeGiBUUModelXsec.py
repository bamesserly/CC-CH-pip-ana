import ROOT
from ROOT import PlotUtils
import os, subprocess,re, sys,glob
from collections import OrderedDict,namedtuple
from array import array

from SelectGiBUUEvents import *

fileList = sys.argv[1]
outname = sys.argv[2]
enu_edges  = array('d',[0., 1.e3, 3.e3, 4.e3, 6.5e3, 9.5e3, 14.e3, 20.e3])
pmu_edges  = array('d',[0., 1.e3,  2.e3, 3.e3, 4.e3, 5.5e3, 7.5e3, 10.e3, 13.e3, 20.e3])
ptmu_edges = array('d',[ 0., 1.e2, 2.e2, 3.e2, 4.e2, 5.e2, 6.e2, 8.e2, 10.e2, 12.5e2, 15.e2, 25.e2])
thmu_edges = array('d',[ 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 14., 16., 20.])
q2_edges   = array('d',[ 0, 0.025e6, 0.05e6, 0.1e6, 0.2e6, 0.3e6, 0.4e6, 0.5e6, 0.7e6,   1.0e6,  1.3e6, 2.0e6, 3.0e6])
W_edges    = array('d',[ 0., 10.e2, 11.e2, 12.e2, 13.e2, 14.e2, 15.e2])
tpi_edges  = array('d',[ 0.,20.,35.,50.,65.,80.,100., 125., 165., 200., 350. ])
thpi_edges = array('d',[ 0., 15., 30., 45, 60, 76., 108., 122., 136., 150., 165., 180.]) 

enu     = PlotUtils.MnvH1D("enu"    ,";E_{#nu};#frac{d#sigma}{dE_{#nu}}"             ,len(enu_edges)-1,enu_edges);
pmu     = PlotUtils.MnvH1D("pmu"    ,";p_{#mu};#frac{d#sigma}{dp_{#mu}}"             ,len(pmu_edges)-1,pmu_edges);
ptmu    = PlotUtils.MnvH1D("ptmu"   ,";p_{#mu,T};#frac{d#sigma}{dp_{#mu,T}}"         ,len(ptmu_edges)-1,ptmu_edges);
pzmu    = PlotUtils.MnvH1D("pzmu"   ,";p_{#mu,||};#frac{d#sigma}{dp_{#mu,||}}"       ,len(pmu_edges)-1,pmu_edges);
thetamu_deg = PlotUtils.MnvH1D("thmu",";#theta_{#mu};#frac{d#sigma}{d#theta_{#mu}}"   ,len(thmu_edges)-1,thmu_edges);
q2         = PlotUtils.MnvH1D("q2"        ,";Q^{2};#frac{d#sigma}{dQ^{2}}"                 ,len(q2_edges)-1,q2_edges);
wexp          = PlotUtils.MnvH1D("wexp"         ,";W_{exp};#frac{d#sigma}{dW_{exp}}"             ,len(W_edges)-1,W_edges);
tpi  = PlotUtils.MnvH1D("tpi" ,";T_{#pi};#frac{d#sigma}{dT_{#pi}}"             ,len(tpi_edges)-1,tpi_edges);
thetapi_deg = PlotUtils.MnvH1D("thpi",";#theta_{#pi};#frac{d#sigma}{d#theta_{#pi}}"   ,len(thpi_edges)-1,thpi_edges);

enu    .SetDirectory(0) 
pmu    .SetDirectory(0) 
ptmu   .SetDirectory(0)
pzmu   .SetDirectory(0)
thetamu_deg.SetDirectory(0)
q2        .SetDirectory(0)
wexp         .SetDirectory(0)
tpi .SetDirectory(0)
thetapi_deg.SetDirectory(0)

nEvents = 0
nSignal = 0
nFiles = 0
with open(fileList) as fnames:
  for fname in fnames:
    print (fname)
    nFiles+=1
    for event in readEvents(fname.rstrip('\n')):
      nEvents +=1
      bSignal,pievtkine = isSignal(event)
      if not bSignal:
        continue

      nSignal+=1
      enu    .Fill( pievtkine.enu*1000, event.wgt )   
      pmu    .Fill( pievtkine.pmu*1000, event.wgt )   
      ptmu   .Fill( pievtkine.pmuT*1000,event.wgt )
      pzmu   .Fill( pievtkine.pmuZ*1000,event.wgt )
      thetamu_deg.Fill( pievtkine.thmu,     event.wgt )
      q2        .Fill( pievtkine.Q*1000000,event.wgt )
      wexp         .Fill( pievtkine.W*1000,   event.wgt )
      tpi .Fill( pievtkine.tpi*1000, event.wgt )
      thetapi_deg.Fill( pievtkine.thpi,     event.wgt )


enu    .Scale(1./nFiles)  
pmu    .Scale(1./nFiles)  
ptmu   .Scale(1./nFiles)
pzmu   .Scale(1./nFiles)
thetamu_deg.Scale(1./nFiles)
q2        .Scale(1./nFiles)
wexp         .Scale(1./nFiles)
tpi .Scale(1./nFiles)
thetapi_deg.Scale(1./nFiles)

#pmu    .Scale(12.0e-38/nFiles,"width")  
#ptmu   .Scale(12.0e-38/nFiles,"width")
#pzmu   .Scale(12.0e-38/nFiles,"width")
#thetamu_deg.Scale(12.0e-38/nFiles,"width")
#q2        .Scale(12.0e-38/nFiles,"width")
#wexp         .Scale(12.0e-38/nFiles,"width")
#tpi .Scale(12.0e-38/nFiles,"width")
#thetapi_deg.Scale(12.0e-38/nFiles,"width")

output_file = ROOT.TFile(outname,"RECREATE");
print ("Writing output to "+outname)
output_file.cd();
enu    .Write() 
pmu    .Write() 
ptmu   .Write()
pzmu   .Write()
thetamu_deg.Write()
q2        .Write()
wexp         .Write()
tpi .Write()
thetapi_deg.Write()

output_file.Close();

