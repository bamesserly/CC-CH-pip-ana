import ROOT
from ROOT import PlotUtils
import os, subprocess,re, sys,glob
from collections import OrderedDict,namedtuple
from array import array

from SelectGiBUUEvents import *

outdir = "/exp/minerva/app/users/granados/cmtuser/MATAna/cc-ch-pip-ana/GiBUU/"
#Two: CH and H2O
#NukeCCPion_GiBUU_tracker_carbon_T0
#NukeCCPion_GiBUU_tracker_carbon_T1
#NukeCCPion_GiBUU_oxygen_T0
#NukeCCPion_GiBUU_oxygen_T1
#
#A scaling
targets = OrderedDict()
targets["tracker_carbon"] = ["tracker",12,1] 
#targets["oxygen"]         = ["water"  ,16,2] 

#Grab hydrogen
f = ROOT.TFile("{0}/CCPion_GiBUU_hydrogen.root".format(outdir))
hyd_muon_p     = f.Get("pmu") 
hyd_muon_pt    = f.Get("ptmu")    
hyd_muon_pz    = f.Get("pzmu") 
hyd_muon_theta = f.Get("thmu")
hyd_q2         = f.Get("q2")
hyd_W          = f.Get("wexp")      
hyd_pion_ekin  = f.Get("tpi")
hyd_pion_theta = f.Get("thpi")
hyd_enu        = f.Get("enu")

hyd_muon_p    .SetDirectory(0)  
hyd_muon_pt   .SetDirectory(0)
hyd_muon_pz   .SetDirectory(0)
hyd_muon_theta.SetDirectory(0)
hyd_q2        .SetDirectory(0)
hyd_W         .SetDirectory(0)
hyd_pion_ekin .SetDirectory(0)
hyd_pion_theta.SetDirectory(0)
hyd_enu       .SetDirectory(0)

for target in targets: 
  for model in ["T0","T1"]:
    #print target, model
    f= ROOT.TFile("{0}/CCPion_GiBUU_{1}_{2}.root".format(outdir,model,target))

    pmu         = f.Get("pmu") 
    ptmu        = f.Get("ptmu")    
    pzmu        = f.Get("pzmu") 
    thetamu_deg = f.Get("thmu")
    q2          = f.Get("q2")
    wexp        = f.Get("wexp")      
    tpi         = f.Get("tpi")
    thetapi_deg = f.Get("thpi")
    enu         = f.Get("enu")

    pmu        .SetDirectory(0)  
    ptmu       .SetDirectory(0)
    pzmu       .SetDirectory(0)
    thetamu_deg.SetDirectory(0)
    q2         .SetDirectory(0)
    wexp       .SetDirectory(0)
    tpi        .SetDirectory(0)
    thetapi_deg.SetDirectory(0)
    enu        .SetDirectory(0)

    pmu        .Scale(targets[target][1])  
    ptmu       .Scale(targets[target][1])
    pzmu       .Scale(targets[target][1])
    thetamu_deg.Scale(targets[target][1])
    q2         .Scale(targets[target][1])
    wexp       .Scale(targets[target][1])
    tpi        .Scale(targets[target][1])
    thetapi_deg.Scale(targets[target][1])
    enu        .Scale(targets[target][1])

    pmu        .Add(hyd_muon_p    ,targets[target][2])  
    ptmu       .Add(hyd_muon_pt   ,targets[target][2])
    pzmu       .Add(hyd_muon_pz   ,targets[target][2])
    thetamu_deg.Add(hyd_muon_theta,targets[target][2])
    q2         .Add(hyd_q2        ,targets[target][2])
    wexp       .Add(hyd_W         ,targets[target][2])
    tpi        .Add(hyd_pion_ekin ,targets[target][2])
    thetapi_deg.Add(hyd_pion_theta,targets[target][2])
    enu        .Add(hyd_enu       ,targets[target][2])

    pmu        .Scale( 1./ (targets[target][1]+targets[target][2]) )  
    ptmu       .Scale( 1./ (targets[target][1]+targets[target][2]) )
    pzmu       .Scale( 1./ (targets[target][1]+targets[target][2]) )
    thetamu_deg.Scale( 1./ (targets[target][1]+targets[target][2]) )
    q2         .Scale( 1./ (targets[target][1]+targets[target][2]) )
    wexp       .Scale( 1./ (targets[target][1]+targets[target][2]) )
    tpi        .Scale( 1./ (targets[target][1]+targets[target][2]) )
    thetapi_deg.Scale( 1./ (targets[target][1]+targets[target][2]) )
    enu        .Scale( 1./ (targets[target][1]+targets[target][2]) )

    print (targets[target][0], targets[target][1]+targets[target][2])

    output_file= ROOT.TFile("{0}/CCPion_GiBUU_{1}_{2}.root".format(outdir,model,targets[target][0]),"RECREATE")
    output_file.cd()

    pmu        .Write() 
    ptmu       .Write()
    pzmu       .Write()
    thetamu_deg.Write()
    q2         .Write()
    wexp       .Write()
    tpi        .Write()
    thetapi_deg.Write()
    enu        .Write()

    output_file.Close();

