import ROOT
from ROOT import PlotUtils

variables = ["h_muon_p","h_muon_pt","h_muon_pz","h_muon_theta","h_q2","h_W","h_pion_ekin","h_pion_theta"]

filename = "/minerva/data/users/abercell/hists/Models/NukeCCPion_GiBUU_hydrogen.root"
fin = ROOT.TFile(filename)

hvars=[]
for var in variables:
  hvars.append(fin.Get(var))
  hvars[-1].SetDirectory(0)
  hvars[-1].Scale(1./70000)

fin.Close()

fout = ROOT.TFile(filename,"RECREATE")
fout.cd()
for h in hvars:
  h.Write()
fout.Close()
