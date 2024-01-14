# ===============================================================================
# Use this script to make predictions for your analysis from different
# generators and their events created using the MINERvA flux.
#
# This script is originally from Dan. And he already went to the trouble of
# creating these events (root files)) for several models. Use those files as
# input to this script.
#
# E.g. nohup python -u make_ccpi_prediction_nuisance.py \\
# /pnfs/minerva/persistent/Models/GENIE/Medium_Energy/FHC/v3_0_6/tracker/
#   G18_02b_02_11a/CH/flat_GENIE_G18_02b_02_11a_50Mcombined.root
# nuisance_test_out.root > nuisance_test.log&
#
# Dan's files are located here: /pnfs/minerva/persistent/Models/
# ===============================================================================
import ROOT
import array
import sys
import math

M_PIP = 139.57039 # MeV

def isclose(a, b, rel_tol=1e-9, abs_tol=0.0):
    """
    Determine if two numbers are close in value.

    Parameters:
    - a, b: Numbers to compare.
    - rel_tol: Relative tolerance.
    - abs_tol: Absolute tolerance.

    Returns:
    - True if the numbers are close, False otherwise.
    """
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def isCC1Pi(mytree):
    # need 1 muon and 1 meson

    nfsp = mytree.nfsp
    pdg = mytree.pdg

    n_muon = 0
    n_meson = 0
    n_pip = 0

    badMesons = [321, 323, 111, 130, 310, 311, 313, 411, 421]

    for p in range(0, nfsp):
        if pdg[p] == 13:
            n_muon += 1
        if abs(pdg[p]) in badMesons:
            n_meson += 1
        if pdg[p] == -211:
            n_meson += 1
        if pdg[p] == 211:
            n_pip += 1

    return n_muon == 1 and n_meson == 0 and n_pip == 1


def getTpi(mytree):
    # highest energy pion

    nfsp = mytree.nfsp
    pdg = mytree.pdg
    Efsp = mytree.E
    px = mytree.px
    py = mytree.py
    pz = mytree.pz

    pip_mom = ROOT.TVector3()
    tpi = 0

    for p in range(0, nfsp):
        E = Efsp[p]*1000.
        if pdg[p] == 211 and E - M_PIP > tpi:
            pip_mom.SetX(px[p]*1000.)
            pip_mom.SetY(py[p]*1000.)
            pip_mom.SetZ(pz[p]*1000.)
            tpi = E - M_PIP
            try:
                assert(isclose(math.sqrt(pip_mom.Mag()**2 + M_PIP**2), E, rel_tol=0.001))
            except AssertionError:
                print("inconsistent", ecalc, E)
          
    return tpi


def isCCQELike(mytree):
    # need 1 muon and no mesons and photons > 10 MeV in final state

    nfsp = mytree.nfsp
    Efsp = mytree.E
    pdg = mytree.pdg

    n_muon = 0
    n_meson = 0
    n_baryon = 0
    n_gamma = 0

    badMesons = [211, 321, 323, 111, 130, 310, 311, 313, 411, 421]
    badBaryons = [3112, 3122, 3212, 3222, 4112, 4122, 4212, 4222]

    for p in range(0, nfsp):
        if pdg[p] == 13:
            n_muon += 1
        if abs(pdg[p]) in badMesons:
            n_meson += 1
        if abs(pdg[p]) in badBaryons:
            n_baryon += 1
        if pdg[p] == 22 and Efsp[p] > 0.010:
            n_gamma += 1

    if n_muon == 1 and n_meson == 0 and n_gamma == 0 and n_baryon == 0:
        return True
    else:
        return False


def getProtonMomentum(mytree):
    nfsp = mytree.nfsp
    Efsp = mytree.E
    pdg = mytree.pdg
    px = mytree.px
    py = mytree.py
    pz = mytree.pz

    prot_mom = ROOT.TVector3()
    prot_e = 0

    for p in range(0, nfsp):
        if pdg[p] == 2212 and Efsp[p] > prot_e:  # highest energy proton
            prot_mom.SetX(px[p])
            prot_mom.SetY(py[p])
            prot_mom.SetZ(pz[p])
    return prot_mom


def getMuonMomentum(mytree):
    muon_mom = ROOT.TVector3()

    nfsp = mytree.nfsp
    Efsp = mytree.E
    pdg = mytree.pdg
    px = mytree.px
    py = mytree.py
    pz = mytree.pz

    for p in range(0, nfsp):
        if pdg[p] == 13:
            muon_mom.SetX(px[p]*1000.)
            muon_mom.SetY(py[p]*1000.)
            muon_mom.SetZ(pz[p]*1000.)
            break
    return muon_mom


def isInclusiveFHC(mytree):
    nu_pdg = mytree.PDGnu
    lep_pdg = mytree.PDGLep

    if nu_pdg == 14 and lep_pdg == 13:
        return True

    return False


def isTKI(mytree):
    # <70 degree proton
    # <20 degree muon
    # 1.5 to 10 muon mom
    # 0.45 to 1.2 proton mom

    muonmom = getMuonMomentum(mytree)
    protonmom = getProtonMomentum(mytree)

    goodProtonMom = protonmom.Mag() > 0.45 and protonmom.Mag() < 1.2
    goodMuonMom = muonmom.Mag() > 1.5 and muonmom.Mag() < 10

    goodProtonAngle = protonmom.Theta() * 180 / 3.1415
    goodMuonAngle = muonmom.Theta() * 180 / 3.1415

    if goodProtonMom and goodMuonMom and goodProtonAngle and goodMuonAngle:
        return True

    return False


def getPn(e):
    promom = getProtonMomentum(e)
    mumom = getMuonMomentum(e)

    protonPtVect = getProtonMomentum(e)
    muonPtVect = getMuonMomentum(e)

    protonPtVect.SetZ(0)
    muonPtVect.SetZ(0)

    pn = -99
    Mn = 939.5654133 / 1000
    Mp = 938.2720813 / 1000
    Mmu = 105.6583745 / 1000
    MA = 6 * Mn + 6 * Mp - 92.162 / 1000
    Bin = 27.13 / 1000
    MAstar = MA - Mn + Bin

    Epprim = ROOT.TMath.Sqrt(Mp * Mp + promom.Mag2())
    kprimL = mumom.Z()
    pprimL = promom.Z()
    Eprim = ROOT.TMath.Sqrt(Mmu * Mmu + mumom.Mag2())

    factor = MA - Eprim - Epprim + kprimL + pprimL
    pT = (protonPtVect + muonPtVect).Mag()
    pL = -(MAstar * MAstar + pT * pT - factor * factor) / 2.0 / factor

    pn = ROOT.TMath.Sqrt(pL * pL + pT * pT)
    return pn


inputFile = ROOT.TFile(sys.argv[1])

mytree = inputFile.Get("FlatTree_VARS")


tpibins = [0., 20, 45., 60., 75., 100., 125., 166., 200., 350.] # MeV
mytpi = ROOT.TH1D("tpi", "tpi", len(tpibins) - 1, array.array("d", tpibins))

'''
  ptbins = [
      0,
      0.075,
      0.15,
      0.25,
      0.325,
      0.4,
      0.475,
      0.55,
      0.7,
      0.85,
      1,
      1.25,
      1.5,
      2.5,
      4.5,
  ]
  pzbins = [
      1.5,
      2.0,
      2.5,
      3.0,
      3.5,
      4.0,
      4.5,
      5.0,
      6.0,
      7.0,
      8.0,
      9.0,
      10.0,
      15.0,
      20.0,
      40.0,
      60.0,
  ]

  q3bins = [0, 0.2, 0.3, 0.4, 0.6, 0.9, 1.2]
  eabins = [0, 0.04, 0.08, 0.12, 0.16, 0.24, 0.32, 0.4, 0.6, 0.8, 1.0]

  myptpz = ROOT.TH2D(
      "ptpz",
      "ptpz",
      len(pzbins) - 1,
      array.array("d", pzbins),
      len(ptbins) - 1,
      array.array("d", ptbins),
  )
  mypt = ROOT.TH1D("pt", "pt", len(ptbins) - 1, array.array("d", ptbins))
  mypz = ROOT.TH1D("pz", "pz", len(pzbins) - 1, array.array("d", pzbins))
  myptpz_rate = ROOT.TH2D(
      "ptpz_rate",
      "ptpz_rate",
      len(pzbins) - 1,
      array.array("d", pzbins),
      len(ptbins) - 1,
      array.array("d", ptbins),
  )
  mytheta = ROOT.TH1D("theta", "theta", 360, 0, 180)

  myptpz_bymode = []
  for i in range(0, 60):
      tmphist = ROOT.TH2D(
          "ptpz-%d" % i,
          "ptpz-%d" % i,
          len(pzbins) - 1,
          array.array("d", pzbins),
          len(ptbins) - 1,
          array.array("d", ptbins),
      )
      myptpz_bymode.append(tmphist)

  myptpz_qelike = ROOT.TH2D(
      "ptpz_qelike",
      "ptpz_qelike",
      len(pzbins) - 1,
      array.array("d", pzbins),
      len(ptbins) - 1,
      array.array("d", ptbins),
  )
'''
counter_tot = 0
counter_sig = 0
print("looping events")
for e in mytree:

    counter_tot += 1
    #if(counter_tot > 10000):
    #    break

    if not isInclusiveFHC(e):
        continue

    # muon angle 20 deg
    coslep = e.CosLep
    if coslep < 0.93969262078:
        continue

    # muon momentum
    if not (1500. < getMuonMomentum(e).Mag() < 10000.):
        continue
      
    # CC1pi final state
    if not isCC1Pi(e):
        continue

    # pion tpi cut
    tpi = getTpi(e)
    if tpi > 350.:
        continue

    '''
      elep = e.ELep

      P = ROOT.TMath.Sqrt(elep * elep - 0.105 * 0.105)
      Pl = coslep * P
      Pt = ROOT.TMath.Sqrt(1 - coslep * coslep) * P

      #    print P,Pl,Pt,ROOT.TMath.Sqrt(Pl*Pl+Pt*Pt)/P,e.Mode

      myptpz.Fill(Pl, Pt, fScaleFactor)
      mypt.Fill(Pt, fScaleFactor)
      mypz.Fill(Pl, fScaleFactor)
      myptpz_rate.Fill(Pl, Pt)
      myptpz_bymode[e.Mode].Fill(Pl, Pt, fScaleFactor)
      if Pl > 1.5 and Pl < 60 and Pt > 0 and Pt < 4.5:
          mytheta.Fill(ROOT.TMath.ACos(coslep) * 180.0 / 3.14159, fScaleFactor)
      if isCCQELike(e):
          myptpz_qelike.Fill(Pl, Pt, fScaleFactor)
    '''

    counter_sig += 1

    fScaleFactor = e.fScaleFactor
    mytpi.Fill(tpi,fScaleFactor)
    
print("done looping events")
print(counter_sig, "signal events.", counter_tot, "total events.")

# myptpz.Scale(1,"width")
# mypt.Scale(1,"width")
# mypz.Scale(1,"width")
# mytpi.Scale(1,"width")

'''
myptpz.GetXaxis().SetTitle("Muon p_{||} (GeV)")
myptpz.GetYaxis().SetTitle("Muon p_{t} (GeV)")

mypz.GetXaxis().SetTitle("Muon p_{||} (GeV)")
mypz.GetYaxis().SetTitle("d#sigma/dp_{||} cm^2/GeV/nucleon")

mypt.GetXaxis().SetTitle("Muon p_{t} (GeV)")
mypt.GetYaxis().SetTitle("d#sigma/dp_{t} cm^2/GeV/nucleon")
'''

mytpi.GetXaxis().SetTitle("T_{#pi} (MeV)")
mytpi.GetYaxis().SetTitle("d#sigma/dT_{#pi} cm^{2}/MeV/nucleon")


myoutput = ROOT.TFile(sys.argv[2], "RECREATE")
#myptpz.Write()
#mypz.Write()
#mypt.Write()
#myptpz_rate.Write()
#mytheta.Write()
mytpi.Write()
#for i in range(60):
#    myptpz_bymode[i].Write()

can = ROOT.TCanvas()
mytpi.Draw("COLZ")

# can2 = ROOT.TCanvas()
# mypt.Draw()

# can3 = ROOT.TCanvas()
# mypz.Draw()


# raw_input("done")
