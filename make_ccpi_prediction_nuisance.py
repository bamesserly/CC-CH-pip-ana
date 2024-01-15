# ==============================================================================
# Use this script to make predictions for your analysis from different
# generators and their events created using the MINERvA flux.
#
# This script is originally from Dan. And he already went to the trouble of
# creating these events (root files) for several models. Use those files as
# input to this script.
#
# E.g. nohup python -u make_ccpi_prediction_nuisance.py \\
# /pnfs/minerva/persistent/Models/GENIE/Medium_Energy/FHC/v3_0_6/tracker/
#   G18_02b_02_11a/CH/flat_GENIE_G18_02b_02_11a_50Mcombined.root
# nuisance_test_out.root > nuisance_test.log&
#
# Dan's files are located here: /pnfs/minerva/persistent/Models/
#
# I believe the fScaleFactor is what allows us to fill the xsec histogram
# directly, taking care of flux and targets norm.
#
# My policies/conventions:
# - Everything in MeV as soon as it's accessed. Nothing in GeV ever.
# - BWN at the plotting stage. Never before. So do no BWN here.
# ==============================================================================
import ROOT
import array
import sys
import math

M_PIP = 139.57039  # MeV/c^2
M_MU = 105.6583745  # MeV/c^2
M_P = 938.2720813  # MeV/c^2


# ==============================================================================
# utility functions
# ==============================================================================
def isClose(a, b, rel_tol=1e-9, abs_tol=0.0):
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


def radToDeg(rad):
    return 180.0 * rad / math.pi


# ==============================================================================
# General analysis
# ==============================================================================
def isInclusiveFHC(mytree):
    nu_pdg = mytree.PDGnu
    lep_pdg = mytree.PDGLep

    if nu_pdg == 14 and lep_pdg == 13:
        return True

    return False


# return TVector3 in MeV/c
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
            muon_mom.SetX(px[p] * 1000.0)
            muon_mom.SetY(py[p] * 1000.0)
            muon_mom.SetZ(pz[p] * 1000.0)
            break

    try:
        elep = mytree.ELep * 1000.0
        pmu_calc = ROOT.TMath.Sqrt(elep**2 - M_MU**2)
        assert isClose(muon_mom.Mag(), pmu_calc, rel_tol=0.001)
    except AssertionError:
        print("inconsistent", muon_mom.Mag(), pmu_calc)

    return muon_mom


# ==============================================================================
# CCPip functions
# ==============================================================================
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


# MeV/c
def getPpi(mytree):
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
        E = Efsp[p] * 1000.0
        if pdg[p] == 211 and E - M_PIP > tpi:
            pip_mom.SetX(px[p] * 1000.0)
            pip_mom.SetY(py[p] * 1000.0)
            pip_mom.SetZ(pz[p] * 1000.0)
            tpi = E - M_PIP
            try:
                assert isClose(
                    math.sqrt(pip_mom.Mag() ** 2 + M_PIP**2), E, rel_tol=0.001
                )
            except AssertionError:
                print("inconsistent", ecalc, E)

    return pip_mom


# MeV
def getTpi(mytree):
    # highest energy pion

    nfsp = mytree.nfsp
    pdg = mytree.pdg
    Efsp = mytree.E
    px = mytree.px
    py = mytree.py
    pz = mytree.pz

    pip_mom = getPpi(mytree)

    return math.sqrt(pip_mom.Mag() ** 2 + M_PIP**2) - M_PIP


# inputs, as always, should be in MeV and rad, output in MeV^2/c^2
def calcQ2(Enu, Emu, Thetamu):
    Q2 = (
        2.0 * Enu * (Emu - math.sqrt(Emu**2 - M_MU**2) * math.cos(Thetamu))
        - M_MU**2
    )
    assert Q2 > 0.0
    return Q2


def getQ2CCPi(mytree):
    return calcQ2(
        mytree.Enu_true * 1000.0, mytree.ELep * 1000.0, math.acos(mytree.CosLep)
    )


def calcWexp(Q2, Ehad):
    W2 = M_P**2 - Q2 + 2.0 * M_P * Ehad
    assert W2 > 0.0
    return math.sqrt(W2)


def getWexpCCPi(mytree):
    return calcWexp(mytree.Q2 * 1e6, mytree.Enu_true * 1000.0 - mytree.ELep * 1000.0)


# ==============================================================================
# CCQElike
# ==============================================================================
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


# ==============================================================================
# main
# ==============================================================================
inputFile = ROOT.TFile(sys.argv[1])

mytree = inputFile.Get("FlatTree_VARS")

# MeV
tpibins = [0.0, 20, 45.0, 60.0, 75.0, 100.0, 125.0, 166.0, 200.0, 350.0]
mytpi = ROOT.TH1D("tpi", "tpi", len(tpibins) - 1, array.array("d", tpibins))

# deg
thpibins = [0.0, 15.0, 30.0, 45, 60, 76.0, 108.0, 122.0, 136.0, 150.0, 165.0, 180.0]
mythpi = ROOT.TH1D("thpi", "thpi", len(thpibins) - 1, array.array("d", thpibins))

# MeV^2/c^2
q2bins = [
    0,
    0.025e6,
    0.05e6,
    0.1e6,
    0.2e6,
    0.3e6,
    0.4e6,
    0.5e6,
    0.7e6,
    1.0e6,
    1.3e6,
    2.0e6,
    3.0e6,
]
myq2 = ROOT.TH1D("q2", "q2", len(q2bins) - 1, array.array("d", q2bins))

# MeV/c
pmubins = [
    0.0,
    1.0e3,
    2.0e3,
    3.0e3,
    4.0e3,
    5.5e3,
    7.5e3,
    10.0e3,
    13.0e3,
    20.0e3,
    30.0e3,
]
mypmu = ROOT.TH1D("pmu", "pmu", len(pmubins) - 1, array.array("d", pmubins))

# MeV
enubins = [0.0, 1.0e3, 3.0e3, 4.0e3, 6.5e3, 9.5e3, 14.0e3, 30.0e3]
myenu = ROOT.TH1D("enu", "enu", len(enubins) - 1, array.array("d", enubins))

# MeV/c^2
wexpbins = [0.0, 10.0e2, 11.0e2, 12.0e2, 13.0e2, 14.0e2, 15.0e2]
mywexp = ROOT.TH1D("wexp", "wexp", len(wexpbins) - 1, array.array("d", wexpbins))

# MeV/c
ptmubins = [
    0.0,
    1.0e2,
    2.0e2,
    3.0e2,
    4.0e2,
    5.0e2,
    6.0e2,
    8.0e2,
    10.0e2,
    12.5e2,
    15.0e2,
    25.0e2,
]
myptmu = ROOT.TH1D("ptmu", "ptmu", len(ptmubins) - 1, array.array("d", ptmubins))

# MeV/c
pzmubins = [
    0.0,
    1.0e3,
    2.0e3,
    3.0e3,
    4.0e3,
    5.0e3,
    6.0e3,
    8.0e3,
    10.0e3,
    15.0e3,
    20.0e3,
]
mypzmu = ROOT.TH1D("pzmu", "pzmu", len(pzmubins) - 1, array.array("d", pzmubins))

counter_tot = 0
counter_sig = 0
print("looping events")
for e in mytree:
    counter_tot += 1
    # if(counter_tot > 1000):
    #    break

    # CC numu
    if not isInclusiveFHC(e):
        continue

    # muon angle 20 deg
    coslep = e.CosLep
    if coslep < 0.93969262078:
        continue

    # muon momentum 1.5-10 GeV
    pmu = getMuonMomentum(e)
    if not (1500.0 < pmu.Mag() < 10000.0):
        continue

    # CC1pi final state
    if not isCC1Pi(e):
        continue

    # tpi < 350 MeV
    tpi = getTpi(e)
    if tpi > 350.0:
        continue

    q2 = e.Q2 * 1.0e6
    q2_calc = getQ2CCPi(e)
    assert isClose(q2, q2_calc, rel_tol=0.01)

    wexp = getWexpCCPi(e)

    if wexp > 1400.0:
        continue

    counter_sig += 1

    fScaleFactor = e.fScaleFactor

    pzmu = coslep * pmu.Mag()
    ptmu = ROOT.TMath.Sqrt(1 - coslep * coslep) * pmu.Mag()
    thpi = radToDeg(getPpi(e).Theta())

    mytpi.Fill(tpi, fScaleFactor)
    mythpi.Fill(thpi, fScaleFactor)
    mypmu.Fill(pmu.Mag(), fScaleFactor)
    myptmu.Fill(ptmu, fScaleFactor)
    mypzmu.Fill(pzmu, fScaleFactor)
    myq2.Fill(q2, fScaleFactor)
    myenu.Fill(mytree.Enu_true * 1000.0, fScaleFactor)
    mywexp.Fill(wexp, fScaleFactor)

print("done looping events")
print(counter_sig, "signal events.", counter_tot, "total events.")

# mytpi.Scale(1,"width") # no BWN at this stage!

mytpi.GetXaxis().SetTitle("T_{#pi} (MeV)")
mytpi.GetYaxis().SetTitle("d#sigma/dT_{#pi} cm^{2}/MeV/nucleon")

mythpi.GetXaxis().SetTitle("#theta_{#pi} (MeV)")
mythpi.GetYaxis().SetTitle("d#sigma/d#theta_{#pi} cm^{2}/deg/nucleon")

mypmu.GetXaxis().SetTitle("Muon p (MeV/c)")
mypmu.GetYaxis().SetTitle("d#sigma/dp cm^{2}/MeV/nucleon")

myptmu.GetXaxis().SetTitle("Muon p_{t} (MeV/c)")
myptmu.GetYaxis().SetTitle("d#sigma/dp_{t} cm^2/MeV/nucleon")

mypzmu.GetXaxis().SetTitle("Muon p_{||} (MeV/c)")
mypzmu.GetYaxis().SetTitle("d#sigma/dp_{||} cm^2/MeV/nucleon")

myq2.GetXaxis().SetTitle("Q^2 (MeV^2)")
myq2.GetYaxis().SetTitle("d#sigma/dQ^} cm^{2}/MeV^2/nucleon")

myenu.GetXaxis().SetTitle("E_{#nu} (MeV)")
myenu.GetYaxis().SetTitle("d#sigma/dE_{#nu} cm^{2}/MeV/nucleon")

mywexp.GetXaxis().SetTitle("W_{exp} (MeV/c^2)")
mywexp.GetYaxis().SetTitle("d#sigma/dW_{exp} cm^{2}/MeV/nucleon")

myoutput = ROOT.TFile(sys.argv[2], "RECREATE")

mytpi.Write()
mythpi.Write()
mypmu.Write()
myptmu.Write()
mypzmu.Write()
myq2.Write()
myenu.Write()
mywexp.Write()


"""

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
"""

"""
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
"""

"""
  for i in range(60):
      myptpz_bymode[i].Write()

  can = ROOT.TCanvas()
  mytpi.Draw("COLZ")

   can2 = ROOT.TCanvas()
   mypt.Draw()

   can3 = ROOT.TCanvas()
   mypz.Draw()


   raw_input("done")
"""
