import ROOT
import array
import sys
import math

def isCCQELike(mytree):

    #need 1 muon and no mesons and photons > 10 MeV in final state

    nfsp = mytree.nfsp
    Efsp = mytree.E
    pdg  = mytree.pdg

    n_muon   = 0
    n_meson  = 0
    n_baryon = 0
    n_gamma  = 0

    badMesons = [211,321,323,111,130,310,311,313,411,421]
    badBaryons = [3112,3122,3212,3222,4112,4122,4212,4222]

    for p in range(0,nfsp):
        if(pdg[p]==13): n_muon+=1
        if(abs(pdg[p]) in badMesons): n_meson+=1
        if(abs(pdg[p]) in badBaryons): n_baryon+=1
        if(pdg[p]==22 and Efsp[p] > 0.010): n_gamma+=1

    if n_muon == 1 and n_meson==0 and n_gamma==0 and n_baryon==0 : return True
    else: return False

def is1piplus(mytree):

    #need 1 muon and no mesons and photons > 10 MeV in final state

    nfsp = mytree.nfsp
    Efsp = mytree.E
    pdg  = mytree.pdg

    n_muon   = 0
    n_piplus = 0
    n_meson  = 0
    n_baryon = 0
    n_gamma  = 0
    tpi = 0

    badMesons = [211,321,323,111,130,310,311,313,411,421,431]
    badBaryons = [3112,3122,3212,3222,4112,4122,4212,4222]

    for p in range(0,nfsp):
        tpi = Efsp[p]-0.139569
        if(pdg[p]==13): n_muon+=1
        if(pdg[p]==211 and tpi > 0 and tpi < 0.35): n_piplus+=1
        if(abs(pdg[p]) in badMesons): n_meson+=1
        if(abs(pdg[p]) in badBaryons): n_baryon+=1
        if(pdg[p]==22 and Efsp[p] > 0.010): n_gamma+=1
        if(abs(pdg[p])==11): n_gamma+=1

    if n_muon == 1 and n_piplus==1 and n_meson==1: return True
    else: return False

def getProtonMomentum(mytree):
    nfsp = mytree.nfsp
    Efsp = mytree.E
    pdg  = mytree.pdg
    px = mytree.px
    py = mytree.py
    pz = mytree.pz
    
    prot_mom = ROOT.TVector3()
    prot_e = 0

    for p in range(0,nfsp):
        if(pdg[p]==2212 and Efsp[p]>prot_e):#highest energy proton
            prot_mom.SetX(px[p])
            prot_mom.SetY(py[p])
            prot_mom.SetZ(pz[p])
    return prot_mom

def getPionMomentum(mytree):
    nfsp = mytree.nfsp
    Efsp = mytree.E
    pdg  = mytree.pdg
    px = mytree.px
    py = mytree.py
    pz = mytree.pz
    
    pion_mom = ROOT.TVector3()
    pion_e = 0

    for p in range(0,nfsp):
        if(pdg[p]==211 and Efsp[p]>pion_e):#highest energy pion
            pion_mom.SetX(px[p])
            pion_mom.SetY(py[p])
            pion_mom.SetZ(pz[p])
            pion_e = Efsp[p]
    return pion_mom

def GetTpi(mytree):
    nfsp = mytree.nfsp
    Efsp = mytree.E
    pdg  = mytree.pdg
    px = mytree.px
    py = mytree.py
    pz = mytree.pz
    
    pion_mom = ROOT.TVector3()
    pion_e = 0
    tpi = 0
    for p in range(0,nfsp):
        if(pdg[p]==211 and Efsp[p]>pion_e):#highest energy pion
            tpi=Efsp[p]-0.139569 #Returns tpi in GeV
	    pion_e = Efsp[p]
    return tpi

def getMuonMomentum(mytree):
    muon_mom = ROOT.TVector3()

    nfsp = mytree.nfsp
    Efsp = mytree.E
    pdg  = mytree.pdg
    px = mytree.px
    py = mytree.py
    pz = mytree.pz
    
    for p in range(0,nfsp):
        if(pdg[p]==13):
            muon_mom.SetX(px[p])
            muon_mom.SetY(py[p])
            muon_mom.SetZ(pz[p])
            break
    return muon_mom

def GetElepTrue(mytree) : # // return Emu in GeV
    muon_mom = ROOT.TVector3()
    muon_E = 0
    nfsp = mytree.nfsp
    Efsp = mytree.E
    pdg  = mytree.pdg
    px = mytree.px
    py = mytree.py
    pz = mytree.pz
    
    for p in range(0,nfsp):
        if(pdg[p]==13):
            muon_E = Efsp[p]
            break
    return muon_E
 
def GetEnuTrue(mytree) : 
    enu = mytree.Enu_true
    return enu

def GetEhadTrue(mytree) : #  // return Ehad in GeV
    return GetEnuTrue(mytree) - GetElepTrue(mytree);

def GetQ2True(mytree) :
    q2=mytree.Q2
    return q2

def CalcWexp(Q2, Ehad) : # return Wexp in GeV
    W = pow(0.9383, 2.0) - Q2 + 2.0 * (0.9383) * Ehad;
    return W;

def GetWexpTrue(mytree) :  # return Wexp in GeV
    return CalcWexp(GetQ2True(mytree), GetEhadTrue(mytree)); 

def isCCInclusiveFHC(mytree):
    nu_pdg = mytree.PDGnu
    lep_pdg = mytree.PDGLep
    isCC = mytree.flagCCINC
    if(nu_pdg == 14 and lep_pdg == 13 and isCC == 1):
        return True
    
    return False

def isTKI(mytree):
    #<70 degree proton
    #<20 degree muon
    #1.5 to 10 muon mom
    #0.45 to 1.2 proton mom

    muonmom = getMuonMomentum(mytree)
    protonmom = getProtonMomentum(mytree)

    goodProtonMom = protonmom.Mag()>0.45 and protonmom.Mag()<1.2
    goodMuonMom   = muonmom.Mag()>1.5 and muonmom.Mag()<20

    goodProtonAngle = protonmom.Theta()*180/3.1415
    goodMuonAngle = muonmom.Theta()*180/3.1415

    if goodMuonMom and goodProtonAngle and goodMuonAngle: return True
    
    return False

def passMuonCuts(mytree):
    Pmu = getMuonMomentum(mytree)

    goodpmu = Pmu.Mag()>1.5 and Pmu.Mag()<20
    goodthetamu = Pmu.Theta() < 0.34906585
    if goodpmu and goodthetamu: return True
    return False

def passWexpCut(mytree):
    wexp = GetWexpTrue(mytree)
    if wexp > 0 and wexp < 1.4: return True
    return False

def IsSignal(mytree): 
    isSignal = True
    isSignal = isSignal and isCCInclusiveFHC(mytree)
    isSignal = isSignal and passMuonCuts(mytree)
    isSignal = isSignal and passWexpCut(mytree)
    isSignal = isSignal and is1piplus(mytree)
    return isSignal
   
def getPn(e):

    promom = getProtonMomentum(e)
    mumom = getMuonMomentum(e)

    protonPtVect = getProtonMomentum(e)
    muonPtVect = getMuonMomentum(e)

    protonPtVect.SetZ(0)
    muonPtVect.SetZ(0)


    pn = -99;
    Mn = 939.5654133/1000;
    Mp = 938.2720813/1000;
    Mmu = 105.6583745/1000;
    MA = 6*Mn + 6*Mp - 92.162/1000;
    Bin=27.13/1000;
    MAstar = MA - Mn + Bin;
    
    Epprim = ROOT.TMath.Sqrt(Mp*Mp + promom.Mag2());
    kprimL = mumom.Z();
    pprimL = promom.Z();
    Eprim = ROOT.TMath.Sqrt(Mmu*Mmu+mumom.Mag2());
    
    factor = MA - Eprim - Epprim + kprimL + pprimL;
    pT = (protonPtVect+muonPtVect).Mag();
    pL = -(MAstar*MAstar + pT*pT-factor*factor)/2.0/factor;
    

    pn = ROOT.TMath.Sqrt(pL*pL + pT*pT);
    return pn

inputFile = ROOT.TFile(sys.argv[1])

mytree = inputFile.Get("FlatTree_VARS")

tpibins = [0.,20,45.,60.,75.,100.,125.,166.,200.,350.]
pmubins = [0.,1000,2000,3000,4000,5500, 7500, 10000, 13000, 20000, 30000]
ptmubins = [0., 1.e2, 2.e2, 3.e2, 4.e2, 5.e2, 6.e2, 8.e2, 10.e2, 12.5e2, 15.e2, 25.e2]
enubins = [0., 1000, 3000, 4000, 6500, 9500, 14000, 30000]
thetamubins = [0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,14.,16.,20.]
pzmubins = [0., 1.e3, 2.e3, 3.e3, 4.e3, 5.e3, 6.0e3, 8.e3, 10.e3, 15.e3, 20.e3]
q2bins = [0, 0.025e6, 0.05e6, 0.1e6, 0.2e6, 0.3e6, 0.4e6, 0.5e6, 0.7e6, 1.0e6, 1.3e6, 2.0e6, 3.0e6]
thetapibins=[0., 15., 30., 45, 60, 76., 108., 122., 136., 150., 165., 180.]
wexpbins=[0., 10.e2, 11.e2, 12.e2, 13.e2, 14.e2, 15.e2]

q3bins = [0,0.2,0.3,0.4,0.6,0.9,1.2]
eabins = [0,0.04,0.08,0.12,0.16,0.24,0.32,0.4,0.6,0.8,1.0]


mytpi = ROOT.TH1D("tpi","tpi",len(tpibins)-1,array.array("d",tpibins))
mypmu = ROOT.TH1D("pmu","pmu",len(pmubins)-1,array.array("d",pmubins))
myptmu = ROOT.TH1D("ptmu","ptmu",len(ptmubins)-1,array.array("d",ptmubins))
mypzmu = ROOT.TH1D("pzmu","pzmu",len(pzmubins)-1,array.array("d",pzmubins))
myenu = ROOT.TH1D("enu","enu",len(enubins)-1,array.array("d",enubins))
mythetamu = ROOT.TH1D("thetamu","thetamu",len(thetamubins)-1,array.array("d",thetamubins))
mythetapi = ROOT.TH1D("thetapi","thetapi",len(thetapibins)-1,array.array("d",thetapibins))
myq2 = ROOT.TH1D("q2","q2",len(q2bins)-1,array.array("d",q2bins))
mywexp = ROOT.TH1D("wexp","wexp",len(wexpbins)-1,array.array("d",wexpbins))


#myptpz = ROOT.TH2D("ptpz","ptpz",len(pzbins)-1, array.array("d",pzbins),len(ptbins)-1,array.array("d",ptbins))
#mypt = ROOT.TH1D("pt","pt",len(ptbins)-1,array.array("d",ptbins))
#mypz = ROOT.TH1D("pz","pz",len(pzbins)-1,array.array("d",pzbins))
#myptpz_rate = ROOT.TH2D("ptpz_rate","ptpz_rate",len(pzbins)-1, array.array("d",pzbins),len(ptbins)-1,array.array("d",ptbins))
#mytheta = ROOT.TH1D("theta","theta",360,0,180)

#myptpz_bymode = []
#for i in range(0,60):
#    tmphist = ROOT.TH2D("ptpz-%d"%i,"ptpz-%d"%i,len(pzbins)-1, array.array("d",pzbins),len(ptbins)-1,array.array("d",ptbins))
#    myptpz_bymode.append(tmphist)

#myptpz_qelike = ROOT.TH2D("ptpz_qelike","ptpz_qelike",len(pzbins)-1, array.array("d",pzbins),len(ptbins)-1,array.array("d",ptbins))
counter = 0
for e in mytree:
   
    if IsSignal(e):
       fScaleFactor = e.fScaleFactor
       tpi = GetTpi(e)
       pmu = getMuonMomentum(e).Mag()
       thetamu = getMuonMomentum(e).Theta()
       ptmu = math.sin(thetamu)*pmu
       pzmu = math.cos(thetamu)*pmu
       enu = GetEnuTrue(e)
       thetapi = getPionMomentum(e).Theta()
       q2 = GetQ2True(e)
       wexp = GetWexpTrue(e)

       mytpi.Fill(tpi*1000,fScaleFactor)
       mypmu.Fill(pmu*1000,fScaleFactor)
       mythetamu.Fill(thetamu*180.0/3.14159,fScaleFactor)
       myptmu.Fill(ptmu*1000,fScaleFactor)
       mypzmu.Fill(pzmu*1000,fScaleFactor)
       myenu.Fill(enu*1000,fScaleFactor)
       mythetapi.Fill(thetapi*180.0/3.14159,fScaleFactor)
       myq2.Fill(q2*1000000,fScaleFactor)
       mywexp.Fill(wexp*1000,fScaleFactor)

       counter+=1
#    if coslep<0.93969262078: continue
#    elep= e.ELep
#    fScaleFactor = e.fScaleFactor
  
#    P = ROOT.TMath.Sqrt(elep*elep-0.105*0.105)
#    Pl = coslep*P
#    Pt = ROOT.TMath.Sqrt(1-coslep*coslep)*P

#    print P,Pl,Pt,ROOT.TMath.Sqrt(Pl*Pl+Pt*Pt)/P,e.Mode

#    myptpz.Fill(Pl,Pt,fScaleFactor)
#    mypt.Fill(Pt,fScaleFactor)
#    mypz.Fill(Pl,fScaleFactor)
#    myptpz_rate.Fill(Pl,Pt)
#    myptpz_bymode[e.Mode].Fill(Pl,Pt,fScaleFactor)
#    if(Pl>1.5 and Pl<60 and Pt>0 and Pt<4.5): mytheta.Fill(ROOT.TMath.ACos(coslep)*180.0/3.14159,fScaleFactor)
#    if(counter==100):
#        break
#    if isCCQELike(e): myptpz_qelike.Fill(Pl,Pt,fScaleFactor)

#myptpz.Scale(1,"width")
#mypt.Scale(1,"width")
#mypz.Scale(1,"width")

#myptpz.GetXaxis().SetTitle("Muon p_{||} (GeV)")
#myptpz.GetYaxis().SetTitle("Muon p_{t} (GeV)")

#mypz.GetXaxis().SetTitle("Muon p_{||} (GeV)")
#mypz.GetYaxis().SetTitle("d#sigma/dp_{||} cm^2/GeV/nucleon")

#mypt.GetXaxis().SetTitle("Muon p_{t} (GeV)")
#mypt.GetYaxis().SetTitle("d#sigma/dp_{t} cm^2/GeV/nucleon")


myoutput = ROOT.TFile(sys.argv[2],"RECREATE")
mytpi.Write()
mypmu.Write()
mythetamu.Write()
myptmu.Write()
mypzmu.Write()
myenu.Write()
mythetapi.Write()
myq2.Write()
mywexp.Write()




#myptpz.Write();
#mypz.Write()
#mypt.Write()
#myptpz_rate.Write()
#mytheta.Write()
for i in range(60):
    myptpz_bymode[i].Write()

#can = ROOT.TCanvas()
#myptpz.Draw("COLZ")

#can2 = ROOT.TCanvas()
#mypt.Draw()

#can3 = ROOT.TCanvas()
#mypz.Draw()


#raw_input("done")

    
