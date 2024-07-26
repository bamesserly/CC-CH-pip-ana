import os, subprocess,re, sys,glob,csv
from collections import OrderedDict,namedtuple
import ROOT

datfields = ["Run","Event","ID","Charge","perweight","x","y","z","E","px","py","pz","history","production_ID","enu"]

#All is GeV
#https://gibuu.hepforge.org/trac/wiki/ParticleIDs
#Hadron
partmass = dict()
partmass[1 ] = 0.938
partmass[2 ] = 1.232
partmass[3 ] = 1.462
partmass[4 ] = 1.534
partmass[5 ] = 1.659
partmass[6 ] = 1.928
partmass[7 ] = 1.524
partmass[8 ] = 1.737
partmass[9 ] = 1.804
partmass[10] = 1.676
partmass[11] = 2.127
partmass[12] = 1.717
partmass[13] = 1.885
partmass[14] = 1.717
partmass[15] = 1.879
partmass[16] = 1.684
partmass[17] = 1.903
partmass[18] = 2.086
partmass[19] = 1.672
partmass[20] = 1.920
partmass[21] = 1.762
partmass[22] = 2.057
partmass[23] = 1.956
partmass[24] = 2.171
partmass[25] = 1.744
partmass[26] = 1.882
partmass[27] = 1.706
partmass[28] = 2.014
partmass[29] = 1.752
partmass[30] = 1.881
partmass[31] = 1.945
partmass[32] = 1.116
partmass[33] = 1.189
partmass[34] = 1.385
partmass[35] = 1.405
partmass[36] = 1.520
partmass[37] = 1.600
partmass[38] = 1.670
partmass[39] = 1.690
partmass[40] = 1.810
partmass[41] = 1.820
partmass[42] = 1.830
partmass[43] = 1.670
partmass[44] = 1.775
partmass[45] = 2.030
partmass[46] = 1.800
partmass[47] = 1.890
partmass[48] = 2.100
partmass[49] = 2.110
partmass[50] = 1.660
partmass[51] = 1.750
partmass[52] = 1.915
partmass[53] = 1.315
partmass[54] = 1.530
partmass[55] = 1.672
partmass[56] = 2.285
partmass[57] = 2.452
partmass[58] = 2.520
partmass[59] = 2.466
partmass[60] = 2.645
partmass[61] = 2.697

#Meson
partmass[101] = 0.1380
partmass[102] = 0.5478
partmass[103] = 0.7755
partmass[104] = 0.8000
partmass[105] = 0.7826
partmass[106] = 0.9578
partmass[107] = 1.0194
partmass[108] = 2.9800
partmass[109] = 3.0969
partmass[110] = 0.4960
partmass[111] = 0.4960
partmass[112] = 0.8920
partmass[113] = 0.8920
partmass[114] = 1.8670
partmass[115] = 1.8670
partmass[116] = 2.0070
partmass[117] = 2.0070
partmass[118] = 1.9690
partmass[119] = 1.9690
partmass[120] = 2.1120
partmass[121] = 2.1120
partmass[122] = 1.2754

#Lepton
partmass[901] = 511e-6
partmass[902] = 0.1056 
partmass[903] = 1.7769 
partmass[911] = 0
partmass[912] = 0
partmass[913] = 0
partmass[997] = 80.4
partmass[998] = 91.2
partmass[999] = 0   

#Want ID, Charge, weight, momenta, Enu
eventtuple = namedtuple("eventtuple", "event enu wgt Q W pid charge threemom energy")
pievttuple = namedtuple("pievttuple", "pmu pmuZ pmuT thmu Q W tpi thpi enu")
pievtvars=["pmu","pmuZ","pmuT","thmu","Q","W","tpi","thpi","enu"]

def calcQ2( Enu, lepE, lepPz, mass ):
  return 2*Enu*(lepE-lepPz) - mass*mass

def calcWexp(Enu,lepE, Q2):
  Ehad = Enu - lepE 
  W2 = partmass[1]*partmass[1] - Q2 + 2*partmass[1]*Ehad
  return W2**.5 if W2 > 0 else -1

def isSignal(event):
  #What is the truth signal? True Variables only
  # 1 muon, 1 piplus, no other mesons, no gammas, any other of baryons

  # W < 1.4 GeV
  # Tpi between 35 and 350 MeV

  # I think natural units are GeV
  # pmu between 1.5 and 20 GeV
  # theta mu < 0.22689280275 radians (13 deg)

  # Calculate W, W < 1.4 GeV
  if event.W < 0 or 1.4 < event.W:
    return False,None

  nMuons = 0
  nPiPlusRange = 0
  nPiPlus = 0
  nMesons = 0
  nBaryons = 0
  nOther = 0

  pievtkine = {}
  for iPart in range(len(event.pid)):
    #Mu- 
    if event.pid[iPart] == 902 and event.charge[iPart] == -1:
      p3mu = event.threemom[iPart]           

      if p3mu.Mag() < 1.5 or 20 < p3mu.Mag():
        return False,None
      if p3mu.Theta() > 0.3491:
        return False,None

      nMuons += 1

      pievtkine["pmu"]  = p3mu.Mag() 
      pievtkine["pmuT"] = p3mu.Perp()
      pievtkine["pmuZ"] = p3mu.Pz()  
      pievtkine["thmu"] = p3mu.Theta()*180/3.14159265359
      continue
 
    #other leptons
    if event.pid[iPart] > 900:
      #Ignore gammas with energies < 10 MeV
      if event.pid[iPart] == 999 and event.energy[iPart] <= 0.01:
        continue
      nOther += 1
      continue

    if event.pid[iPart] < 100:
      nBaryons += 1
      continue

    if event.pid[iPart] > 100 and event.pid[iPart]<200:
      nMesons  += 1

    if event.pid[iPart] == 101 and event.charge[iPart] == 1:
      nPiPlus += 1
      tpi = event.energy[iPart]-partmass[101]

      if 0. < tpi and tpi < 0.35:
        nPiPlusRange += 1
        
        p3pi = event.threemom[iPart]

        pievtkine["tpi"]  = tpi
        pievtkine["thpi"] = p3pi.Theta()*180/3.14159265359
  
  if nMuons != 1:
    return False, None
  if nPiPlusRange != 1:
    return False, None
  if nPiPlus != 1:
    return False, None
  if nMesons != 1:
    return False, None
  if nOther != 0:
    return False, None

  pievtkine["Q"] = event.Q
  pievtkine["W"] = event.W
  pievtkine["enu"] = event.enu
 
  pivars = pievttuple(*[ pievtkine[_] for _ in pievtvars])
  return True, pivars
  

def readEvents(filename):
  iLine = 0
  with open(filename) as csvfile:
    reader = csv.DictReader(csvfile,fieldnames =datfields, delimiter=' ', skipinitialspace=True)
    oldEvent = 0
    #Setup
    pid      = []
    charge   = []
    threemomenta  = []
    energy   = []
   
    bFindLep = False
    lepE, lepPz = 0,0
 
    for row in reader:
      if iLine == 0:
        iLine += 1
        continue
      if iLine == 1:
        iLine += 1
        oldEvent = int(row["Event"])
        #Set things up
        event = int(row["Event"])
        wgt   = float(row["perweight"])
        enu   = float(row["enu"])

      #Debug
      #if iLine > 96:
      #  break
      #iLine += 1
      
      if int(row["Event"]) != oldEvent:
        #Calculated everything new
        Q2 = -1 
        W  = -1 

        if bFindLep:
          Q2 = calcQ2( enu, lepE, lepPz, partmass[902] )
          W  = calcWexp( enu, lepE, Q2 )

        #Pass the new event
        nextEvent=eventtuple(event,enu,wgt,Q2,W,pid,charge,threemomenta,energy)

        yield nextEvent 

        oldEvent = int(row["Event"])
        event = oldEvent
        bFindLep = False

        wgt   = float(row["perweight"])
        enu   = float(row["enu"])
        pid      = []
        charge   = []
        threemomenta  = []
        energy   = []

      tmpP3 = ROOT.TVector3( float(row["px"]), float(row["py"]), float(row["pz"]) )
      #tmpP3.RotateX( -0.05887 );#numi beam rad

      #Check the weight to find the initial nucleon (should be 0)
      if float(row["perweight"]) < 1e-10:
        continue
      else: 
        pid    .append( int(row["ID"]) )
        charge .append( int(row["Charge"]) ) 
        energy .append( float(row["E"]) ) 
        threemomenta.append( tmpP3 ) 

        #mu-?  If so, store the info
        if pid[-1] == 902 and charge[-1] == -1:
          lepE, lepPz = float(row["E"]), float(row["pz"])
          bFindLep = True
