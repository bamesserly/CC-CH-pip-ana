# Original Author is Heidi Schellman.
# Taken from https://minerva-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=25029
# I believe the log chi2 calc might be a way to handle Peele's Pertinent Puzzle.
import sys, os, string
import math

# version 10 implement new low q2 error band from /minerva/data/users/drut1186/Mateus_Pub_Inputs/Modified_AddedLowQ2Suppression
from ROOT import *

from PlotUtils import *

NOSYS = True

BIN = False

debin = True

Ratio = True

SCALE = 10 ^ 40
SCALE = 1.0

Mateus = False

minb = 0

q2 = "q2"

nounfold = False

basemodel = "notCV-"

varname = {
    "q2": "$Q^2$",
    "ptmu": "$p_\perp$",
    "pzmu": "p_\parallel",
    "pzmu_ptmu": "$p_\parallel - p_\perp$",
    "enu": "E_\nu",
}
basemodel = ""
if nounfold:
    basemodel = "nounfold-" + basemodel


def Debinwidth(h_data):
    h_nobin = h_data.Clone()
    h_nobin.SetDirectory(0)
    n = h_data.GetNbinsX()
    for i in range(0, n + 1):
        wid = h_data.GetBinWidth(i)
        h_nobin.SetBinContent(i, h_data.GetBinContent(i) * wid)
        h_nobin.SetBinError(i, h_data.GetBinError(i) * wid)
    return h_nobin


def Debinwidth2(h_data):
    h_nobin = h_data.Clone()
    h_nobin.SetDirectory(0)
    n = h_data.GetNbinsX()
    m = h_data.GetNbinsY()
    for i in range(0, n + 1):
        for j in range(0, m + 1):
            wid1 = h_data.GetXaxis().GetBinWidth(i)
            wid2 = h_data.GetYaxis().GetBinWidth(j)
            wid = wid1 * wid2
            if wid < 0:
                wid = 0
            h_nobin.SetBinContent(i, j, h_data.GetBinContent(i, j) * wid)
            h_nobin.SetBinError(i, j, h_data.GetBinError(i, j) * wid)
    return h_nobin


def readData(h_data, diag=False):
    # returns a transposed list

    size = h_data.GetNbinsX() * h_data.GetNbinsY()
    compact = TVectorD(size)

    npt = h_data.GetNbinsY() + 2
    npz = h_data.GetNbinsX() + 2

    ibin = 0
    for iy in range(1, npt - 1):
        for ix in range(1, npz - 1):
            compact[ibin] = h_data.GetBinContent(ix, iy)
            ibin += 1

    return compact


def readData1D(h_data, diag=False):
    # returns a transposed list

    size = h_data.GetNbinsX()
    compact = TVectorD(size)

    npz = h_data.GetNbinsX() + 2

    ibin = 0

    for ix in range(1, npz - 1):
        compact[ibin] = h_data.GetBinContent(ix)
        ibin += 1

    return compact


def readMatrix(covmx, h_data, diag=False, log=False):
    npz = h_data.GetNbinsX() + 2
    npt = h_data.GetNbinsY() + 2

    size = (npt - 2) * (npz - 2)

    compact = TMatrixD(size, size)
    ixx = -1
    for ix in range(npz, (npz) * (npt - 1)):
        if ix % (npz) == 0 or ix % (npz) == npz - 1:
            continue

        ixx += 1
        iyy = -1
        for iy in range(npz, (npz) * (npt - 1)):
            if iy % (npz) == 0 or iy % (npz) == npz - 1:
                continue

            iyy += 1

            if not log:
                compact[ixx][iyy] = covmx[ix][iy]
            else:
                compact[ixx][iyy] = covmx[ix][iy]

            if diag and ixx != iyy:
                compact[ixx][iyy] = 0
            # special to deal with completely empty bins

    return compact


def readMatrix1D(covmx, h_data, diag=False, log=False):
    print("READMATRIX")
    # covmx.

    npz = h_data.GetNbinsX() + 2

    size = npz - 2

    compact = TMatrixD(size, size)
    ixx = -1
    for ix in range(1, npz - 1):
        ixx = ix - 1

        for iy in range(1, npz - 1):
            iyy = iy - 1

            if not log:
                compact[ixx][iyy] = covmx[ix][iy]
            else:
                compact[ixx][iyy] = covmx[ix][iy]

            if diag and ixx != iyy:
                compact[ixx][iyy] = 0

    # compact.Print()
    return compact


def compressVector(v, e):
    vnew = TVectorD(len(e))
    for i in e:
        vnew[i] = v[e[i]]

    return vnew


def compressMatrix(v, e):
    vnew = TMatrixD(len(e), len(e))
    for i in e:
        for j in e:
            vnew[i][j] = v[e[i]][e[j]]

    return vnew


if basemodel == "":
    dir = "/Users/schellma/Dropbox/ccqe/PAPER6/CCQENu_v21r1p1_Pub_Sep2019ntuples_YESNue_Nov11_CV"
else:
    dir = "/Users/schellma/Dropbox/ccqe/PAPER6/CCQENu_v21r1p1_Pub_Sep2019ntuples_YESNue_Nov11_CV_RPA_Res_MINOS"

print(dir)
mcdir = "/Users/schellma/Dropbox/ccqe/PAPER6/CCQENu_v21r1p1_Pub_Sep2019ntuples_YESNue_Nov11_MODEL"

if Mateus:
    # dir = dir.replace("/Users/schellma/Dropbox/ccqe/MT/","/pnfs/minerva/persistent/users/mateusc/")
    # mcdir = mcdir.replace("/Users/schellma/Dropbox/ccqe/MT/","/pnfs/minerva/persistent/users/mateusc/")
    dir = dir.replace("MT", "FNAL")

# * CV/ (MinervaGENIE v1: +2p2h+recoil+RPA $\pi$tune)
# * default/  (+2p2h)
# * CV_RPA_Res_MINOS/ (+2p2h+recoil+RPA $\pi$tune+MINOS low q2 sup.)
# * CV_RPA_Res_Nieves/ (+2p2h+recoil+RPA $\pi$tune+Nieves low q2 sup.)
# * $\pi$tune/ (+2p2h $\pi$tune)
# * pion_2p2h/  (+2p2h+recoil fit $\pi$tune)
# * pion_rpa/  (+2p2h+RPA+ NonRES $\pi$tune)

#    For the record following are the root files and relevant histograms for the q2qe
#
# Nuwro:
# /pnfs/minerva/persistent/users/mateusc/CCQENu_v21r1p1_Pub_Sep2019ntuples_YESNue_Nov11_CV/QELike_ME_Mateus_Nuwro_v1902_600M.root
# > nuwroSF_kQsqQE_all
# > nuwroLFG_kQsqQE_all
#
# GIBUU:
# /pnfs/minerva/persistent/users/mateusc/CCQENu_v21r1p1_Pub_Sep2019ntuples_YESNue_Nov11_CV/MMECCQEGiBUU.root
# > q2qeall (inside folder MMECCQEGiBUUnu)
#
# Geniev3:
# /pnfs/minerva/persistent/users/mateusc/CCQENu_v21r1p1_Pub_Sep2019ntuples_YESNue_Nov11_CV/MMECCQEGENIEOOB_30000000Evts_20191003.root
# > q2qeall (inside folder MMECCQEGENIEOOBnu)
#
# We don't currently have text versions of the external results and they are not implemented in the scripts as they still look for the analysis name standards. i can try to add them to the scripts latter today.
var = "pzmu_ptmu"
if len(sys.argv) >= 2:
    var = sys.argv[1]


models = [
    "default",
    "CV",
    "piontune",
    "pion_rpa",
    "CV_RPA_Res_MINOS",
    "CV_RPA_Res_Nieves",
    "pion_2p2h",
    "CV_minerva_joint_lowq2",
]
# models = ["default","CV","CV_RPA_Res_MINOS","piontune","pion_2p2h","pion_rpa","CV_minerva_joint_lowq2"]

# models = ["CV","default"]
n = len(models)
bmodels = models
print(models)
models = []
for model in bmodels[0:n]:
    print("this is a model")
    models.append(model + "-no-2p2h")
for model in bmodels:
    models.append(model)

extramodels = ["GiBUU", "NuWroSF", "NuWroLFG", "Geniev3"]
if var != "enu":
    models = models + extramodels
# if "mu" in var:
# models.remove("GiBUU")
print(models)

translate = {
    "CV": "\qquad+RPA+$\pi$tune+recoil fit (MnvGENIEv1)",
    "CV_RPA_Res_MINOS": "\qquad+recoil fit+RPA+$\pi$tune+MINOS low $Q^2$ sup.",
    "default": "GENIE 2.12.6 + 2p2h",
    "piontune": "\qquad+$\pi$tune",
    "pion_rpa": "\qquad+RPA+$\pi$tune",
    "pion_2p2h": "\qquad+recoil fit+$\pi$tune",
    "CV_RPA_Res_Nieves": "\qquad+recoil fit+RPA+$\pi$tune+Nieves low $Q^2$ sup.",
    "CV-no-2p2h": "%+recoil fit+RPA+$\pi$tune",
    "CV_RPA_Res_MINOS-no-2p2h": "\qquad+RPA+$\pi$tune+MINOS low $Q^2$ sup.",
    "default-no-2p2h": "GENIE 2.12.6",
    "piontune-no-2p2h": "\qquad+$\pi$tune",
    "pion_rpa-no-2p2h": "\qquad+RPA+$\pi$tune",
    "pion_2p2h-no-2p2h": "%+recoil fit+$\pi$tune",
    "CV_RPA_Res_Nieves-no-2p2h": "%+recoil fit+RPA+$\pi$tune+Nieves low $Q^2$ sup.",
    "GiBUU": "GiBUU",
    "NuWroSF": "NuWro SF",
    "NuWroLFG": "NuWro LFG",
    "Geniev3": "GENIE v3",
    "CV_minerva_joint_lowq2": "\qquad+recoil fit+RPA $\pi$tune+MINERvA $\pi$tune+low $Q^2$ sup.",
    "CV_minerva_joint_lowq2-no-2p2h": "%+recoil fit + RPA + NonRES piontune + MINERvA $\pi$tune low $Q^2$ sup.",
}
# models = ["default"]
# models = ["$\pi$tune","pion_rpa","pion_2p2h"]

# if "mu" in var:
#  basehist = "Modified_Pzmu.root"
# else:
#  if "q2" in var:
#    basehist = "Modified_q2.root"
#  else:
#    print("no mod hist for this variable")
#    sys.exit(1)

basehist = "h_pzmu_ptmu_data_nobck_unfold_effcor_cross_section"

mcbasehist = "h_pzmu_ptmu_qelikecomp_cross_section"


if var not in ["q2", "enu", "pzmu", "ptmu", "pzmu_ptmu"]:
    print(" don't recognize input", var)
    sys.exit(0)

if var in ["q2", "enu"]:
    basehist = basehist.replace("pzmu", var)
    mcbasehist = mcbasehist.replace("pzmu", var)


proj = ""
if var in ["q2", "enu", "pzmu"]:
    proj = "_px"
if var in ["ptmu"]:
    proj = "_py"


thefile = "CrossSection_CombinedPlaylists-CombinedPlaylists-me1Bme1Cme1Dme1Eme1Fme1Gme1Lme1Mme1Nme1Ome1P_v3.root"


dataname = os.path.join(dir, thefile.replace("app", "v3"))

dataname = "/Users/schellma/Dropbox/ccqe/PAPER6/lowq2syst/Modified_Pzmu.root"

dataname = "/Users/schellma/Dropbox/ccqe/PAPER7/Modified_AddedLowQ2Suppression/Modified_DataXSecResults_pzptmu_December8th.root"

if not os.path.exists(dataname):
    print("data file does not exist", dataname)
    sys.exit(0)
if var == "q2":
    dataname = dataname.replace("pzptmu", "q2ptmu")
if var == "enu":
    dataname = dataname.replace("pzptmu", "q2ptmu")


rname = var + "_ratios.root"
if BIN:
    bin = "_bin"
else:
    bin = ""
outname = "chi2_" + var + "lowq2" + ".txt"
datafile = TFile.Open(os.path.join(dataname))

# print (basehist)
# datafile.ls()
datahist = MnvH2D()
datahist = datafile.Get(basehist).Clone()
datahist.SetDirectory(0)
datahist.Scale(SCALE)
if BIN:
    datahist.Scale(1.0, "width")
print("var is ", var)
# datadraw=datahist.GetCVHistoWithError()
# datadraw.Scale(1.,"width")
rfile = TFile.Open(rname, "RECREATE")
if proj == "_px":
    oneDhist = datahist.ProjectionX()
    datavals = readData1D(oneDhist)
    matrix = readMatrix1D(oneDhist.GetTotalErrorMatrix(), oneDhist)
    uncov = readMatrix(oneDhist.GetSysErrorMatrix("unfoldingCov"), oneDhist)
    # oneDhist.GetSysErrorMatrix("unfoldingCov").Print()
    # uncov.Print()
    if nounfold:
        matrix -= uncov
    matrix.Print()
    # oneDhist.SetName("data"+oneDhist.GetName())
    oneDhist.Write()
    datadraw = oneDhist.GetCVHistoWithError()
    datadraw.Print("ALL")
    oneDhist.Delete()
if proj == "_py":
    oneDhist = datahist.ProjectionY()
    datavals = readData1D(oneDhist)
    matrix = readMatrix1D(oneDhist.GetTotalErrorMatrix(), oneDhist)
    uncov = readMatrix(oneDhist.GetSysErrorMatrix("unfoldingCov"), oneDhist)
    # uncov.Print()
    if nounfold:
        matrix -= uncov
    matrix.Print()
    # oneDhist.SetName("data"+oneDhist.GetName())
    datadraw = oneDhist.GetCVHistoWithError()
    oneDhist.Write()
    oneDhist.Delete()
if proj == "":
    datavals = readData(datahist)
    matrix = readMatrix(datahist.GetTotalErrorMatrix(), datahist)
    datahist.Write()

# datavals.Print()


# datavals.Print()

goodelements = {}
ok = 0
for i in range(0, len(datavals)):
    if datavals[i] != 0:
        goodelements[ok] = i
        ok += 1
# print (goodelements)

data = compressVector(datavals, goodelements)
logdata = TVectorD(data)
for i in range(0, len(goodelements)):
    if data[i] == 0:
        continue
    logdata[i] = math.log(data[i])


mchist = {}
mcdraw = {}
mc = {}
logmc = {}
ofile = {}
#
widthcorr = True
for model in models:
    mcname = os.path.join(mcdir.replace("MODEL", model), thefile)

    if "no-2p2h" in model:
        mcname = mcname.replace("-no-2p2h", "")
    # print("opening", mcname)
    twoD = True
    if model == "GiBUU":
        mcname = (
            "/Users/schellma/Dropbox/ccqe/PAPER6/MINERvAMEGiBUUMMECCQE_2D_20191114.root"
        )
        mcbasehist = "MINERvAMEGiBUUMMECCQEnu/q2qeall"
        if var != q2:
            mcbasehist = "MINERvAMEGiBUUMMECCQEnu/muonptVSmupz_all"
            # /Users/schellma/Dropbox/ccqe/PAPER6/MINERvAMEGiBUUMMECCQE_2D_20191114.root

    if model == "NuWroSF":
        mcname = (
            "/Users/schellma/Dropbox/ccqe/PAPER5/QELike_ME_Mateus_Nuwro_v1902_600M.root"
        )
        mcbasehist = "nuwroSF_kQsqQE_all"
        if var != q2:
            mcbasehist = "nuwroSF_kMuonPtPz_all"

    if model == "NuWroLFG":
        mcname = (
            "/Users/schellma/Dropbox/ccqe/PAPER5/QELike_ME_Mateus_Nuwro_v1902_600M.root"
        )
        mcbasehist = "nuwroLFG_kQsqQE_all"
        if var != q2:
            mcbasehist = "nuwroLFG_kMuonPtPz_all"

    if model == "Geniev3":
        mcname = "/Users/schellma/Dropbox/ccqe/PAPER5/MMECCQEGENIEOOB_30000000Evts_20191003.root"
        mcbasehist = "MMECCQEGENIEOOBnu/q2qeall"
        if var != q2:
            mcbasehist = "MMECCQEGENIEOOBnu/muonptVSmupzall"

    if var == q2:
        twoD = False
    print(model, mcname, mcbasehist)
    try:
        mcfile = TFile.Open(mcname)
    except:
        print("mcfile does not exist ", mcfile)
        continue
    # print("mcfile is ",type(mcfile))
    if mcfile == None:
        print("mcfile does not exist ", mcfile)
        continue
    mchist[model] = MnvH2D()
    mchist[model] = mcfile.Get(mcbasehist).Clone()
    # mchist[model].SetName(model+"_"+mchist[model].GetTitle())
    mchist[model].Scale(SCALE)
    if BIN and widthcorr:
        mchist[model].Scale(1.0, "width")
    if debin and model in extramodels and not twoD:
        mchist[model] = Debinwidth(mchist[model])
    if debin and model in extramodels and twoD:
        print("should I shrink it ", proj)
        mchist[model] = Debinwidth2(mchist[model])
    mchist[model].Print()
    mchist[model].SetDirectory(0)
    # remove 2p2h from the default model
    if "no-2p2h" in model:
        print("---- remove the 2p2h component ----", model, translate[model])
        h2p2hcorr = MnvH2D()
        h2p2hcorrname = "h_pzmu_ptmu_2p2hcomp_cross_section"
        if var in ["q2", "enu"]:
            h2p2hcorrname = h2p2hcorrname.replace("pzmu", var)
        print("TH1", h2p2hcorrname)
        h2p2hcorr = mcfile.Get(h2p2hcorrname).Clone()
        h2p2hcorr.SetDirectory(0)
        h2p2hcorr.Print()
        mchist[model].Add(h2p2hcorr, -1.0)
        h2p2hcorr.Delete()

        mchist[model].Print()
        print(model, "changed\n\n")

    # mcdraw[model]=mchist[model].GetCVHistoWithError()
    # mcdraw[model].Scale(1.,"width")

    rfile.cd()
    if proj == "_px":
        if not model in extramodels:
            oneDhist = mchist[model].ProjectionX()
            mcdraw[model] = oneDhist.GetCVHistoWithError()

        else:
            oneDhist = mchist[model].Clone()

            mcdraw[model] = oneDhist.Clone()
        oneDhist.SetDirectory(0)
        mcdraw[model].SetDirectory(0)
        # mcdraw[model].Print("ALL")
        mcvals = readData1D(oneDhist)
        oneDhist.SetName(model + oneDhist.GetName())

        for bin in range(0, datadraw.GetNbinsX() + 1):
            mcdraw[model].SetBinError(bin, 0.0)
        oneDhist.Write()
        oneDhist.Delete()
    if proj == "_py":
        if not model in extramodels:
            oneDhist = mchist[model].ProjectionY()
            mcdraw[model] = oneDhist.GetCVHistoWithError()
        else:
            if var == q2:
                oneDhist = mchist[model].Clone()
                mcdraw[model] = oneDhist.Clone()
            else:  # pt
                oneDhist = mchist[model].ProjectionY()
                mcdraw[model] = oneDhist.Clone()
        mcdraw[model].SetDirectory(0)
        oneDhist.SetBinContent(0, 0)
        oneDhist.SetDirectory(0)
        oneDhist.SetName(model + oneDhist.GetName())
        mcvals = readData1D(oneDhist)
        # print (mcvals)
        oneDhist.Print("ALL")

        for bin in range(0, datadraw.GetNbinsX() + 1):
            mcdraw[model].SetBinError(bin, 0.0)
        oneDhist.Write()
        oneDhist.Delete()
    if proj == "":
        if not model in extramodels:
            mcvals = readData(mchist[model].GetCVHistoWithStatError())
        else:
            mcvals = readData(mchist[model])
        mchist[model].Write()
    mc[model] = compressVector(mcvals, goodelements)
    logmc[model] = TVectorD(mc[model])
    for i in range(0, len(goodelements)):
        if mc[model][i] == 0:
            continue
        print(mc[model][i])
        (logmc[model])[i] = math.log((mc[model])[i])
    # mchist[model].Delete()
    # mc[model].Print()
    mcfile.Close()

# type(mc["pion_rpa"])
# mc["double_rpa"]=TVectorD(mc["pion_rpa"])
# logmc["double_rpa"] = TVectorD(mc["double_rpa"])
# rpafactor = TVectorD(mc["pion_rpa"])
# for i in range(0,len(goodelements)):
#    rpafactor[i] /= mc["$\pi$tune"][i]
#    mc["double_rpa"][i]*=rpafactor[i]
#    logmc["double_rpa"][i] = math.log(mc["double_rpa"][i])
# rpafactor.Print()

# mc["pion_rpa"].Print()
# mc["double_rpa"].Print()

# models.append("double_rpa")
print(mcdraw)
maxb = len(data)
sysnames = datahist.GetSysErrorMatricesNames()

syscov = {}

onames = {}
# for model in models:
# ofile[model] = open(basemodel+var+"_"+model+".tex",'w')
oname = basemodel + "%s-lowq2.tex" % (var)
ofile["ALL"] = open(oname, "w")
out = ofile["ALL"]
top = "\\documentclass[12pt]{article}\n \\begin{document} \n"
begin = "\\begin{table}\n\\begin{tabular}{lrrr}\n"
out.write(top)
out.write(begin)

InvertedErrorMatrix = {}
LogInvertedErrorMatrix = {}

sysnames = datahist.GetSysErrorMatricesNames()
# print(sysnames)
# cov = datahist.GetSysErrorMatrix("unfoldingCov")
# covx = datahist.ProjectionX().GetSysErrorMatrix("unfoldingCov")
# covx.Print()
sysnames.push_back("ALL")

if NOSYS:
    sysnames = ["ALL"]


for syst in sysnames:
    covmx = compressMatrix(matrix, goodelements)

    if syst != "ALL":
        fullsysmat = readMatrix(datahist.GetSysErrorMatrix(sys), datahist)
    else:
        fullsysmat = readMatrix(datahist.GetTotalErrorMatrix(), datahist)
    syscov[syst] = compressMatrix(fullsysmat, goodelements)
    #
    #  if syst == "unfoldingCov":
    #    print("unfoldingCov")
    #    test = TMatrixD(fullsysmat)
    #    test.ResizeTo(10,10)
    #    test.Print()
    #
    #  if syst != "ALL" and not NOSYS:
    #    covmx -= syscov[syst]

    logcovmx = compressMatrix(matrix, goodelements)

    for i in range(0, len(goodelements)):
        if data[i] == 0:
            continue
        for j in range(0, len(goodelements)):
            if data[j] == 0:
                continue
            logcovmx[i][j] = covmx[i][j] / data[i] / data[j]
    InvertedError = TDecompSVD(covmx)
    InvertedErrorMatrix[syst] = TMatrixD(covmx)
    InvertedError.Invert(InvertedErrorMatrix[syst])

    LogInvertedError = TDecompSVD(logcovmx)
    LogInvertedErrorMatrix[syst] = TMatrixD(logcovmx)
    LogInvertedError.Invert(LogInvertedErrorMatrix[syst])


test = TMatrixD(LogInvertedErrorMatrix["ALL"])
test.ResizeTo(5, 5)
test.Print()

head = "%s &\t %s &\t %s  \\\\ \n" % ("Model", "$\chi^2$ - linear", "$\chi^2$ - log")
print(head)

c1 = TCanvas("c1", "c1_n4" + model, 0, 0, 900, 750)
if not Ratio:
    c1.SetLogy()
c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.15)
c1.SetBottomMargin(0.15)
c1.SetTopMargin(0.15)
if var in ["ptmu", "q2"]:
    c1.SetLogx()

if var in ["q2", "enu", "ptmu", "pzmu"]:
    mcnorm = mcdraw["default"].Clone()
    mcnorm.SetDirectory(0)
    mcdraw["default"].Print()

    datadraw.Print()
    gStyle.SetOptStat(1)
    # data.SetOptStat(0)

    if Ratio:
        datadraw.Divide(mcnorm)
    datadraw.SetDirectory(0)


ofile["ALL"].write(head)
h_chi2cont = {}
n = covmx.GetNrows()
for model in models:
    h_chi2cont[model] = TH2D("cont" + model, "chi2 cont", n, 1, n, n, 1, n)
    if var in ["q2", "enu", "ptmu", "pzmu"]:
        mcdraw[model].SetDirectory(0)
    chi2cont = TMatrixD(covmx)
    chi2cont.Zero()
    for syst in sysnames:
        totchi2 = 0.0
        totlogchi2 = 0.0
        dof = 0
        residual = TVectorD(data)
        logresidual = TVectorD(residual)
        for i in range(0, covmx.GetNrows()):
            residual[i] = data[i] - mc[model][i]
            logresidual[i] = logdata[i] - logmc[model][i]
            # print(i, data[i],mc[model][i],logdata[i],logmc[model][i])
            dof += 1
        chi2diag = TVectorD(residual)
        chi2diag.Zero()
        for i in range(minb, maxb):
            for j in range(minb, maxb):
                totchi2 += residual[i] * InvertedErrorMatrix[syst][i][j] * residual[j]
                totlogchi2 += (
                    logresidual[i] * LogInvertedErrorMatrix[syst][i][j] * logresidual[j]
                )
                chi2cont[i][j] = (
                    residual[i] * InvertedErrorMatrix[syst][i][j] * residual[j]
                )
                h_chi2cont[model].Fill(i, j, chi2cont[i][j])
                if i == j:
                    chi2diag[i] = chi2cont[i][j]

        # print("chi2", model, syst, totchi2, dof, totlogchi2)
        outs = "%s &\t %10.0f &\t %10.0f \\\\\n" % (
            translate[model],
            totchi2,
            totlogchi2,
        )
        print("dump of contributions", model)
        print(outs)
        # chi2diag.Print()
        # chi2cont.Print()
        ofile["ALL"].write(outs)
        h_chi2cont[model].Write()
        print(" before pdf ", syst, var, model)
        if var != "pzmu_ptmu" and syst == "ALL":
            title = "%s ; %s %6.0f" % (model, var, totchi2)

            gStyle.SetTitleAlign(11)
            # datadraw.GetXaxis().SetRangeUser(0.001,10)
            mcdraw[model].Print()
            if Ratio:
                mcdraw[model].Divide(mcnorm)
            mcdraw[model].SetName(title)
            if Ratio:
                mcdraw[model].SetMaximum(2.0)
                mcdraw[model].SetMinimum(0.0)
            else:
                mcdraw[model].SetMaximum(1.0e-36)
            mcdraw[model].Draw("HIST")
            datadraw.Draw("PE same")

            datadraw.Write()
            pname = basemodel + "%s-%s-lowq2_v7.pdf" % (var, model)
            c1.SaveAs(pname, "pdf")
            print(" Try to save ", pname)
            # mchist[model].Delete()
            try:
                mcdraw[model].Delete()
            except:
                print("failed to delete", model)

# ofile[model].close()
#    mchist[model].Delete()
#    mcfile.Close()

end = "\\end{tabular}\n"
out.write(end)
end = "\\caption{%s chi2  for  %d degrees of freedom }\n\\end{table}" % (
    varname[var],
    dof,
)
out.write(end)
out.write("\\end{document}\n")
out.close()

try:
    datadraw.Delete()
except:
    print("failed in data delete")
    sys.exit(1)
