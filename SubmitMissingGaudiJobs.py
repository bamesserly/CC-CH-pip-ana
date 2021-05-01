################################################################################
# Search an output ProcessAna dir for missing subrun rootfiles and resubmit
################################################################################

import os, sys

plists = ["1A","1B","1C","1D","1E","1F","1G","1L","1M","1N","1O","1P"]
ignore_plists = ["1E","1F","1G","1L","1M","1N","1O","1P"]
data_first_pass = (
    "/pnfs/minerva/persistent/users/bmesserl/pions/20210307/data/ME{plist}/grid/minerva/ana/numibeam/"
)
mc_first_pass = (
    "/pnfs/minerva/persistent/users/bmesserl/pions/20210307/mc/ME{plist}/grid/central_value/minerva/ana/"
)
second_pass = (
    "/pnfs/minerva/persistent/users/bmesserl/pions/20210430/"
)
tarball = "/pnfs/minerva/resilient/tarballs/bmesserl_v22r1p1-20210307_154853.tgz"
dry_run = False
do_data = True
do_mc = True

# Get the run/subrun from a filename.
# Output format [run, subrun]
def FindRunSubFromFilename(f, ftype, isData):
    # if mc 2,4 data 1,3
    run_spot = 2
    sub_spot = 4
    if isData:
        run_spot = 1
        sub_spot = 3
    filename = os.path.basename(f.rstrip("\n")).split("_")
    if ftype == "reco":
        # if mc 2,3 data 1,2
        run_spot = 2
        sub_spot = 3
        if isData:
            run_spot = 1
            sub_spot = 2
    run = int(filename[run_spot])
    subrun = int(filename[sub_spot])
    return [run, subrun]


# Combine missing run/subruns into format processana can use.
# Output format [run, 'subrun1, subrun2, subrun3,']
def GetRunSub(miss):
    myrunsubs = []
    found_run = []
    for i, el in enumerate(miss):
        r = el[0]
        s = el[1]
        if r in found_run:
            continue
        found_run.append(r)
        mysubruns = "%d," % s
        for el2 in miss[i + 1 :]:
            if el2[0] == r:
                mysubruns += "%d," % (el2[1])
        myrunsubs.append([r, mysubruns])
    return myrunsubs


# Call ProcessAna for certain runs/subruns.
def SubmitJobs(miss, pl, out_dir_base, isData):
    myrunssubruns = GetRunSub(miss)
    for r in myrunssubruns:
        myrun = r[0]
        mysubruns = r[1].rstrip(",")  # last element as a comma
        cmd = (
            "python $PRODUCTIONSCRIPTSLITEROOT/ana_scripts/ProcessAna.py --ana_tool CCNuPionInc "
            "--inv Inextinguishable --kludge NX --no_verify_kludge "
            "--usecat "
            "--lifetime 12h "
            "{TARBALL} "
            "--run {RUN} "
            "--subrun {SUBRUN} "
            "--outdir {BASEDIR}/{DATAMC}/ME{PL}_{RUN} "
            "--{DATAMC} "
            "{TRACKER} "
        ).format(
            TARBALL = "--tarball {0}".format(tarball) if tarball else "",
            RUN = myrun,
            SUBRUN = mysubruns,
            BASEDIR = out_dir_base,
            DATAMC = "data" if isData else "mc",
            PL = pl,
            TRACKER = "" if isData else "--tracker"
        )

        print pl, myrun, mysubruns
        #print(cmd)
        if not dry_run:
          os.system(cmd)


if __name__ == "__main__":
  if do_data:
    print("DATA")
    for l in plists:
        if l in ignore_plists:
            continue
        isData = True
        ana_run_sub = []
        mas_run_sub = []
        data_first_pass_pl = data_first_pass.format(plist=l)
        analyzed_list = os.popen("find %s -type f -name \"*.root\"|sort" % (data_first_pass_pl)).readlines()
        master_list = os.popen(
            "samweb list-definition-files rodriges_data_reco_minervame%s_inextinguishable"
            % (l)
        )

        mas_run_sub_samname = {}
        for f in analyzed_list:
            s = FindRunSubFromFilename(f, "ana", isData)
            ana_run_sub.append(s)
        for f in master_list:
            s = FindRunSubFromFilename(f, "reco", isData)
            mas_run_sub.append(s)
            mas_run_sub_samname[tuple(s)] = f.rstrip("\n")

        missing = []
        missing_samname = []
        for el in mas_run_sub:
            if el not in ana_run_sub:
                missing.append(el)
                missing_samname.append(mas_run_sub_samname[tuple(el)])

        #print(",".join(missing_samname))
        SubmitJobs(missing, l, second_pass, isData)
        print("Data ME%s: %d subruns resubmitted" % (l, len(missing)))

  if do_mc:
    print("Tracker")
    for l in plists:
        if l in ignore_plists:
            continue
        isData = False
        ana_run_sub = []
        mas_run_sub = []
        mc_first_pass_pl = mc_first_pass.format(plist=l)
        analyzed_list = os.popen("find %s -type f -name \"*.root\" |sort" % (mc_first_pass_pl)).readlines()
        master_list = os.popen(
            "samweb list-definition-files rodriges_mc_reco_minervame%s_inextinguishable_tracker"
            % (l)
        )

        mas_run_sub_samname = {}
        for f in analyzed_list:
            s = FindRunSubFromFilename(f, "ana", isData)
            ana_run_sub.append(s)
        for f in master_list:
            s = FindRunSubFromFilename(f, "reco", isData)
            mas_run_sub.append(s)
            mas_run_sub_samname[tuple(s)] = f.rstrip("\n")

        missing = []
        missing_samname = []
        for el in mas_run_sub:
            if el not in ana_run_sub:
                missing.append(el)
                missing_samname.append(mas_run_sub_samname[tuple(el)])
        #print(",".join(missing_samname))
        SubmitJobs(missing, l, second_pass, isData)
        print("MC ME%s: %d subruns resubmitted" % (l, len(missing)))

  sys.exit()
