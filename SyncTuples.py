################################################################################
# Consolidate tuples from multiple processing attempts
#
# Technically: copy tuples from a given second pass location into a first pass
# location
################################################################################
import os
import sys
from itertools import product
plists = ["1A","1B","1C","1D","1E","1F","1G","1L","1M","1N","1O","1P"]
ignore_plists = ["1E","1F","1G","1L","1M","1N","1O","1P"]
second_pass_fmt = "/pnfs/minerva/persistent/users/bmesserl/pions/20210417/{DATAMC}/"
first_pass_fmt = "/pnfs/minerva/persistent/users/bmesserl/pions/20210307/{DATAMC}/ME{PLIST}/"

for dmc in ["data", "mc"]:
    print(dmc.upper())
    second_pass_base_dir = second_pass_fmt.format(DATAMC = dmc)
    addeddirs = os.popen("ls %s"%(second_pass_base_dir)).readlines()
    for d,p in product(addeddirs, plists):
        if p in ignore_plists:
          continue
        if p not in d:
          continue
        second_pass = second_pass_base_dir + d.rstrip("\n") + "/"
        first_pass = first_pass_fmt.format(DATAMC = dmc, PLIST = p)
        cmd = (
            "rsync -arvv --progress "
            "--ignore-existing "
            "--remove-source-files "
            "--exclude 'logfiles/' "
            "--exclude 'opts/' "
            "--exclude 'opts_in/' "
            "--exclude 'pot/' "
            "--exclude 'hist/' "
            "--exclude '*.tgz' "
            "{SOURCE} {TARGET}".format(
            SOURCE = second_pass,
            TARGET = first_pass
            )
        )
        print(cmd)
        os.system(cmd)
