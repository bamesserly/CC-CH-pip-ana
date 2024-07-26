import os, subprocess,re, sys,glob
from collections import OrderedDict,namedtuple

outdir = "/exp/minerva/app/users/granados/cmtuser/MATAna/cc-ch-pip-ana/GiBUU"

modeldir = "/pnfs/minerva/persistent/Models/GiBUU/Medium_Energy/FHC/v2021"

filetuple = namedtuple("filetuple","filelist outname")

tarnames=OrderedDict()
#tarnames["Nuclear_C"]  = "carbon"
#tarnames["Nuclear_Fe"] = "iron"
#tarnames["Nuclear_O"]  = "oxygen"
#tarnames["Nuclear_Pb"] = "lead"
tarnames["Tracker_Carbon"]  = "tracker_carbon"
tarnames["Tracker_C"]  = "tracker_carbon"

tarmodels=OrderedDict()
tarmodels["tracker"] = ["Tracker_Carbon_T0", "Tracker_C_T1", "Tracker_Hydrogen"]
#tarmodels["tracker"] = [ "Tracker_C_T1"]

#tarmodels["nuclear"] = ["Nuclear_C_T0",
#"Nuclear_C_T1",
#"Nuclear_Fe_T0",
#"Nuclear_Fe_T1",
#"Nuclear_O_T0",
#"Nuclear_O_T1",
#"Nuclear_Pb_T0",
#"Nuclear_Pb_T1"]

#tarmodels["nuclear"] = ["Nuclear_O_T0","Nuclear_O_T1"]

filelist=[]
outnames = []

#NukeCCPion_<Generator>_<Model>_<Target>

for target in tarmodels:
  for tarmodel in tarmodels[target]:
    if tarmodel.find("Hydrogen")>-1:
      outname = "NukeCCPion_GiBUU_hydrogen"
    regex = re.search("(.*)_T(.*)",tarmodel)
    if regex:
      outname = "NukeCCPion_GiBUU_T{1}_{0}".format(tarnames[regex.group(1)],regex.group(2))
    outnames.append(outname)
    print (outname)

#    print "{0}/{1}/{2}".format(modeldir,target,tarmodel)
#    tmpfiles = [] 
#    for root, dirs, files in os.walk("{0}/{1}/{2}".format(modeldir,target,tarmodel)):
#      for name in files:
#        #print root,name
#        if name.find(".dat")==-1:
#          continue
#        tmpfiles.append("{0}/{1}".format(root,name))
#
#    with open("GiBUUDat/"+outname+".txt","w") as infile:
#      for line in tmpfiles:
#        infile.write(line+"\n")

template_command = "python xsec/MakeGiBUUModelXsec.py xsec/GiBUUDat/{0}.txt {1}/{0}.root -b"
for outname in outnames:
  model_command = template_command.format(outname,outdir)
  subprocess.call(model_command,shell=True)
  
  
#model_command = template_command.format("testList",outdir)
#subprocess.call(model_command,shell=True)


#
#infile = "/pnfs/minerva/persistent/Models/GENIE/Medium_Energy/FHC/v3_0_6/tracker/G18_10a_02_11a/flat_GENIE_G18_10a_02_11a_50M.root"
#outfile = "{0}/NukeCCPion_{1}_{2}_{3}.root".format(outdir,"GENIE","G18_10a_02_11a","tracker")
#
#model_command = template_command.format(infile=infile,tartype="tracker",outfile=outfile)
#subprocess.call(model_command,shell=True)

