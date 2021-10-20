echo "This program shows the number of logfiles that have a Exit status != 0"
echo "you can watch which files have the Exit status != 0 in outFilewithAfiledlogFiles.txt" 

declare -a ME1Amclist=("00" "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31" "32" "33" "34" "35" "36" "37" "38" "39" "40")

declare -a ME1AUpdatalist=("ME1Adatalist60" "ME1Adatalist61" "ME1Adatalist62" "ME1Adatalist80" "ME1Adatalist01")
declare -a ME1Adatalist60=("38" "39" "40" "41" "42" "43" "44" "45" "46" "47" "48" "51" "52" "57" "58" "59" "60" "61" "62" "63" "64" "65" "66" "67" "68" "69" "70" "71" "72" "73" "74" "75" "76" "77" "78" "79" "80" "81" "82" "83" "84" "85" "86" "87" "88" "89" "90" "91" "92" "93" "94" "95")
declare -a ME1Adatalist61=("12" "15" "18" "21" "30" "33" "36" "39" "42" "45" "53" "56" "59" "67" "70" "73" "76" "82" "85" "88" "91" "94" "97" "13" "16" "19" "28" "31" "34" "37" "40" "43" "51" "54" "57" "60" "68" "71" "74" "77" "83" "86" "89" "92" "95" "98" "14" "17" "20" "29" "32" "35" "38" "41" "44" "52" "55" "58" "66" "69" "72" "75" "78" "84" "87" "90" "93" "96" "99")
declare -a ME1Adatalist62=("00" "02" "04" "06" "08" "10" "12" "18" "21" "23" "25" "27" "29" "31" "33" "35" "37" "39" "41" "44" "46" "48" "50" "52" "54" "56" "01" "03" "05" "07" "09" "11" "13" "20" "22" "24" "26" "28" "30" "32" "34" "36" "38" "40" "43" "45" "47" "49" "51" "53" "55")
declare -a ME1Adatalist80=("00" "01" "02" "03" "04" "05" "06" "07" "11" "12" "13" "14" "15" "29" "30" "31" "32" "33" "34")
declare -a ME1Adatalist01=("00" "03" "06" "10" "13" "16" "19" "22" "25" "28" "34" "37" "40" "43" "46" "49" "52" "56" "59" "62" "65" "01" "04" "07" "11" "14" "17" "20" "23" "26" "31" "35" "38" "41" "44" "47" "50" "53" "57" "60" "63" "66" "02" "05" "08" "12" "15" "18" "21" "24" "27" "32" "36" "39" "42" "45" "48" "51" "54" "58" "61" "64")

declare -a ME1Bmclist1110=("00" "01" "02" "03" "04" "05" "06" "07" "08" "09" "10")
declare -a ME1Bdatalist00=("68" "69" "70" "71" "72" "73" "74" "75" "76" "77" "79" "80" "81" "82" "83" "84" "85" "86" "88" "89" "91" "95" "96" "97" "98" "99")
declare -a ME1Bdatalist01=("00" "01" "06" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "28")

#indir="/pnfs/minerva/persistent/users/granados/data_MADME1A_v41/grid/minerva/logfiles/numibeam/v22r1p1/00/00/60/%s/"
#indir="/pnfs/minerva/persistent/users/granados/mc_MADME1B_v41p1/grid/minerva/logfiles/numibeam/v22r1p1/00/11/10/%s/"
indir="/pnfs/minerva/persistent/users/granados/mc_MADME1B_v49/grid/central_value/minerva/logfiles/v22r1p1/00/11/10/%s/"
outfile="textFiles/outFilewithAfiledlogFiles_%s.txt"

rm textFiles/*.txt
rm FailList.txt

for list in "${ME1Bmclist1110[@]}"; do
#  declare -a aux=$list
#  for list1 in "${aux[@]}"; do
#    echo " printing ${list} list ${list1}"
#  done
  infile=$(printf "$indir" $list)
#  rm ${outfile}
  outputfile=$(printf "$outfile" $list)
  zgrep -ce "Exit Status: 1" -ce "Exit Status: 2" -ce "Exit Status: 3" -ce "Exit Status: 4" -ce "Exit Status: 5" -ce "Exit Status:6" ${infile}*.log.gz > ${outputfile}
  echo "Number of fail files in 11/10/${list}:"
  grep -c "log.gz:1" ${outputfile} 
  echo "Number of correct files in 11/10/${list}:"
  grep -c "log.gz:0" ${outputfile}
done
grep -n "log.gz:1" textFiles/*.txt > FailList.txt 
