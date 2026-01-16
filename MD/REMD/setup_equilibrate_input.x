#!/bin/bash

rm -f equilibrate.groupfile

nrep=$(wc -l < temperatures.dat)
echo "N replicas = $nrep"

COUNT=0
while read -r TEMP; do
  COUNT=$((COUNT+1))
  REP=$(printf "%03d" "$COUNT")

  echo "TEMPERATURE: $TEMP K ==> FILE: equilibrate.mdin.$REP"

  sed "s/XXXXX/$TEMP/g" equilibrate.mdin \
    | sed "s/RANDOM_NUMBER/$RANDOM/g" \
    > equilibrate.mdin.$REP

  echo "-O -rem 0 -i equilibrate.mdin.$REP \
-o equilibrate.mdout.$REP \
-c min1.rst \
-r equilibrate.rst.$REP \
-x equilibrate.mdcrd.$REP \
-inf equilibrate.mdinfo.$REP \
-p design_1.parm7" >> equilibrate.groupfile

done < temperatures.dat

echo "#" >> equilibrate.groupfile
echo "Done."
