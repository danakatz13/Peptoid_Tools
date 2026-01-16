#!/usr/bin/env bash
set -euo pipefail

TFILE="temperatures.dat"
TEMPLATE="prod.in"
GROUP="remd.groupfile"

rm -f "$GROUP"

nrep=$(wc -l < "$TFILE")
echo "N REPLICAS = $nrep"

COUNT=0
while read -r TEMP; do
  [[ -z "${TEMP// }" ]] && continue
  [[ "$TEMP" =~ ^# ]] && continue

  COUNT=$((COUNT+1))
  REP=$(printf "%03d" "$COUNT")

  OUTMDIN="remd.mdin.$REP"
  echo "TEMPERATURE: $TEMP K ==> FILE: $OUTMDIN"

  sed -e "s/XXXXX/${TEMP}/g" -e "s/RANDOM_NUMBER/${RANDOM}/g" "$TEMPLATE" > "$OUTMDIN"

  if ! grep -q "$TEMP" "$OUTMDIN"; then
    echo "WARNING: '$TEMP' not found in $OUTMDIN. Check your token in $TEMPLATE." >&2
  fi

  echo "-O -rem 1 -remlog rem.log -i $OUTMDIN -o prod.out.$REP -c equilibrate.rst.$REP -r remd.rst.$REP -x remd.mdcrd.$REP -inf remd.mdinfo.$REP -p __PARM7__" >> "$GROUP"
done < "$TFILE"

echo "#" >> "$GROUP"
echo "Done."

