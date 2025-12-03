#! /bin/bash

explst=""
ofn=""
i=0
for dir in ../Experiments/286*Venus_cloudtop*newedd_newnetS/; do
  if [ -e "${dir}/Output/Photo/atm_composition.dat" ]; then
    explst="$explst ${dir:15:-1}"
    ofn="${ofn}${dir:15:8}"
    i=$((i+1))

    if [ $i -gt 9 ]; then
      python compareExp.py "$explst" "O2 H2O H2 CO SO2 SO OCS HCl NO" -data="Venus/obs2" \
             -addh=58 -ylim="58 120" -height -ofn=$ofn".pdf" -save -sqaxis
      i=0
      explst=""
      ofn=""
    fi
  fi
done

python compareExp.py "$explst" "O2 H2O H2 CO SO2 SO OCS HCl NO" -data="Venus/obs2" \
                     -addh=58 -ylim="58 120" -height -ofn=$ofn".pdf" -save -sqaxis
