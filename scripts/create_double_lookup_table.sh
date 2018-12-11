#!/usr/bin/env bash

echo "double double_logsum_lookup[] = {" > ../src/double_logsum_lookup.h
  for i in {0..36044}
  do
      echo " l( 1.0 + e(-1 * $i / 1000))" | bc -l | awk -v num=$i '{\
  if(num == 36044){\
  printf "%.*f\
  };\n",14,$0}else{\
  printf "%.*f,\n",14,$0\
  }}'>> ../src/double_logsum_lookup.h
  done
