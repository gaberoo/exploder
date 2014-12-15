#!/bin/bash

# Calculates the likelihood surface for an interval

minN=20
maxN=400
incN=2

minBeta=0.1
maxBeta=4.0
incBeta=0.05

fixMu=0.1

minPsi=0.01
maxPsi=0.51
incPsi=0.02

rho=0

times=$1

#for N in `seq $minN $incN $maxN`
#do
N=1000
  for beta in `seq $minBeta $incBeta $maxBeta`
  do
    for psi in `seq $minPsi $incPsi $maxPsi`
    do
      ./calc_likelihood -l 0 -N $N -b $beta -R $fixMu -s $psi -o $rho $times
    done
  done
#done

