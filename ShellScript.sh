#!/bin/bash
seed=7
nRuns=300
Rwt=2.9
Dres=0
DDR=0.1
mu=0.00001
mu_after_on=1
mig=0
Kmain=10000
Krefu=0
V=1

echo "
				  $seed 
				  $nRuns
				  $Rwt
				  $Dres
				  $DDR
				  $mu
				  $mu_after_on
				  $mig
				  $Kmain
				  $Krefu
				  $V
				  " | ./HIVevolution
