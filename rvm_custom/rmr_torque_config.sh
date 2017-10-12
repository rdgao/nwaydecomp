  1 #!/bin/sh
  2 #
  3 #PBS -q hotel
  4 #PBS -N upennram
  5 #PBS -l nodes=1
  6 #PBS -l walltime=150:00:00
  7 #PBS -l mem=25gb
  8 #PBS -o upennram.o
  9 #PBS -e upennram.e
 10 #PBS -V
 11 #PBS -M rvandermeij@ucsd.edu
 12 #PBS -m ae
 13 cd ~
 14 /opt/matlab/2015a/bin/matlab -singleCompThread -nodesktop -nojvm -r rmr_upennram_getspace
 15
 16 # get out of here
 17 exit