#!/bin/bash

################################################################
# INITIAL BOOKKEEPING AND ENVIRONMENT VARIABLES
################################################################
# Author: Connor Masterson
# This script tests the effects of differnt loop chunk sizes
# in the inner loop of Gauss-Seidell style pagerank.

GAPBS_DIR=/soe/ccmaster/gapbs
GRAPH_DIR=/lscratch/graphs
OUT_DIR=$GAPBS_DIR/output/loopchunk.sh-results
#​
# 1 socket - 0-23
# 2 socket - 0-47
#​
THREADS=24
# THREADS=48
THREAD_LIMIT=$THREADS-1
export OMP_NUM_THREADS=${THREADS}
export GOMP_CPU_AFFINITY="0-23"
# export GOMP_CPU_AFFINITY="0-47"
PRE="numactl -N 0 -m 0"
# PRE="numactl --interleave=all"

################################################################
# TESTS FOR REAL WORLD GRAPHS 
################################################################
gsout="$OUT_DIR/GS-rw-$THREADS.out"
if [ -e $gskout ]; then
  rm $gsout
fi                                                                              
touch $gsout

echo "Starting runs on Real-World Graphs"
echo "All results will be in"
echo "$gsout"
echo ""
echo "Now doing..."

for graph in twitter.sg road.sg web.sg; do
  for cs in 4 8 16 32 64 128 256 512 1024 2048 8192 16384 32768 65536 131072 262144 524288 1048576; do 
    echo "GRAPH: $graph SIZE: $cs"
    echo "GRAPH: $graph SIZE: $cs" >> $gsout
    g++ -fopenmp -Dchunksize=$cs -std=c++11 -O3 -Wall src/pr.cc -o prgss
    $GAPBS_DIR/prgss -f $GRAPH_DIR/$graph -n 3 -i 100 | grep -B 2 "Average Time" >> $gsout
  done
  echo "" >> $gsout
done

