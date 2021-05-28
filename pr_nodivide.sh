#!/bin/bash

################################################################
# INITIAL BOOKKEEPING AND ENVIRONMENT VARIABLES
################################################################
# Author: Connor Masterson
# This script tests the perf stat output from running PageRank
# in both the original Gauss-Seidell implementation and the
# version with no division.

GAPBS_DIR=/soe/ccmaster/gapbs
GRAPH_DIR=/lscratch/graphs
OUT_DIR=$GAPBS_DIR/output/pr_nodivide-results
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
# KRON
################################################################
g='kron'
out_nd="$OUT_DIR/nodivide_$g.out"
if [ -e $out_nd ]; then 
  rm $out_nd
fi
touch $out_nd
out_o="$OUT_DIR/orig_$g.out"
if [ -e $out_o ]; then 
  rm $out_o
fi
touch $out_o
out_nd_r="$OUT_DIR/nodivide_${g}_raw.out"
if [ -e $out_nd_r ]; then 
  rm $out_nd_r
fi
touch $out_nd_r
out_o_r="$OUT_DIR/orig_${g}_raw.out"
if [ -e $out_o_r ]; then 
  rm $out_o_r
fi
touch $out_o_r

echo ""
echo "Starting runs on Kron graphs"
echo "Now doing..."

for scale in 18 19 20 21 22 23 24 25 26 27; do
  echo "SCALE $scale on pr_orig"
  echo "SCALE $scale on pr_orig" >>$out_o
  perf stat ./pr_orig  -g $scale -n 3 2>>$out_o  >>$out_o_r
  echo "SCALE $scale on pr_recip"
  echo "SCALE $scale on pr_recip" >>$out_nd
  perf stat ./pr_recip -g $scale -n 3 2>>$out_nd >>$out_nd_r
done


################################################################
# URAND 
################################################################
g='urand'
out_nd="$OUT_DIR/nodivide_$g.out"
if [ -e $out_nd ]; then 
  rm $out_nd
fi
touch $out_nd
out_o="$OUT_DIR/orig_$g.out"
if [ -e $out_o ]; then 
  rm $out_o
fi
touch $out_o
out_nd_r="$OUT_DIR/nodivide_${g}_raw.out"
if [ -e $out_nd_r ]; then 
  rm $out_nd_r
fi
touch $out_nd_r
out_o_r="$OUT_DIR/orig_${g}_raw.out"
if [ -e $out_o_r ]; then 
  rm $out_o_r
fi
touch $out_o_r

echo ""
echo "Starting runs on Urand graphs"
echo "Now doing..."

for scale in 18 19 20 21 22 23 24 25 26 27; do
  echo "SCALE $scale on pr_orig"
  echo "SCALE $scale on pr_orig" >>$out_o
  perf stat ./pr_orig  -u $scale -n 3 2>>$out_o  >$out_o_r
  echo "SCALE $scale on pr_recip"
  echo "SCALE $scale on pr_recip" >>$out_nd
  perf stat ./pr_recip -u $scale -n 3 2>>$out_nd >$out_nd_r
done

################################################################
# REAL WORLD GRAPHS
################################################################

echo ""
echo "Starting runs on real world graphs"
echo "Now doing..."

for g in twitter road web; do
  echo "$g.sg"
  out_nd="$OUT_DIR/nodivide_$g.out"
  if [ -e $out_nd ]; then 
    rm $out_nd
  fi
  touch $out_nd
  out_o="$OUT_DIR/orig_$g.out"
  if [ -e $out_o ]; then 
    rm $out_o
  fi
  touch $out_o
  out_nd_r="$OUT_DIR/nodivide_${g}_raw.out"
  if [ -e $out_nd_r ]; then 
    rm $out_nd_r
  fi
  touch $out_nd_r
  out_o_r="$OUT_DIR/orig_${g}_raw.out"
  if [ -e $out_o_r ]; then 
  rm $out_o_r
  fi
  touch $out_o_r
  echo "$g on pr_orig"
  echo "$g on pr_orig" >>$out_o
  perf stat ./pr_orig  -f $GRAPH_DIR/$g.sg -n 3 2>>$out_o  >$out_o_r
  echo "$g on pr_recip"
  echo "$g on pr_recip" >>$out_nd
  perf stat ./pr_recip -f $GRAPH_DIR/$g.sg -n 3 2>>$out_nd >$out_nd_r
done
