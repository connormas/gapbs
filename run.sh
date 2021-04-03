################################################################
# INITIAL BOOKKEEPING AND ENVIRONMENT VARIABLES
################################################################
# Author: Connor Masterson
# This script runs both Gauss-Seidell and Jacobi style pagerank
# on Kron and Urand graphs at varying scales.

GAPBS_DIR=/soe/ccmaster/gapbs
GRAPH_DIR=/lscratch/graphs
OUT_DIR=$GAPBS_DIR/output/run.sh-results
#​
# 1 socket - 0-23
# 2 socket - 0-47
#​
THREADS=24
THREAD_LIMIT=$THREADS-1
export OMP_NUM_THREADS=${THREADS}
export GOMP_CPU_AFFINITY="0-23"
PRE="numactl -N 0 -m 0"
# PRE="numactl --interleave=all"

################################################################
# TEST FOR KRON GRAPHS
################################################################

gskout="$OUT_DIR/GS-kron-$THREADS.out"
jkout="$OUT_DIR/J-kron-$THREADS.out"
if [ -e $gskout ]; then
  rm $gskout
fi                                                                              
touch $gskout
if [ -e $jkout ]; then
  rm $jkout
fi                                                                              
touch $jkout

echo "Starting KRON runs scale 17-28"
echo "All results will be in"
echo "$gskout" 
echo "$jkout"

echo "Now doing..."
for SCALE in 17 18 19 20 21 22 23 24 25 26 27 28; do
  echo "kron-$SCALE"  
  echo "SCALE: $SCALE" >> $gskout
  echo "SCALE: $SCALE" >> $jkout
  $GAPBS_DIR/prgs -g $SCALE -n 3 -i 100 | grep -B 2 "Average Time" >> $gskout
  $GAPBS_DIR/prj  -g $SCALE -n 3 -i 100 | grep -B 2 "Average Time" >> $jkout
done

################################################################
# TEST FOR URAND GRAPHS
################################################################

gskout="$OUT_DIR/GS-urand-$THREADS.out"
jkout="$OUT_DIR/J-urand-$THREADS.out"
if [ -e $gskout ]; then
  rm $gskout
fi                                                                              
touch $gskout
if [ -e $jkout ]; then
  rm $jkout
fi                                                                              
touch $jkout

echo "Starting URAND runs scale 17-28"
echo "All results will be in"
echo "$gskout" 
echo "$jkout"

echo "Now doing..."
for SCALE in 17 18 19 20 21 22 23 24 25 26 27 28; do
  echo "urand-$SCALE"  
  echo "SCALE: $SCALE" >> $gskout
  echo "SCALE: $SCALE" >> $jkout
  $GAPBS_DIR/prgs -u $SCALE -n 3 -i 100 | grep -B 2 "Average Time" >> $gskout
  $GAPBS_DIR/prj  -u $SCALE -n 3 -i 100 | grep -B 2 "Average Time" >> $jkout
done

################################################################
# TEST FOR REAL WORLD GRAPHS
################################################################

# bookkeeping
GRAPHS_DIR=/lscratch/graphs
gsout="$OUT_DIR/GS-rwg-$THREADS.out"
jout="$OUT_DIR/J-rwg-$THREADS.out"
if [ -e $gsout ]; then
  rm $gsout
fi                                                                              
touch $gsout
if [ -e $jout ]; then
  rm $jout
fi                                                                              
touch $jout

echo "Now doing real world graphs..."
for graph in road.sg twitter.sg web.sg; do
  echo "$GRAPHS_DIR/$graph"
  echo "$graph" >> $gsout
  echo "$graph" >> $jout 
  $GAPBS_DIR/prgs -f $GRAPHS_DIR/$graph -n 3 -i 100 | grep -B 2 "Average Time" >> $gsout
  $GAPBS_DIR/prj  -f $GRAPHS_DIR/$graph -n 3 -i 100 | grep -B 2 "Average Time" >> $jout
done
