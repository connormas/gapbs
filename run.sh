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

# make output files
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
for SCALE in 17 18 19 20  21 22 23 24 25 26 27 28; do
  
  echo "kron-$SCALE"
  
  echo "SCALE: $SCALE" >> $gskout
  echo "SCALE: $SCALE" >> $jkout
  $GAPBS_DIR/prgs -g $SCALE -i 5 | grep -B 2 "Average Time" >> $gskout
  $GAPBS_DIR/prj  -g $SCALE -i 5 | grep -B 2 "Average Time" >> $jkout

done

