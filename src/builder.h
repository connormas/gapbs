// copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef BUILDER_H_
#define BUILDER_H_

#include <algorithm>
#include <cinttypes>
#include <fstream>
#include <functional>
#include <type_traits>
#include <utility>

#include "command_line.h"
#include "generator.h"
#include "graph.h"
#include "platform_atomics.h"
#include "pvector.h"
#include "reader.h"
#include "timer.h"
#include "util.h"

/*
GAP Benchmark Suite
Class:  BuilderBase
Author: Scott Beamer

Given arguements from the command line (cli), returns a built graph
 - MakeGraph() will parse cli and obtain edgelist and call
   MakeGraphFromEL(edgelist) to perform actual graph construction
 - edgelist can be from file (reader) or synthetically generated (generator)
 - Common case: BuilderBase typedef'd (w/ params) to be Builder (benchmark.h)
*/


template <typename NodeID_, typename DestID_ = NodeID_,
          typename WeightT_ = NodeID_, bool invert = true>
class BuilderBase {
  typedef EdgePair<NodeID_, DestID_> Edge;
  typedef pvector<Edge> EdgeList;

  const CLBase &cli_;
  bool symmetrize_;
  bool needs_weights_;
  bool inPlace_ = false;
  int64_t num_nodes_ = -1;


 public:
  explicit BuilderBase(const CLBase &cli) : cli_(cli) {
    symmetrize_ = cli_.symmetrize();
    needs_weights_ = !std::is_same<NodeID_, DestID_>::value;
  }

  DestID_ GetSource(EdgePair<NodeID_, NodeID_> e) {
    return e.u;
  }

  DestID_ GetSource(EdgePair<NodeID_, NodeWeight<NodeID_, WeightT_>> e) {
    return NodeWeight<NodeID_, WeightT_>(e.u, e.v.w);
  }

  NodeID_ FindMaxNodeID(const EdgeList &el) {
    NodeID_ max_seen = 0;
    #pragma omp parallel for reduction(max : max_seen)
    for (auto it = el.begin(); it < el.end(); it++) {
      Edge e = *it;
      max_seen = std::max(max_seen, e.u);
      max_seen = std::max(max_seen, (NodeID_) e.v);
    }
    return max_seen;
  }

  pvector<NodeID_> CountDegrees(const EdgeList &el, bool transpose) {
    pvector<NodeID_> degrees(num_nodes_, 0);
    #pragma omp parallel for
    for (auto it = el.begin(); it < el.end(); it++) {
      Edge e = *it;
      if (symmetrize_ || (!symmetrize_ && !transpose))
        fetch_and_add(degrees[e.u], 1);
      if (symmetrize_ || (!symmetrize_ && transpose))
        fetch_and_add(degrees[(NodeID_) e.v], 1);
    }
    return degrees;
  }

  static
  pvector<SGOffset> PrefixSum(const pvector<NodeID_> &degrees) {
    pvector<SGOffset> sums(degrees.size() + 1);
    SGOffset total = 0;
    for (size_t n=0; n < degrees.size(); n++) {
      sums[n] = total;
      total += degrees[n];
    }
    sums[degrees.size()] = total;
    return sums;
  }

  static
  pvector<SGOffset> ParallelPrefixSum(const pvector<NodeID_> &degrees) {
    const size_t block_size = 1<<20;
    const size_t num_blocks = (degrees.size() + block_size - 1) / block_size;
    pvector<SGOffset> local_sums(num_blocks);
    #pragma omp parallel for
    for (size_t block=0; block < num_blocks; block++) {
      SGOffset lsum = 0;
      size_t block_end = std::min((block + 1) * block_size, degrees.size());
      for (size_t i=block * block_size; i < block_end; i++)
        lsum += degrees[i];
      local_sums[block] = lsum;
    }
    pvector<SGOffset> bulk_prefix(num_blocks+1);
    SGOffset total = 0;
    for (size_t block=0; block < num_blocks; block++) {
      bulk_prefix[block] = total;
      total += local_sums[block];
    }
    bulk_prefix[num_blocks] = total;
    pvector<SGOffset> prefix(degrees.size() + 1);
    #pragma omp parallel for
    for (size_t block=0; block < num_blocks; block++) {
      SGOffset local_total = bulk_prefix[block];
      size_t block_end = std::min((block + 1) * block_size, degrees.size());
      for (size_t i=block * block_size; i < block_end; i++) {
        prefix[i] = local_total;
        local_total += degrees[i];
      }
    }
    prefix[degrees.size()] = bulk_prefix[num_blocks];
    return prefix;
  }

  // Removes self-loops and redundant edges
  // Side effect: neighbor IDs will be sorted
  void SquishCSR(const CSRGraph<NodeID_, DestID_, invert> &g, bool transpose,
                 DestID_*** sq_index, DestID_** sq_neighs) {
    pvector<NodeID_> diffs(g.num_nodes());
    DestID_ *n_start, *n_end;
    #pragma omp parallel for private(n_start, n_end)
    for (NodeID_ n=0; n < g.num_nodes(); n++) {
      if (transpose) {
        n_start = g.in_neigh(n).begin();
        n_end = g.in_neigh(n).end();
      } else {
        n_start = g.out_neigh(n).begin();
        n_end = g.out_neigh(n).end();
      }
      std::sort(n_start, n_end);
      DestID_ *new_end = std::unique(n_start, n_end);
      new_end = std::remove(n_start, new_end, n);
      diffs[n] = new_end - n_start;
    }
    pvector<SGOffset> sq_offsets = ParallelPrefixSum(diffs);
    *sq_neighs = new DestID_[sq_offsets[g.num_nodes()]];
    *sq_index = CSRGraph<NodeID_, DestID_>::GenIndex(sq_offsets, *sq_neighs);
    #pragma omp parallel for private(n_start)
    for (NodeID_ n=0; n < g.num_nodes(); n++) {
      if (transpose)
        n_start = g.in_neigh(n).begin();
      else
        n_start = g.out_neigh(n).begin();
      std::copy(n_start, n_start+diffs[n], (*sq_index)[n]);
    }
  }

  CSRGraph<NodeID_, DestID_, invert> SquishGraph(
      const CSRGraph<NodeID_, DestID_, invert> &g) {
    DestID_ **out_index, *out_neighs, **in_index, *in_neighs;
    SquishCSR(g, false, &out_index, &out_neighs);
    if (g.directed()) {
      if (invert)
        SquishCSR(g, true, &in_index, &in_neighs);

      return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), out_index,
                                                out_neighs, in_index,
                                                in_neighs);
    } else {
      return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), out_index,
                                                out_neighs);
    }
  }

  /*
  This function takes an EdgeList Obj by reference. It assumes this is an edgelist that
  has been partially overwritten with the outneighs of the original EdgeList data. This code
  makes and returns a pvector called inoffsets that will contain the offsets for writing the
  inneighs.
  NOTE: not currently being used, may be cut out eventually
  */
  pvector<SGOffset> MakeOffsetsFromOutNeighs(const EdgeList &el, int elLength){
    pvector<SGOffset> inoffsets(num_nodes_ + 1, 0);
    DestID_* N = (DestID_*)el.data();
    int total = 0;
    for(int i = 0; i < (int)elLength; i++, N++){
      inoffsets[*N + 1]++;
    }
    for(int i = 0; i < (int)inoffsets.size(); i++){
      total += inoffsets[i];
      inoffsets[i] = total;
    }
    return inoffsets;
  }

  /*
  Code that attempts to build in-place and keep peak mem
  usage low. We do this by repurposing the EdgeList pvector and
  using it as the new neighbors array.
  */
  void MakeCSRInPlace(const EdgeList &el, bool transpose, DestID_*** index, DestID_** neighs,
                      DestID_*** inv_index, DestID_** inv_neighs){

    // VARIABLE/OBJECT DECLARATIONS
    std::sort(el.begin(), el.end());
    pvector<NodeID_> degrees = CountDegrees(el, false);
    pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
    pvector<NodeID_> indegrees = CountDegrees(el, true);
    //pvector<SGOffset> inoffsets = ParallelPrefixSum(indegrees);
    DestID_* overWriteEL = (DestID_*)el.data();
    *neighs = (DestID_*)el.data();
    int elLength = el.size();
    *inv_neighs = overWriteEL;
    // OUT GOING NEIGHBORS
    for(int i = 0; i < offsets.size(); i++){
      std::cout << offsets[i] << " ";
    }
    std::cout << std::endl;
    for(auto it = el.begin(); it < el.end(); it++){
      Edge e = *it;
      std::cout << "(" << e.u << ", " << e.v << ") ";
    }
    std::cout << std::endl;

    //overwrite EdgeList memory
    for(auto it = el.begin(); it < el.end(); it++){
      Edge e = *it;
      auto ev = e.v;
      /*for(int i = 0; i < offsets.size(); i++){
        std::cout << offsets[i] << " ";
      }
      std::cout << "(" << e.u << ", " << e.v << ") <-- ";
      for(auto itt = el.begin(); itt < el.end(); itt++){
        Edge ee = *itt;
        std::cout << "(" << ee.u << ", " << ee.v << ") ";
      }
      std::cout << std::endl;*/

      if (symmetrize_ || (!symmetrize_ && !transpose)){
        (*neighs)[fetch_and_add(offsets[e.u], 1)] = ev;
      }
    }

    //printing for debugging
    for(int i = 0; i < offsets.size(); i++){
      std::cout << offsets[i] << " ";
    }
    std::cout << std::endl;
    for(auto it = el.begin(); it < el.end(); it++){
      Edge e = *it;
      std::cout << "(" << e.u << ", " << e.v << ") ";
    }
    std::cout << std::endl;

    //revert offsets
    std::cout << std::endl;
    for(int i = offsets.size(); i >= 0; i--){
        offsets[i] = i != 0 ? offsets[i-1] : offsets[i] = 0;
    }

    // realloc to proper size
    // make sure can handle weighted edges too
    if(!symmetrize_)
      *neighs = (DestID_*)std::realloc((DestID_*)el.data(), (elLength * sizeof(DestID_)));
    else
      *neighs = (DestID_*)std::realloc((DestID_*)el.data(), 2 * (elLength * sizeof(DestID_)));

    std::cout << "calling GenIndex for offsets/neighs\n";
    for(int i = 0; i < degrees.size(); i++){
      std::cout << degrees[i] << " ";
    }
    std::cout << std::endl;
    *index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, *neighs);
    pvector<SGOffset> inoffsets = ParallelPrefixSum(indegrees);
    if((symmetrize_ || (!symmetrize_ && transpose))){//(!symmetrize_){
      std::cout << "calling GenIndex for inoffsets/inv_neighs\n";
      *inv_index = CSRGraph<NodeID_, DestID_>::GenIndex(inoffsets, *inv_neighs);
    }

    // INCOMING NEIGHBORS
    *inv_neighs = new DestID_[inoffsets[num_nodes_]];
    if (!symmetrize_) {
      // write in-neighs to new malloc'd memory
      *inv_index = CSRGraph<NodeID_, DestID_>::GenIndex(inoffsets, *inv_neighs);
      auto deg = degrees.data();
      DestID_* N = (DestID_*)el.data();
      int neighbor = 0;
      for(int i = 0; i < (int)degrees.size(); i++, deg++, neighbor++){
        for(int j = 0; j < (int)(*deg); j++, N++){
          (*inv_neighs)[fetch_and_add(inoffsets[*N], 1)] = neighbor;
        }
      }
    }
    // printing out final result should be outneighs then inneighs
    // will be removed later
    /*DestID_* n = (DestID_*)el.data();
    for(int i = 0; i < (int)elLength; i++, n++) {
      std::cout << "MakeCSRInPlace " << i << ": " << *n << "\n";
    }
    if(!symmetrize_){
      n = inv_neighs[0];
      for(int i = 0; i < (int)elLength; i++, n++) {
        std::cout << "MakeCSRInPlace " << i << ": " << *n << "\n";
      }
    }*/
  }

  /*
  Graph Bulding Steps (for CSR):
    - Read edgelist once to determine vertex degrees (CountDegrees)
    - Determine vertex offsets by a prefix sum (ParallelPrefixSum)
    - Allocate storage and set points according to offsets (GenIndex)
    - Copy edges into storage
  */
  void MakeCSR(const EdgeList &el, bool transpose, DestID_*** index,
               DestID_** neighs) {
    std::sort(el.begin(), el.end());
    pvector<NodeID_> degrees = CountDegrees(el, transpose);
    pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
    *neighs = new DestID_[offsets[num_nodes_]];
    *index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, *neighs);
    #pragma omp parallel for
    for (auto it = el.begin(); it < el.end(); it++) {
      Edge e = *it;
      if (symmetrize_ || (!symmetrize_ && !transpose)){
        (*neighs)[fetch_and_add(offsets[e.u], 1)] = e.v;
      }
      if (symmetrize_ || (!symmetrize_ && transpose)){
        (*neighs)[fetch_and_add(offsets[static_cast<NodeID_>(e.v)], 1)] =
            GetSource(e);
      }
    }
    // print back out the neighbors we just wrote
    for(int i = 0; i < (int)el.size(); i++){
      std::cout << "MakeCSR " << i << " " << ((*neighs)[i]) << "\n";
    }
  }

  CSRGraph<NodeID_, DestID_, invert> MakeGraphFromEL(EdgeList &el) {
    DestID_ **index = nullptr, **inv_index = nullptr;
    DestID_ *neighs = nullptr, *inv_neighs = nullptr;
    Timer t;
    t.Start();
    if (num_nodes_ == -1)
      num_nodes_ = FindMaxNodeID(el)+1;
    if (needs_weights_)
      Generator<NodeID_, DestID_, WeightT_>::InsertWeights(el);
    if (inPlace_){
      MakeCSRInPlace(el, false, &index, &neighs, &inv_index, &inv_neighs);
    } else {
      MakeCSR(el, false, &index, &neighs);
      if (!symmetrize_ && invert){
        MakeCSR(el, true, &inv_index, &inv_neighs);
      }
    }
    t.Stop();
    PrintTime("Build Time", t.Seconds());
    if (symmetrize_){
      std::cout << "returning graph w/ out inv index/neighs\n";
      return CSRGraph<NodeID_, DestID_, invert>(num_nodes_, index, neighs);
    } else {
      std::cout << "returning graph with inv index/neighs\n";
      return CSRGraph<NodeID_, DestID_, invert>(num_nodes_, index, neighs,
                                                inv_index, inv_neighs);
    }
  }

  CSRGraph<NodeID_, DestID_, invert> MakeGraph(bool mFlag=false) {
    CSRGraph<NodeID_, DestID_, invert> g;
    {  // extra scope to trigger earlier deletion of el (save memory)
      inPlace_ = mFlag;
      EdgeList el;
      if (cli_.filename() != "") {
        Reader<NodeID_, DestID_, WeightT_, invert> r(cli_.filename());
        if ((r.GetSuffix() == ".sg") || (r.GetSuffix() == ".wsg")) {
          return r.ReadSerializedGraph();
        } else {
          el = r.ReadFile(needs_weights_);
        }
      } else if (cli_.scale() != -1) {
        Generator<NodeID_, DestID_> gen(cli_.scale(), cli_.degree());
        el = gen.GenerateEL(cli_.uniform());
      }
      std::cout << "calling MakeGraphFromEL\n";
      g = MakeGraphFromEL(el);
      std::cout << "done w/ MakeGraphFromEL and ABOUT TO FREE EL MEMORY\n";
    }
    //can be taken out later, or modify squishgraph()
    std::cout << "calling squishgraph\n";
    return SquishGraph(g);
  }

  // Relabels (and rebuilds) graph by order of decreasing degree
  static
  CSRGraph<NodeID_, DestID_, invert> RelabelByDegree(
      const CSRGraph<NodeID_, DestID_, invert> &g) {
    if (g.directed()) {
      std::cout << "Cannot relabel directed graph" << std::endl;
      std::exit(-11);
    }
    Timer t;
    t.Start();
    typedef std::pair<int64_t, NodeID_> degree_node_p;
    pvector<degree_node_p> degree_id_pairs(g.num_nodes());
    #pragma omp parallel for
    for (NodeID_ n=0; n < g.num_nodes(); n++)
      degree_id_pairs[n] = std::make_pair(g.out_degree(n), n);
    std::sort(degree_id_pairs.begin(), degree_id_pairs.end(),
              std::greater<degree_node_p>());
    pvector<NodeID_> degrees(g.num_nodes());
    pvector<NodeID_> new_ids(g.num_nodes());
    #pragma omp parallel for
    for (NodeID_ n=0; n < g.num_nodes(); n++) {
      degrees[n] = degree_id_pairs[n].first;
      new_ids[degree_id_pairs[n].second] = n;
    }
    pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
    DestID_* neighs = new DestID_[offsets[g.num_nodes()]];
    DestID_** index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs);
    #pragma omp parallel for
    for (NodeID_ u=0; u < g.num_nodes(); u++) {
      for (NodeID_ v : g.out_neigh(u))
        neighs[offsets[new_ids[u]]++] = new_ids[v];
      std::sort(index[new_ids[u]], index[new_ids[u]+1]);
    }
    t.Stop();
    PrintTime("Relabel", t.Seconds());
    return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), index, neighs);
  }
};

#endif  // BUILDER_H_
