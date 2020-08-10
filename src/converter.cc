// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <iostream>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "reader.h"
#include "writer.h"

using namespace std;

int main(int argc, char* argv[]) {
  CLConvert cli(argc, argv, "converter");
  cli.ParseArgs();
  
  /*int* arr = new (int[67108864]);
  std::cout << "size of arr: " << sizeof(*arr)  << " addresss: " << arr << std::endl;
  for(int i = 0; i < (67108864); i++){
    arr[i] = 0;
  }
  std::cout << "size of arr: " << sizeof(*arr) << " addresss: " << arr << std::endl;
  std::realloc(arr, 67108864/2);
  if(arr == nullptr){
    std::cout << "its null! realloc failed\n";
  }
  for(int i = 0; i < (67108864/2); i++){
    arr[i] = 0;
  }
  std::cout << "size of arr: " << sizeof(*arr) << " addresss: " << arr << std::endl;
  delete arr;
  exit(0);*/

  if (cli.out_weighted()) {
    WeightedBuilder bw(cli);
    WGraph wg = bw.MakeGraph(cli.lowmem());
    wg.PrintStats();
    WeightedWriter ww(wg);
    ww.WriteGraph(cli.out_filename(), cli.out_sg());
  } else {
    Builder b(cli);
    std::cout << "calling makegraph\n";
    Graph g = b.MakeGraph(cli.lowmem());
    std::cout << "done w squishgraph and GRAPH IS MADE\n";
    g.PrintStats();
    Writer w(g);
    w.WriteGraph(cli.out_filename(), cli.out_sg());
    std::cout << "DONE WITH CONVERTER\n";
  }
  return 0;
}
