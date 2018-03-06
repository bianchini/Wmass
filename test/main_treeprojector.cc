#include "../interface/TreeProjector.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[]){
  TreeProjector* tp = new TreeProjector( "/scratch/bianchini/tree_test.root", "tree", "./test/test_tp.root", atoi(argv[1]));
  tp->run_pt_bias_vs_qt();
  //tp->run_test();
  return 0;
}
