#include <stdio.h>
#include <iostream>
#include "Master.h"
#include "getopt.h"
using namespace std;

int main(int argc, char* argv[]) {
  int k, b;
  double h;
  // string alg, seed, mergeMethod, expandMethod;
  string alg;
  TStr dataset;
  int o;
  string vcc_data;
  bool dis_cck = false;
  bool dis_cmv = false;
  const char* optstring = "a:d:b:k:c:s:m:e:t:h:u:";
  while ((o = getopt(argc, argv, optstring)) != -1) {
    switch (o) {
      case 'a':
        alg = optarg;
        printf("opt is a, oprarg is: %s\n", optarg);
        break;
      case 'd':
        dataset = optarg;
        printf("opt is d, oprarg is: %s\n", optarg);
        break;
      case 'c':
        vcc_data = optarg;
        printf("opt is c, oprarg is %s\n", optarg);
      case 'b':
        b = atoi(optarg);
        printf("opt is b, oprarg is: %s\n", optarg);
        break;
      case 'k':
        k = atoi(optarg);
        printf("opt is k, oprarg is: %s\n", optarg);
        break;
      case 'h':
        h = atoi(optarg);
        printf("opt is h, oprarg is: %s\n", optarg);
        break;
      case 'u': {
        std::string disable_opt = optarg;
        if (disable_opt == "cck") {
          dis_cck = true;
        } else if (disable_opt == "cmv") {
          dis_cmv = true;
        } else if (disable_opt == "all") {
          dis_cck = true;
          dis_cmv = true;
        } else {
          std::cerr << "错误：无效的 -u 选项 " << disable_opt << std::endl;
          std::cerr << "用法：-u [cck|cmv|all]" << std::endl;
          return 1;
        }
        break;
      }
      case '?':
        printf("Correct Usage:\n");
        return 0;
    }
  }

  PUNGraph G = TSnap::LoadEdgeList<PUNGraph>(dataset + ".txt", 0, 1);

  printf("G: \nnode_nums = %d, edge_nums = %d\n", G->GetNodes(), G->GetEdges());

  Master master(G, k, b, dis_cck, dis_cmv);
  master.Anchoring(alg, vcc_data, h);
  // master.Exact_Anchoring(alg, vcc_data);
}
