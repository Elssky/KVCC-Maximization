#include <stdio.h>
#include <iostream>
#include "Master.h"
#include "getopt.h"
using namespace std;

int main(int argc, char* argv[]) {
  int k, b;
  // string alg, seed, mergeMethod, expandMethod;
  string alg;
  TStr dataset;
  int o;
  string vcc_data;
  const char* optstring = "a:d:b:k:c:s:m:e:t:";
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
      case '?':
        printf("Correct Usage:\n");
        return 0;
    }
  }

  PUNGraph G = TSnap::LoadEdgeList<PUNGraph>(dataset + ".txt", 0, 1);

  printf("G: \nnode_nums = %d, edge_nums = %d\n", G->GetNodes(), G->GetEdges());

  Master master(G, k, b);
  master.Anchoring(alg, vcc_data);
}

// #include "maxMatch.h"
// #include <iostream>
// int main() {
//     // 创建一个无向图
//     PUNGraph Graph = TUNGraph::New();

//     // 添加节点
//     for (int i = 0; i < 6; i++) {
//         Graph->AddNode(i);
//     }

//     // 定义左部顶点
//     vector<int> LeftVertices = {0, 1, 2};

//     // 添加边
//     Graph->AddEdge(0, 3);
//     Graph->AddEdge(0, 4);
//     Graph->AddEdge(1, 4);
//     Graph->AddEdge(1, 5);
//     Graph->AddEdge(2, 4);

//     // 存储匹配的顶点对
//     std::vector<std::pair<int, int>> matches;

//     // 计算最大匹配
//     int MaxMatch = HopcroftKarp(Graph, LeftVertices, matches);

//     cout << "最大匹配数: " << MaxMatch << endl;
//     cout << "匹配的顶点对:" << endl;
//     for (const auto& match : matches) {
//         cout << "(" << match.first << ", " << match.second << ")" << endl;
//     }

//     return 0;
// }

// #include "utility.h"
// #include <iostream>

// int main() {
//     // 创建一个大图 G
//     PUNGraph G = TUNGraph::New();
//     for (int i = 0; i < 10; ++i) {
//         G->AddNode(i);
//     }
//     G->AddEdge(0, 1); // S_L 内部边
//     G->AddEdge(0, 5); // 跨部边
//     G->AddEdge(2, 6); // 跨部边
//     G->AddEdge(3, 4); // S_L 内部边
//     G->AddEdge(5, 6); // S_R 内部边
//     G->AddEdge(7, 8); // S_R 内部边

//     // 定义 S_L 和 S_R
//     TIntV S_L;
//     S_L.Add(0);
//     S_L.Add(1);
//     S_L.Add(2);
//     S_L.Add(3);
//     S_L.Add(4);

//     TIntV S_R;
//     S_R.Add(5);
//     S_R.Add(6);
//     S_R.Add(7);
//     S_R.Add(8);
//     S_R.Add(9);

//     // 移除内部边
//     PUNGraph NewSubGraph = RemoveInternalEdges(G, S_L, S_R);

//     // 输出新子图的边
//     std::cout << "新子图的边:" << std::endl;
//     for (TUNGraph::TEdgeI EI = NewSubGraph->BegEI(); EI < NewSubGraph->EndEI(); EI++) {
//         std::cout << "(" << EI.GetSrcNId() << ", " << EI.GetDstNId() << ")" << std::endl;
//     }

//     return 0;
// }