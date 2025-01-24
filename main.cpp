// #include <stdio.h>
// #include <iostream>
// #include "Master.h"
// #include "getopt.h"
// using namespace std;

// int main(int argc, char* argv[]) {
//   int k, b;
//   // string alg, seed, mergeMethod, expandMethod;
//   string alg;
//   TStr dataset;
//   int o;
//   string vcc_data;
//   const char* optstring = "a:d:b:k:c:s:m:e:t:";
//   while ((o = getopt(argc, argv, optstring)) != -1) {
//     switch (o) {
//       case 'a':
//         alg = optarg;
//         printf("opt is a, oprarg is: %s\n", optarg);
//         break;
//       case 'd':
//         dataset = optarg;
//         printf("opt is d, oprarg is: %s\n", optarg);
//         break;
//       case 'c':
//         vcc_data = optarg;
//         printf("opt is c, oprarg is %s\n", optarg);
//       case 'b':
//         b = atoi(optarg);
//         printf("opt is b, oprarg is: %s\n", optarg);
//         break;
//       case 'k':
//         k = atoi(optarg);
//         printf("opt is k, oprarg is: %s\n", optarg);
//         break;
//       case '?':
//         printf("Correct Usage:\n");
//         return 0;
//     }
//   }

//   PUNGraph G = TSnap::LoadEdgeList<PUNGraph>(dataset + ".txt", 0, 1);

//   printf("G: \nnode_nums = %d, edge_nums = %d\n", G->GetNodes(), G->GetEdges());

//   Master master(G, k, b);
//   master.Anchoring(alg, vcc_data);
// }

#include <iostream>
#include <vector>
#include "maxMatch.h"
#include "Snap.h"

int main() {
    // 创建一个无向图
    PUNGraph Graph = TUNGraph::New();

    // 添加节点
    for (int i = 0; i < 8; i++) {
        Graph->AddNode(i);
    }

    // 定义左部顶点
    vector<int> LeftVertices = {0, 1, 2, 3};

    // 添加边
    Graph->AddEdge(0, 4);
    Graph->AddEdge(0, 5);
    Graph->AddEdge(1, 5);
    Graph->AddEdge(1, 6);
    Graph->AddEdge(2, 6);
    Graph->AddEdge(2, 7);
    Graph->AddEdge(3, 7);
    Graph->AddEdge(3, 4);

    // 计算最大匹配
    int MaxMatch = HopcroftKarp(Graph, LeftVertices);

    cout << "最大匹配数: " << MaxMatch << endl;

    return 0;
}