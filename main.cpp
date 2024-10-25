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
    const char* optstring = "a:d:b:k:s:m:e:t:";
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

    PUNGraph G =
        TSnap::LoadEdgeList<PUNGraph>("./dataset/" + dataset + ".txt", 0, 1);

    printf("G: \nnode_nums = %d, edge_nums = %d\n", G->GetNodes(),
        G->GetEdges());

    Master master(G, k, b);
    master.Anchoring(alg);
}