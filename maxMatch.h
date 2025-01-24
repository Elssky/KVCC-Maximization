#ifndef MAXMATCH_H
#define MAXMATCH_H

#include <vector>
#include "Snap.h"
using namespace std;

typedef TPt<TUNGraph> PUNGraph;

// 声明 Hopcroft-Karp 相关函数
int HopcroftKarp(const PUNGraph& Graph, const vector<int>& LeftVertices);

#endif // MAXMATCH_H