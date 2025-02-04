#ifndef MAXMATCH_H
#define MAXMATCH_H

#include <vector>
#include <utility>
// 假设这是 snap 库的头文件路径
#include "../../snap-core/Snap.h"

typedef TPt<TUNGraph> PUNGraph;

// 声明 Hopcroft-Karp 相关函数
int HopcroftKarp(const PUNGraph& Graph, const std::vector<int>& LeftVertices, const std::vector<int>& RightVertices, std::vector<std::pair<int, int>>& matches);

#endif // MAXMATCH_H