#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include "Snap.h"

using namespace std;

const int INF = std::numeric_limits<int>::max();

// BFS 函数
bool BFS(const PUNGraph& Graph, const std::vector<int>& LeftVertices, std::vector<int>& PairU, std::vector<int>& PairV, std::vector<int>& Dist) {
    std::queue<int> Q;
    for (int u : LeftVertices) {
        if (PairU[u] == -1) {
            Dist[u] = 0;
            Q.push(u);
        } else {
            Dist[u] = INF;
        }
    }
    bool foundAugmentingPath = false;

    while (!Q.empty()) {
        int u = Q.front();
        Q.pop();
        if (Dist[u] < INF) {
            TUNGraph::TNodeI NI = Graph->GetNI(u);
            for (int i = 0; i < NI.GetOutDeg(); ++i) {
                int v = NI.GetNbrNId(i);
                if (PairV[v] == -1) {
                    foundAugmentingPath = true;
                } else if (Dist[PairV[v]] == INF) {
                    Dist[PairV[v]] = Dist[u] + 1;
                    Q.push(PairV[v]);
                }
            }
        }
    }
    return foundAugmentingPath;
}

// DFS 函数
bool DFS(const PUNGraph& Graph, int u, std::vector<int>& PairU, std::vector<int>& PairV, std::vector<int>& Dist) {
    if (u != -1) {
        TUNGraph::TNodeI NI = Graph->GetNI(u);
        for (int i = 0; i < NI.GetOutDeg(); ++i) {
            int v = NI.GetNbrNId(i);
            if (PairV[v] == -1 || (Dist[PairV[v]] == Dist[u] + 1 && DFS(Graph, PairV[v], PairU, PairV, Dist))) {
                PairV[v] = u;
                PairU[u] = v;
                return true;
            }
        }
        Dist[u] = INF;
        return false;
    }
    return true;
}

// Hopcroft-Karp 算法实现
int HopcroftKarp(const PUNGraph& Graph, const std::vector<int>& LeftVertices) {
    int MaxMatch = 0;
    int N = Graph->GetNodes();
    std::vector<int> PairU(N, -1);  // 左部顶点的匹配
    std::vector<int> PairV(N, -1);  // 右部顶点的匹配
    std::vector<int> Dist(N);

    // 主循环
    while (BFS(Graph, LeftVertices, PairU, PairV, Dist)) {
        for (int u : LeftVertices) {
            if (PairU[u] == -1 && DFS(Graph, u, PairU, PairV, Dist)) {
                MaxMatch++;
            }
        }
    }

    return MaxMatch;
}