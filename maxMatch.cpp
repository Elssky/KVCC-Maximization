#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <unordered_map>
#include "maxMatch.h"

using namespace std;

const int INF = std::numeric_limits<int>::max();

// BFS 函数
bool BFS(const PUNGraph& Graph, const std::vector<int>& LeftVertices, const std::vector<int>& RightVertices, std::vector<int>& PairU, std::vector<int>& PairV, std::vector<int>& Dist, const std::unordered_map<int, int>& leftVertexToIndex, const std::unordered_map<int, int>& rightVertexToIndex) {
    std::queue<int> Q;
    for (int u : LeftVertices) {
        int uIndex = leftVertexToIndex.at(u);
        if (PairU[uIndex] == -1) {
            Dist[uIndex] = 0;
            Q.push(u);
        } else {
            Dist[uIndex] = INF;
        }
    }
    bool foundAugmentingPath = false;

    while (!Q.empty()) {
        int u = Q.front();
        Q.pop();
        int uIndex = leftVertexToIndex.at(u);
        if (Dist[uIndex] < INF) {
            TUNGraph::TNodeI NI = Graph->GetNI(u);
            for (int i = 0; i < NI.GetOutDeg(); ++i) {
                int v = NI.GetNbrNId(i);
                if (rightVertexToIndex.count(v) == 0) continue; // 跳过不在右部顶点集合中的顶点
                int vIndex = rightVertexToIndex.at(v);
                if (PairV[vIndex] == -1) {
                    foundAugmentingPath = true;
                } else if (leftVertexToIndex.count(PairV[vIndex]) > 0) {
                    int pairVIndex = leftVertexToIndex.at(PairV[vIndex]);
                    if (Dist[pairVIndex] == INF) {
                        Dist[pairVIndex] = Dist[uIndex] + 1;
                        Q.push(PairV[vIndex]);
                    }
                }
            }
        }
    }
    return foundAugmentingPath;
}

// DFS 函数
bool DFS(const PUNGraph& Graph, int u, std::vector<int>& PairU, std::vector<int>& PairV, std::vector<int>& Dist, const std::unordered_map<int, int>& leftVertexToIndex, const std::unordered_map<int, int>& rightVertexToIndex) {
    if (u != -1) {
        int uIndex = leftVertexToIndex.at(u);
        TUNGraph::TNodeI NI = Graph->GetNI(u);
        for (int i = 0; i < NI.GetOutDeg(); ++i) {
            int v = NI.GetNbrNId(i);
            if (rightVertexToIndex.count(v) == 0) continue; // 跳过不在右部顶点集合中的顶点
            int vIndex = rightVertexToIndex.at(v);
            if (PairV[vIndex] == -1 || 
                (leftVertexToIndex.count(PairV[vIndex]) > 0 && Dist[leftVertexToIndex.at(PairV[vIndex])] == Dist[uIndex] + 1 && 
                 DFS(Graph, PairV[vIndex], PairU, PairV, Dist, leftVertexToIndex, rightVertexToIndex))) {
                PairV[vIndex] = u;
                PairU[uIndex] = v;
                return true;
            }
        }
        Dist[uIndex] = INF;
        return false;
    }
    return true;
}

// Hopcroft-Karp 算法实现
int HopcroftKarp(const PUNGraph& Graph, const std::vector<int>& LeftVertices, const std::vector<int>& RightVertices, std::vector<std::pair<int, int>>& matches) {
    int MaxMatch = 0;
    int N = Graph->GetNodes();

    // 建立左部顶点编号到连续索引的映射
    std::unordered_map<int, int> leftVertexToIndex;
    int leftIndex = 0;
    for (int vertex : LeftVertices) {
        leftVertexToIndex[vertex] = leftIndex++;
    }

    // 建立右部顶点编号到连续索引的映射
    std::unordered_map<int, int> rightVertexToIndex;
    int rightIndex = 0;
    for (int vertex : RightVertices) {
        rightVertexToIndex[vertex] = rightIndex++;
    }

    std::vector<int> PairU(LeftVertices.size(), -1);  // 左部顶点的匹配
    std::vector<int> PairV(RightVertices.size(), -1);  // 右部顶点的匹配
    std::vector<int> Dist(LeftVertices.size());

    // std::cout << "节点数: " << N << std::endl;
    int iteration = 1;
    // 主循环
    while (BFS(Graph, LeftVertices, RightVertices, PairU, PairV, Dist, leftVertexToIndex, rightVertexToIndex)) {
        // std::cout << "第 " << iteration << " 次 BFS 找到增广路径，开始 DFS 扩展匹配..." << std::endl;
        for (int u : LeftVertices) {
            int uIndex = leftVertexToIndex.at(u);
            if (PairU[uIndex] == -1 && DFS(Graph, u, PairU, PairV, Dist, leftVertexToIndex, rightVertexToIndex)) {
                MaxMatch++;
                // std::cout << "第 " << iteration << " 次迭代中，左部顶点 " << u << " 匹配成功，当前最大匹配数: " << MaxMatch << std::endl;
            }
        }
        iteration++;
    }
    // std::cout << "Hopcroft-Karp 算法结束，未找到更多增广路径" << std::endl;

    // 记录匹配的顶点对
    // std::cout << "记录最终匹配的顶点对..." << std::endl;
    for (int u : LeftVertices) {
        int uIndex = leftVertexToIndex.at(u);
        if (PairU[uIndex] != -1) {
            int v = PairU[uIndex];
            if (rightVertexToIndex.count(v) > 0) {
                matches.emplace_back(u, v);
                // std::cout << "匹配对: (" << u << ", " << v << ")" << std::endl;
            }
        }
    }

    return MaxMatch;
}