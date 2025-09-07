#include "utility.h"

PUNGraph RemoveInternalEdges(const PUNGraph& G, const std::vector<int>& S_L, const std::vector<int>& S_R) {
    // 创建子图节点集合
    TIntV SubGraphNodes;
    for (int node : S_L) {
        SubGraphNodes.Add(node);
    }
    for (int node : S_R) {
        SubGraphNodes.Add(node);
    }

    // 创建子图
    PUNGraph SubGraph = TSnap::GetSubGraph(G, SubGraphNodes);

    // 遍历子图的边
    TIntPrV EdgesToRemove;
    for (TUNGraph::TEdgeI EI = SubGraph->BegEI(); EI < SubGraph->EndEI(); EI++) {
        int src = EI.GetSrcNId();
        int dst = EI.GetDstNId();
        bool bothInSL = false;
        bool bothInSR = false;

        // 检查两个端点是否都在 S_L 中
        for (int node : S_L) {
            if (src == node) {
                for (int otherNode : S_L) {
                    if (dst == otherNode) {
                        bothInSL = true;
                        break;
                    }
                }
            }
            if (bothInSL) break;
        }

        // 检查两个端点是否都在 S_R 中
        for (int node : S_R) {
            if (src == node) {
                for (int otherNode : S_R) {
                    if (dst == otherNode) {
                        bothInSR = true;
                        break;
                    }
                }
            }
            if (bothInSR) break;
        }

        // 如果两个端点都在 S_L 或者都在 S_R 中，标记该边为要移除的边
        if (bothInSL || bothInSR) {
            EdgesToRemove.Add(TIntPr(src, dst));
        }
    }

    // 移除标记的边
    for (int i = 0; i < EdgesToRemove.Len(); ++i) {
        SubGraph->DelEdge(EdgesToRemove[i].Val1, EdgesToRemove[i].Val2);
    }

    return SubGraph;
}