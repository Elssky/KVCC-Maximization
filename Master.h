#ifndef MASTER_H
#include <algorithm>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include "../../snap-core/Snap.h"
#include "cliques.h"
#include "maxMatch.h"
#include "utility.h"
#include <chrono>
using namespace std;

typedef TVec<TInt> TIntV;
typedef TVec<PUNGraph> TUNGraV;
typedef TPt<TUNGraph> PUNGraph;
typedef TVec<TIntV> TIntVIntV;
typedef TVec<TPair<TIntV, TIntV>> TPrVIntV;
typedef THash<TIntV, TIntV> TIntVTIntVH;

// Define a hash function for std::pair<int, int>
struct pair_hash {
  template <class T1, class T2>
  std::size_t operator()(const std::pair<T1, T2>& pair) const {
    return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
  }
};

struct Compare {
  bool operator()(const std::pair<std::pair<int, int>, double>& a,
                  const std::pair<std::pair<int, int>, double>& b) {
    return a.second < b.second;  // 按照 double 值从大到小排序
  }
};

struct Compare_new {
  //按照 double 值(r)从大到小排序，当 double 值相同时再按照第三个 int
  //值(cost)从小到大排序
  bool operator()(const std::tuple<int, int, int, double>& a,
                  const std::tuple<int, int, int, double>& b) const {
    // 先比较 double 值
    if (std::get<3>(a) != std::get<3>(b)) {
      return std::get<3>(a) < std::get<3>(b);
    }
    // 如果 double 值相同，比较第三个 int 值
    return std::get<2>(a) > std::get<2>(b);
  }
};

// use priority queue to select the clique with max size + level (|c| + i)
auto comp = [](const std::pair<std::vector<int>, std::pair<int, int>>& a,
               const std::pair<std::vector<int>, std::pair<int, int>>& b) {
  const auto& a_first = a.second.first;
  const auto& b_first = b.second.first;
  const auto& a_second = a.second.second;
  const auto& b_second = b.second.second;

  return (a_first == b_first) ? (a_second > b_second) : (a_first < b_first);
};

class Master {
 public:
  int nfs;    // the number of all followers;
  int acost;  // the number of inserted edges;
  int k;
  int b;
  PUNGraph G;
  Master(PUNGraph G, int k, int b);
  void Anchoring(string alg, string vcc_data);
  void GroupSelection_together(
      TIntV& G_S, TIntV& delta_S, TIntV& delta_S_bar,
      unordered_set<pair<int, int>, pair_hash>& Inserted_Edge,
      TIntV& Expanded_Vertex);
  void GroupSelection_multi_vertex(
      TIntV& G_S, TIntV& delta_S, TIntV& delta_S_bar,
      unordered_set<pair<int, int>, pair_hash>& Inserted_Edge,
      TIntV& Expanded_Vertex);
  void GroupSelection_single_vertex(
      TIntV& G_S, TIntV& delta_S, TIntV& delta_S_bar,
      unordered_set<pair<int, int>, pair_hash>& Inserted_Edge,
      TIntV& Expanded_Vertex);
  void Merge_overlap_vcc(
      TIntVIntV& VCCs, unordered_set<pair<int, int>, pair_hash>& Inserted_Edge,
      TIntV& Expanded_Vertex);
  void Merge_adjacent_vcc(
      TIntVIntV& VCCs, unordered_set<pair<int, int>, pair_hash>& Inserted_Edge,
      TIntV& Expanded_Vertex);
  TIntV GetBoundary(TIntV G_S, TIntV& delta_S_bar);
  void Load_kvcc(TIntVIntV& kvcc_array, string vcc_data);
  void update_neighbour(TIntVIntV& S, TIntIntVH& in_neighs,
                        TIntIntVH& out_neighs, int v, TIntV& res, int& level);
  void sort_by_deg(TIntV& Vcc);
  vector<int> get_common_set(TIntV& Vcc_1, TIntV& Vcc_2);
  TIntV get_difference_set(TIntV& Vcc_1, TIntV& Vcc_2);
  int get_degree_in_set(int vertex, const TIntV& set);
  TIntV getNeighborhood(const TIntV& set, const TIntV& targetSet);
  TIntV getIntersection(const TIntV& set1, const TIntV& set2);
  std::pair<int, int> find_max_in_com_neigh(
      const std::unordered_map<std::pair<int, int>, std::vector<int>,
                               pair_hash>& com_neigh);
  std::pair<int, int> find_max_in_gamma(
      const std::unordered_map<std::pair<int, int>, int, pair_hash>& gamma);

  void CalConnectKVcc(TIntVIntV& VCCs, int& i_star, int& j_star,
                      vector<vector<int>>& t_star, double& r_star,
                      int& cost_star, double& best_gain);

  void CalMulVerices(TIntVIntV& VCCs, int& i_star, vector<int>& mc_star,
                     double& r_star, int& cost_star, double& best_gain);

  void ExpSinVertices(TIntVIntV& VCCs, double& r,
                      unordered_set<pair<int, int>, pair_hash>& Inserted_Edge, int &v,  int& vcc_idx, double& best_gain);

  void MerConnectKVcc(TIntV& VCC_i, TIntV& VCC_j, vector<int> t_ij,
                      vector<int> t_ji, double r_ij, int cost,
                      unordered_set<pair<int, int>, pair_hash>& Inserted_Edge);

  void ExpMulVertices(TIntV& VCC_i, TIntV& mc_j, double r_ij,
                      unordered_set<pair<int, int>, pair_hash>& Inserted_Edge);
  void InsertEdgesIntoGraph(
      const std::unordered_set<std::pair<int, int>, pair_hash>& edgeSet) {
    for (const auto& edge : edgeSet) {
      int src = edge.first;
      int dst = edge.second;
      if (!G->IsNode(src)) {
        G->AddNode(src);
      }
      if (!G->IsNode(dst)) {
        G->AddNode(dst);
      }
      if (!G->IsEdge(src, dst)) {
        G->AddEdge(src, dst);
      }
    }
  }
};
#endif