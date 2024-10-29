#include <algorithm>
#include <iostream>
#include <queue>
#include <string>
#include <unordered_set>
#include <vector>
#include "../../snap-core/Snap.h"
#include "cliques.h"
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

class Master {
public:
   int nfs;    // the number of all followers;
   int acost;  // the number of inserted edges;
   int k;
   int b;
   PUNGraph G;
   Master(PUNGraph G, int k, int b);
   void Anchoring(string alg, string vcc_data);
   vector<double> GroupSelection_together(TIntV& G_S, TIntV& delta_S,
      TIntV& delta_S_bar,
      unordered_set<pair<int, int>, pair_hash>& Inserted_Edge,
      TIntV& Expanded_Vertex);
   vector<double> GroupSelection_multi_vertex(TIntV& G_S, TIntV& delta_S,
      TIntV& delta_S_bar,
      unordered_set<pair<int, int>, pair_hash>& Inserted_Edge,
      TIntV& Expanded_Vertex);
   vector<double> GroupSelection_single_vertex(TIntV& G_S, TIntV& delta_S,
      TIntV& delta_S_bar,
      unordered_set<pair<int, int>, pair_hash>& Inserted_Edge,
      TIntV& Expanded_Vertex);
   TIntV GetBoundary(TIntV G_S, TIntV& delta_S_bar);
   void Load_kvcc(TIntVIntV& kvcc_array, string vcc_data);
   void update_neighbour(TIntVIntV& S, TIntIntVH& in_neighs,
      TIntIntVH& out_neighs, int v, TIntV& res, int& level);
};