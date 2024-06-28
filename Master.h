#include "../../snap-core/Snap.h"
#include "cliques.h"
#include <vector>
#include <algorithm>
#include <queue>
#include <iostream>
using namespace std;

typedef TVec<TInt> TIntV;
typedef TVec<PUNGraph> TUNGraV;
typedef TPt<TUNGraph> PUNGraph;
typedef TVec<TIntV> TIntVIntV;
typedef TVec <TPair<TIntV, TIntV>> TPrVIntV;

class Master {
    public:
    	int nfs;// the number of all followers;
	    int acost; // the number of inserted edges;
        int k; 
        int b;
        PUNGraph G;
        Master(PUNGraph G, int k, int b);
        void Anchoring();
        vector<double> GroupSelection(TIntV& G_S, TIntV& delta_S, TIntV& delta_S_bar);
        vector<double> GroupSelection_test(TIntV& G_S, TIntV& delta_S, TIntV& delta_S_bar, TIntV& Expanded_Vertex);
        TIntV GetBoundary(TIntV G_S, TIntV &delta_S_bar);
        void Load_kvcc(TIntVIntV& kvcc_array);
        void update_neighbour(TIntVIntV& S, TIntIntVH& in_neighs, TIntIntVH& out_neighs, int v, TIntV& res);
};