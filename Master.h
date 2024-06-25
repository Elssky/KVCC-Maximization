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
        PUNGraph G;
        Master(PUNGraph G, int k);
        void Anchoring(int b);
        vector<double> GroupSelection(int b, TIntV& G_S, TIntV& delta_S, TIntV& delta_S_bar);
        TIntV GetBoundary(TIntV G_S, TIntV &delta_S_bar);
};