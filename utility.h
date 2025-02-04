#ifndef UTILITY_H
#define UTILITY_H

#include <vector>
#include "../../snap-core/Snap.h"

typedef TPt<TUNGraph> PUNGraph;
typedef TIntV TIntVector;
PUNGraph RemoveInternalEdges(const PUNGraph& G, const std::vector<int>& S_L, const std::vector<int>& S_R);

#endif // UTILITY_H