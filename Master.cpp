#include "Master.h"

#define _DEBUG
#ifdef _DEBUG
#define debug_print(msg) std::cout << msg
#else
#define debug_print(msg)
#endif

Master::Master(PUNGraph G, int k, int b) {
  this->G = G;
  acost = 0;
  nfs = 0;
  this->b = b;
  this->k = k;
}

void Master::Anchoring(string alg, string vcc_data) {
  double t_begin = (double)clock();
  int round = 0;
  double group_anchor_time = 0.0, vertex_anchor_time = 0.0;
  TIntVIntV kvcc_array;
  TIntV kvcc, delta_S, delta_S_bar;
  TIntV Expanded_Vertex;
  Load_kvcc(kvcc_array, vcc_data);
  // select kvcc to expand
  // kvcc = kvcc_array[1];
  // unordered_set<pair<int, int>, pair_hash> Inserted_Edge;
  vector<int> Io, Ie, J;
  vector<double> Ro, Re;
  vector<vector<int>> T, MC;
  double r;
  unordered_set<pair<int, int>, pair_hash> Inserted_Edge;

  // CalConnectKVcc(kvcc_array, Io, J, T, Ro);
  CalMulVerices(kvcc_array, Ie, MC, Re);

  // double ra = *Ro.begin();
  double rb = *Re.begin();
  ExpSinVertices(kvcc_array, r, Inserted_Edge);

  // while (acost < b) {
  //   if (!Ro.empty() && !Re.empty()) {
  //     int ra = *Ro.begin();
  //     int rb = *Re.begin();
  //     ExpSinVertices(kvcc_array, r, Inserted_Edge);
  //     // std::cout << "ra: " << ra << ", rb: " << rb << std::endl;
  //   } else {
  //     break;
  //   }
  // }

  // while (acost < b) {
  //   cout << " -- Anchoring round: " << round++ << endl;
  //   double vertex_begin = (double)clock();
  //   double node_score = 0;

  //   // Compute by Multiple Vertex Anchoring
  //   // vector<double> group;

  // }

  // cout << "acost: " << acost << endl;
  // cout << "gain: " << Expanded_Vertex.Len() << endl;
  // cout << "Expanded_Vertex:";
  // for (TIntV::TIter NI = Expanded_Vertex.BegI(); NI < Expanded_Vertex.EndI();
  //      NI++) {
  //   cout << *NI << " ";
  // }
  // cout << endl;
  // double t_end = (double)clock();
  // cout << "the anchoring time is:" << (t_end - t_begin) / CLOCKS_PER_SEC <<
  // "s."
  //      << endl;
}

// void Master::Anchoring(string alg, string vcc_data) {
//   double t_begin = (double)clock();
//   int round = 0;
//   double group_anchor_time = 0.0, vertex_anchor_time = 0.0;
//   TIntVIntV kvcc_array;
//   TIntV kvcc, delta_S, delta_S_bar;
//   TIntV Expanded_Vertex;
//   Load_kvcc(kvcc_array, vcc_data);
//   // select kvcc to expand
//   kvcc = kvcc_array[1];
//   unordered_set<pair<int, int>, pair_hash> Inserted_Edge;
//   while (acost < b) {
//     cout << " -- Anchoring round: " << round++ << endl;
//     double vertex_begin = (double)clock();
//     double node_score = 0;

//     // Compute by Multiple Vertex Anchoring
//     // vector<double> group;

//     int insEdge_size = Inserted_Edge.size();
//     if (alg == string("t")) {
//       GroupSelection_together(kvcc, delta_S, delta_S_bar, Inserted_Edge,
//                               Expanded_Vertex);
//     } else if (alg == string("m")) {
//       GroupSelection_multi_vertex(kvcc, delta_S, delta_S_bar, Inserted_Edge,
//                                   Expanded_Vertex);
//     } else if (alg == string("s")) {
//       GroupSelection_single_vertex(kvcc, delta_S, delta_S_bar, Inserted_Edge,
//                                    Expanded_Vertex);
//     } else if (alg == string("mo")) {
//       Merge_overlap_vcc(kvcc_array, Inserted_Edge, Expanded_Vertex);
//     } else if (alg == string("ma")) {
//       Merge_adjacent_vcc(kvcc_array, Inserted_Edge, Expanded_Vertex);
//     } else {
//       cout << "wrong alg parameter" << endl;
//       return;
//     }
//     // if no edge can be inserted, stop
//     if (insEdge_size == Inserted_Edge.size()) {
//       break;
//     }
//   }
//   cout << "acost: " << acost << endl;
//   cout << "gain: " << Expanded_Vertex.Len() << endl;
//   cout << "Expanded_Vertex:";
//   for (TIntV::TIter NI = Expanded_Vertex.BegI(); NI < Expanded_Vertex.EndI();
//        NI++) {
//     cout << *NI << " ";
//   }
//   cout << endl;
//   double t_end = (double)clock();
//   cout << "the anchoring time is:" << (t_end - t_begin) / CLOCKS_PER_SEC <<
//   "s."
//        << endl;
// }

/* together 和 multi的区别在于同时考虑所有不同S_i层级的clique，
而multi是i从大到小考虑
*/
void Master::GroupSelection_together(
    TIntV& G_S, TIntV& delta_S, TIntV& delta_S_bar,
    unordered_set<pair<int, int>, pair_hash>& Inserted_Edge,
    TIntV& Expanded_Vertex) {
  cout << "kvcc: ";
  for (TIntV::TIter NI = G_S.BegI(); NI < G_S.EndI(); NI++) {
    cout << *NI << " ";
  }
  cout << endl;

  // calcuate one-hop neighbor of Graph G as candidate set
  delta_S = GetBoundary(G_S, delta_S_bar);
  TIntV G_sub = delta_S;
  G_sub.AddVMerged(delta_S_bar);
  PUNGraph G_Cand = TSnap::GetSubGraph(G, G_sub);

  // calculate the edge num from vertex in delta_S_bar to k-vcc,
  // divide it into S_(k), S_(k-1), S_(k-2), ...
  TIntVIntV S;
  for (int i = 0; i <= k; i++) {
    S.Add({});
  }
  TIntV res = {};
  TIntV nb_u1 = {}, nb_u2 = {};

  // in_neighs stores the neighbors of each candidate node belonging to k-vcc
  // out_neighs stores the neighbors of each candidate node belonging to
  // delta_S_bar
  TIntIntVH in_neighs, out_neighs;
  for (TIntV::TIter TI = delta_S_bar.BegI(); TI < delta_S_bar.EndI(); TI++) {
    nb_u1.Clr();
    nb_u2.Clr();
    TUNGraph::TNodeI v = G_Cand->GetNI(*TI);
    for (int i = 0; i < v.GetInDeg(); ++i) {
      if (delta_S.IsIn(v.GetInNId(i)))  // only consider neighbors in
        // k-vcc
        nb_u1.AddMerged(v.GetInNId(i));
      else
        nb_u2.AddMerged(v.GetInNId(i));
    }
    in_neighs.AddDat(*TI, nb_u1);
    out_neighs.AddDat(*TI, nb_u2);
    int index = nb_u1.Len();

    if (index >= k) index = k;
    S[index].AddMerged(*TI);
  }

  // calculate maximal cliques MC in each S
  int num = 0;
  TIntV S_total;
  TIntVIntV cliques;
  TIntVTIntVH cliques_neighs_union;

  priority_queue<pair<vector<int>, pair<int, int>>,
                 vector<pair<vector<int>, pair<int, int>>>, decltype(comp)>
      clique_in_neighs(comp);

  int level = 0;
  for (int i = k - 1; i > 0; i--) {
    debug_print("S[" << i << "]: ");
    for (int j = 0; j < S[i].Len(); j++) {
      int v = S[i][j];
      debug_print(v << " ");
    }
    debug_print(endl);
    TIntV S_i = S[i];
    TIntV neigh_Union;
    TIntV S_i_temp;

    S_total.Clr();
    S_total.AddVMerged(S[i]);
    // to be considered how to deal with?
    // if (S_total.Len() < k - i + 1) continue;
    PUNGraph sub_G = TSnap::GetSubGraph(G_Cand, S_total);
    // PUNGraph sub_core = TSnap::GetKCore(sub_G, k - i);
    // to be considered how to deal with clique size smaller than k - i + 1
    TCliqueOverlap::GetMaxCliques(sub_G, 0, cliques);
    // for (int clique_idx = 0; clique_idx < cliques.Len(); clique_idx++) {
    for (TIntVIntV::TIter cq = cliques.BegI(); cq < cliques.EndI(); cq++) {
      neigh_Union.Clr();
      int flag = 0;
      vector<int> vec_cq;
      for (TIntV::TIter TI = cq->BegI(); TI < cq->EndI(); TI++) {
        // seems can be deleted, all vertices idx greater than 0?
        // int idx = in_neighs.GetDat(*TI).Len();
        // if (idx == 0) {
        //     flag = 1;
        //     break;
        // }
        vec_cq.push_back(*TI);
        neigh_Union.AddV(in_neighs.GetDat(*TI));
      }
      neigh_Union.Merge();
      cliques_neighs_union.AddDat(*cq, neigh_Union);
      clique_in_neighs.push({vec_cq, pair<int, int>(i + vec_cq.size(), i)});
      vec_cq.clear();
    }
  }

  while (!clique_in_neighs.empty()) {
    vector<int> c = clique_in_neighs.top().first;
    int i = clique_in_neighs.top().second.second;
    TIntV c_tmp;
    int flag = 0;
    for (auto& v : c) {
      c_tmp.Add(v);
    }
    clique_in_neighs.pop();
    TIntV c_neighs = cliques_neighs_union.GetDat(c_tmp);
    for (auto& v : c) {
      if (flag != 0) break;
      int need = max(1, k - i - static_cast<int>(c.size()) +
                            1);  // eahc v need how much edges to insert,
      // 至少需要插入1条边
      debug_print(v << " " << endl);
      for (TIntV::TIter TI = G_S.BegI(); TI < G_S.EndI(); TI++) {
        if (!in_neighs.GetDat(v).IsIn(*TI)) {
          // union neighs 数量达标, 可以随便选点
          if (c_neighs.Len() >= k) {
            if (acost >= b) {
              // TODO: 如果当前clique的budget不够，应该考虑所需budget更小的
              // cout << "acost: " << acost << endl;
              // cout << "gain: " << Expanded_Vertex.Len() << endl;
              flag = 1;
              break;
            }
            if (need == 0) {
              flag = 2;  // 满足条件，已经选完了
              break;
            }
            if (Inserted_Edge.find({*TI, v}) != Inserted_Edge.end()) continue;
            Inserted_Edge.insert({*TI, v});
            debug_print("(insert: " << *TI << " " << v << ") " << endl);
            need--;
            acost++;
          }
          // union neighs 数量不达标，必须从kvcc的其他点中选
          else if (!c_neighs.IsIn(*TI)) {
            if (acost >= b) {
              // cout << "acost: " << acost << endl;
              // cout << "gain: " << Expanded_Vertex.Len() << endl;
              flag = 1;
              break;
            }
            if (need == 0) {
              break;
            }
            // insert_edges.Add({*TI, v});
            if (Inserted_Edge.find({*TI, v}) != Inserted_Edge.end()) continue;
            Inserted_Edge.insert({*TI, v});
            debug_print("(insert: " << *TI << " " << v << ") " << endl);
            // b--;
            need--;
            acost++;
            c_neighs.Add(*TI);
          }
        }
      }
    }
    if (flag != 1 ||
        (flag == 1 && k - i - static_cast<int>(c.size()) + 1 <= 0 &&
         c_neighs.Len() >= k)) {
      for (auto& v : c) {
        // v have k neighbors in kvcc
        Expanded_Vertex.Add(v);
        debug_print("expanded: " << v << endl);
        // v is expanded, update its neighbors
        update_neighbour(S, in_neighs, out_neighs, v, Expanded_Vertex, level);
      }
    }
    if (flag == 1) {
      G_S.AddVMerged(Expanded_Vertex);
      return;
    }
    cliques.Clr();
  }
  G_S.AddVMerged(Expanded_Vertex);
  return;
  // how to select edge to insert
}

void Master::GroupSelection_multi_vertex(
    TIntV& G_S, TIntV& delta_S, TIntV& delta_S_bar,
    unordered_set<pair<int, int>, pair_hash>& Inserted_Edge,
    TIntV& Expanded_Vertex) {
  cout << "kvcc: ";
  for (TIntV::TIter NI = G_S.BegI(); NI < G_S.EndI(); NI++) {
    cout << *NI << " ";
  }
  cout << endl;
  // calcuate one-hop neighbor of Graph G as candidate set
  delta_S = GetBoundary(G_S, delta_S_bar);
  TIntV G_sub = delta_S;
  G_sub.AddVMerged(delta_S_bar);
  PUNGraph G_Cand = TSnap::GetSubGraph(G, G_sub);

  // calculate the edge num from vertex in delta_S_bar to k-vcc,
  // divide it into S_(k), S_(k-1), S_(k-2), ...
  TIntVIntV S;
  for (int i = 0; i <= k; i++) {
    S.Add({});
  }
  TIntV res = {};
  TIntV nb_u1 = {}, nb_u2 = {};

  // in_neighs stores the neighbors of each candidate node belonging to k-vcc
  // out_neighs stores the neighbors of each candidate node belonging to
  // delta_S_bar
  TIntIntVH in_neighs, out_neighs;
  for (TIntV::TIter TI = delta_S_bar.BegI(); TI < delta_S_bar.EndI(); TI++) {
    nb_u1.Clr();
    nb_u2.Clr();
    TUNGraph::TNodeI v = G_Cand->GetNI(*TI);
    for (int i = 0; i < v.GetInDeg(); ++i) {
      if (delta_S.IsIn(v.GetInNId(i)))  // only consider neighbors in
        // k-vcc
        nb_u1.AddMerged(v.GetInNId(i));
      else
        nb_u2.AddMerged(v.GetInNId(i));
    }
    in_neighs.AddDat(*TI, nb_u1);
    out_neighs.AddDat(*TI, nb_u2);
    int index = nb_u1.Len();

    // impossible situation
    // if (index >= k) index = k;
    S[index].AddMerged(*TI);
  }

  // calculate maximal cliques MC in each S
  int num = 0;
  TIntV S_total;
  TIntVIntV cliques;
  TIntVTIntVH cliques_neighs_union;

  priority_queue<pair<vector<int>, pair<int, int>>,
                 vector<pair<vector<int>, pair<int, int>>>, decltype(comp)>
      clique_in_neighs(comp);

  int level = 0;
  for (int i = k - 1; i > 0; i--) {
    debug_print("S[" << i << "]: ");
    TIntV S_i = S[i];
    TIntV neigh_Union;
    TIntV S_i_temp;

    S_total.Clr();
    S_total.AddVMerged(S[i]);
    // to be considered how to deal with?
    // if (S_total.Len() < k - i + 1) continue;
    PUNGraph sub_G = TSnap::GetSubGraph(G_Cand, S_total);
    // PUNGraph sub_core = TSnap::GetKCore(sub_G, k - i);
    // to be considered how to deal with clique size smaller than k - i + 1
    TCliqueOverlap::GetMaxCliques(sub_G, 0, cliques);
    // for (int clique_idx = 0; clique_idx < cliques.Len(); clique_idx++) {
    for (TIntVIntV::TIter cq = cliques.BegI(); cq < cliques.EndI(); cq++) {
      neigh_Union.Clr();
      int flag = 0;
      vector<int> vec_cq;
      for (TIntV::TIter TI = cq->BegI(); TI < cq->EndI(); TI++) {
        // seems can be deletd, all vertices idx greater than 0?
        // int idx = in_neighs.GetDat(*TI).Len();
        // if (idx == 0) {
        //     flag = 1;
        //     break;
        // }
        vec_cq.push_back(*TI);
        neigh_Union.AddV(in_neighs.GetDat(*TI));
      }
      neigh_Union.Merge();
      cliques_neighs_union.AddDat(*cq, neigh_Union);
      clique_in_neighs.push({vec_cq, pair<int, int>(i + vec_cq.size(), i)});
      vec_cq.clear();
    }

    // (clique_in_neighs.empty()) continue;
    if (clique_in_neighs.empty()) continue;
    vector<int> c = clique_in_neighs.top().first;
    TIntV c_tmp;
    int flag = 0;
    for (auto& v : c) {
      c_tmp.Add(v);
    }
    clique_in_neighs.pop();
    TIntV c_neighs = cliques_neighs_union.GetDat(c_tmp);
    for (auto& v : c) {
      if (flag == 1) break;
      int need = max(1, k - i - static_cast<int>(c.size()) +
                            1);  // eahc v need how much edges to insert,
      // 至少需要插入1条边
      debug_print(v << " " << endl);
      for (TIntV::TIter TI = G_S.BegI(); TI < G_S.EndI(); TI++) {
        if (!in_neighs.GetDat(v).IsIn(*TI)) {
          // union neighs 数量达标, 可以随便选点
          if (c_neighs.Len() >= k) {
            if (acost >= b) {
              // cout << "acost: " << acost << endl;
              // cout << "gain: " << Expanded_Vertex.Len() << endl;
              flag = 1;
              break;
            }
            if (need == 0) {
              break;
            }
            debug_print("(insert: " << *TI << " " << v << ") " << endl);
            need--;
            acost++;

          }
          // union neighs 数量不达标，必须从kvcc的其他点中选
          else if (!c_neighs.IsIn(*TI)) {
            if (acost >= b) {
              // cout << "acost: " << acost << endl;
              // cout << "gain: " << Expanded_Vertex.Len() << endl;
              flag = 1;
              break;
            }
            if (need == 0) {
              break;
            }
            // insert_edges.Add({*TI, v});
            if (Inserted_Edge.find({*TI, v}) != Inserted_Edge.end()) continue;
            Inserted_Edge.insert({*TI, v});
            debug_print("(insert: " << *TI << " " << v << ") " << endl);
            // b--;
            need--;
            acost++;
            c_neighs.Add(*TI);
          }
        }
      }
    }
    if (flag == 0 || flag == 1 && k - i - static_cast<int>(c.size()) + 1 <= 0 &&
                         c_neighs.Len() >= k) {
      for (auto& v : c) {
        // v have k neighbors in kvcc
        Expanded_Vertex.Add(v);
        debug_print("expanded: " << v << endl);
        // v is expanded, update its neighbors
        update_neighbour(S, in_neighs, out_neighs, v, Expanded_Vertex, level);
      }
    }
    if (flag == 1) {
      Expanded_Vertex.Merge();
      // cout << "acost: " << acost << endl;
      // cout << "gain: " << Expanded_Vertex.Len() << endl;
      G_S.AddVMerged(Expanded_Vertex);
      return;
    }
  }
  cliques.Clr();
  G_S.AddVMerged(Expanded_Vertex);
  return;
  // how to select edge to insert
}

// delta_S_bar 是 S 的邻居集合，但不包括 S 中的顶点
TIntV Master::GetBoundary(TIntV G_S, TIntV& delta_S_bar) {
  TIntV delta_S;
  delta_S.Clr();
  delta_S_bar.Clr();

  for (TIntV::TIter NI = G_S.BegI(); NI < G_S.EndI(); NI++) {
    TUNGraph::TNodeI Node = G->GetNI(*NI);
    for (int i = 0; i < Node.GetInDeg(); i++) {
      if (!G_S.IsIn(Node.GetNbrNId(i))) {
        delta_S.Add(Node.GetId());
        delta_S_bar.Add(Node.GetNbrNId(i));
      }
    }
  }
  delta_S.Merge();
  delta_S_bar.Merge();
  return delta_S;
}

std::string createString(const std::string& vcc_data, int k) {
  std::ostringstream oss;
  oss << vcc_data << "_k=" << k
      << "_algorithm=VCCE.kvcc";  // 使用ostringstream进行拼接
  return oss.str();               // 返回拼接后的字符串
}

void Master::Load_kvcc(TIntVIntV& kvcc_array, string vcc_data) {
  try {
    string kvcc_data_name = createString(vcc_data, this->k);
    cout << kvcc_data_name << endl;
    TFIn inFile(kvcc_data_name.c_str());
    kvcc_array.Load(inFile);
    cout << kvcc_array.Len() << endl;
  } catch (TPt<TExcept>) {
    cout << endl << "***kvcc result does not exist.***" << endl;
  }
  for (TIntVIntV::TIter TI = kvcc_array.BegI(); TI < kvcc_array.EndI(); TI++) {
    PUNGraph G_kvcc = TSnap::GetSubGraph(G, *TI);
    // debug_print("kvcc_nodes: " << G_kvcc->GetNodes() << " kvcc_edges: "
    //                            << G_kvcc->GetEdges() << endl);
  }
}

void Master::GroupSelection_single_vertex(
    TIntV& G_S, TIntV& delta_S, TIntV& delta_S_bar,
    unordered_set<pair<int, int>, pair_hash>& Inserted_Edge,
    TIntV& Expanded_Vertex) {
  cout << "kvcc: ";
  for (TIntV::TIter NI = G_S.BegI(); NI < G_S.EndI(); NI++) {
    cout << *NI << " ";
  }
  cout << endl;
  delta_S = GetBoundary(G_S, delta_S_bar);
  TIntV G_sub = delta_S;
  G_sub.AddVMerged(delta_S_bar);
  PUNGraph G_Cand = TSnap::GetSubGraph(G, G_sub);
  // cout << delta_S_bar.Len() << endl;

  TIntVIntV S;
  for (int i = 0; i <= k; i++) {
    S.Add({});
  }
  TIntV res = {};
  TIntV nb_u1 = {}, nb_u2 = {};

  TIntIntVH in_neighs, out_neighs;
  for (TIntV::TIter TI = delta_S_bar.BegI(); TI < delta_S_bar.EndI(); TI++) {
    nb_u1.Clr();
    nb_u2.Clr();
    TUNGraph::TNodeI v = G_Cand->GetNI(*TI);
    for (int i = 0; i < v.GetInDeg(); ++i) {
      if (delta_S.IsIn(v.GetInNId(i)))  // only consider neighbors in
        // k-vcc
        nb_u1.AddMerged(v.GetInNId(i));
      else
        nb_u2.AddMerged(v.GetInNId(i));
    }
    in_neighs.AddDat(*TI, nb_u1);
    out_neighs.AddDat(*TI, nb_u2);
    int index = nb_u1.Len();

    // impossible situation
    if (index >= k) index = k;
    S[index].AddMerged(*TI);
  }
  int level;
  for (int i = k; i > 0; i--) {
    debug_print("S[" << i << "]: ");
    level = i;
    for (int j = 0; j < S[i].Len(); j++) {
      int v = S[i][j];
      debug_print(v << " ");
    }
    debug_print(endl);
  }
  for (int i = k; i > 0;) {
    debug_print("S[" << i << "]: " << endl);
    level = i;
    for (int j = 0; j < S[i].Len(); j++) {
      int v = S[i][j];
      debug_print("Selected Vertex:" << v << " " << endl);
      int need = k - i;  // v need how much edges to insert
      for (TIntV::TIter TI = G_S.BegI(); TI < G_S.EndI(); TI++) {
        if (!in_neighs.GetDat(v).IsIn(*TI)) {
          if (acost >= b) {
            // cout << "acost: " << acost << endl;
            // cout << "gain: " << Expanded_Vertex.Len() << endl;
            return;
          }
          // insert_edges.Add({*TI, v});
          if (Inserted_Edge.find({*TI, v}) != Inserted_Edge.end()) continue;
          Inserted_Edge.insert({*TI, v});
          debug_print("(insert: " << *TI << " " << v << ") " << endl);
          // b--;

          need--;
          acost++;
          if (need == 0) {
            // v have k neighbors in kvcc
            Expanded_Vertex.Add(v);
            debug_print("expanded: " << v << endl);
            // v is expanded, update its neighbors
            update_neighbour(S, in_neighs, out_neighs, v, Expanded_Vertex,
                             level);
            break;
          }
        }
      }
      // cout<<S[i][j]<<" ";
    }
    // exist a vertex in S[level] that is not expanded
    if (level > i) {
      i = level;
    } else {
      i--;
    }
  }
  // cout << "acost: " << acost << endl;
  G_S.AddVMerged(Expanded_Vertex);
  return;
}

void Master::update_neighbour(TIntVIntV& S, TIntIntVH& in_neighs,
                              TIntIntVH& out_neighs, int v, TIntV& res,
                              int& level) {
  // res: new expanded vertices
  TIntV neigh = out_neighs.GetDat(v);
  // cout << neigh.Len() << endl;

  // cout << level << endl;
  for (TIntV::TIter NI = neigh.BegI(); NI < neigh.EndI(); NI++) {
    /*cout << v <<" " << *NI << endl;*/
    if (in_neighs.GetDat(*NI).Len() != 0) {
      out_neighs.GetDat(*NI).DelIfIn(v);
      in_neighs.GetDat(*NI).Add(v);
    }
  }
  for (TIntV::TIter NI = neigh.BegI(); NI < neigh.EndI(); NI++) {
    int idx = in_neighs.GetDat(*NI).Len() -
              1;  //上面的循环已经更新过，所以这里为获得idx初始值应该-1
    /*cout << idx << endl;*/
    // cout << "updated_idx: " << idx + 1 << " v:" << *NI << endl;
    // 回到上一层
    //避免越界
    idx = min(k - 1, idx);
    if (idx + 1 > level) {
      level = idx + 1;
    }
    if (idx < 1) continue;
    if (idx == k - 1) {
      S[idx].DelIfIn(*NI);
      // S[idx + 1].AddMerged(*NI);
      res.AddMerged(*NI);
      debug_print("expanded by update: " << *NI << endl);
      in_neighs.GetDat(*NI) = {};
      update_neighbour(S, in_neighs, out_neighs, *NI, res, level);
    } else if (idx < k - 1) {
      S[idx].DelIfIn(*NI);
      S[idx + 1].AddUnique(*NI);
    }
  }
}
vector<int> Master::get_common_set(TIntV& Vcc_1, TIntV& Vcc_2) {
  vector<int> common_neighbors;

  auto it1 = Vcc_1.BegI();
  auto it2 = Vcc_2.BegI();

  // both vcc_1 and vcc_2 are sorted
  while (it1 != Vcc_1.EndI() && it2 != Vcc_2.EndI()) {
    if (*it1 == *it2) {
      common_neighbors.push_back(*it1);
      ++it1;
      ++it2;
    } else if (*it1 < *it2) {
      ++it1;
    } else {
      ++it2;
    }
  }
  return common_neighbors;
}

std::pair<int, int> Master::find_max_in_com_neigh(
    const std::unordered_map<std::pair<int, int>, std::vector<int>, pair_hash>&
        com_neigh) {
  int max_value = std::numeric_limits<int>::min();  // Initialize to the
                                                    // smallest integer value
  std::pair<int, int> max_key;  // To record the key with the maximum value

  for (const auto& entry : com_neigh) {
    const std::pair<int, int>& key = entry.first;
    const int local_max = entry.second.size();

    // Update the global maximum value and corresponding key
    if (local_max > max_value) {
      max_value = local_max;
      max_key = key;
    }
  }

  // Return the key corresponding to the maximum value
  return max_key;
}

std::pair<int, int> Master::find_max_in_gamma(
    const std::unordered_map<std::pair<int, int>, int, pair_hash>& gamma) {
  int max_value = std::numeric_limits<int>::min();  // Initialize to the
                                                    // smallest integer value
  std::pair<int, int> max_key;  // To record the key with the maximum value

  for (const auto& entry : gamma) {
    const std::pair<int, int>& key = entry.first;
    const int local_max = entry.second;
    // Update the global maximum value and corresponding key
    if (local_max > max_value) {
      max_value = local_max;
      max_key = key;
    }
  }

  // Return the key corresponding to the maximum value
  return max_key;
}

TIntV Master::get_difference_set(TIntV& Vcc_1, TIntV& Vcc_2) {
  TIntV diff;
  for (int i = 0; i < Vcc_1.Len(); i++) {
    if (!Vcc_2.IsIn(Vcc_1[i])) {
      diff.Add(Vcc_1[i]);
    }
  }
  return diff;
}

// int Master::get_min_deg_vertex(TIntV& Vcc) {
//   int min_deg = std::numeric_limits<int>::max();
//   for (int i = 0; i < Vcc.size(); i++) {
//     int deg = G->GetNI(Vcc[i]).GetDeg();
//     if (deg < min_deg) {
//       min_deg = deg;
//     }
//   }
//   return min_deg;
// }

void Master::sort_by_deg(TIntV& Vcc) {
  sort(Vcc.BegI(), Vcc.EndI(), [this](int a, int b) {
    return G->GetNI(a).GetDeg() < G->GetNI(b).GetDeg();
  });
}

void Master::Merge_overlap_vcc(
    TIntVIntV& VCCs, unordered_set<pair<int, int>, pair_hash>& Inserted_Edge,
    TIntV& Expanded_Vertex) {
  std::unordered_map<pair<int, int>, vector<int>, pair_hash> com_neigh;
  cout << "VCCs: " << VCCs.Len() << endl;
  for (int i = 0; i < VCCs.Len(); i++) {
    for (int j = i + 1; j < VCCs.Len(); j++) {
      TIntV* VCC_i = &VCCs[i];
      TIntV* VCC_j = &VCCs[j];
      if (com_neigh.find(make_pair(i, j)) == com_neigh.end()) {
        com_neigh[make_pair(i, j)] = get_common_set(*VCC_i, *VCC_j);
      }
    }
    // cout << "com_neigh: " << com_neigh.size() << endl;
  }

  while (acost < b) {
    auto key = find_max_in_com_neigh(com_neigh);
    int i = key.first;
    int j = key.second;
    int p = com_neigh[key].size();
    cout << "i: " << i << " j: " << j << " p: " << p << endl;
    TIntV S_1 = get_difference_set(VCCs[i], VCCs[j]);
    TIntV S_2 = get_difference_set(VCCs[j], VCCs[i]);
    sort_by_deg(S_1);
    sort_by_deg(S_2);

    if (acost + k - p > b) {
      com_neigh.erase(key);
      continue;
    }
    int blank = k - p;

    int s2_idx = 0;
    for (int s1_idx = 0; s1_idx < S_1.Len() && blank > 0; s1_idx++) {
      int v = S_1[s1_idx];
      int u = S_2[s2_idx++];
      if (Inserted_Edge.find({v, u}) != Inserted_Edge.end() &&
          Inserted_Edge.find({u, v}) != Inserted_Edge.end())
        continue;
      Inserted_Edge.insert({v, u});
      blank -= 1;
      acost += 1;
      // S_2.DelIfIn(u);
      if (blank == 0) {
        com_neigh[key] = {};
        break;
      }
    }
  }
  cout << "acost: " << acost << " b: " << b << endl;
  debug_print("Inserted_Edge: " << Inserted_Edge.size() << endl);
  for (auto& edge : Inserted_Edge) {
    debug_print(edge.first << " " << edge.second << endl);
  }
  // debug_print("Expanded_Vertex: " << Expanded_Vertex.Len() << endl);
  return;
}

void Master::Merge_adjacent_vcc(
    TIntVIntV& VCCs, unordered_set<pair<int, int>, pair_hash>& Inserted_Edge,
    TIntV& Expanded_Vertex) {
  std::unordered_map<pair<int, int>, int, pair_hash> gamma;
  cout << "VCCs: " << VCCs.Len() << endl;
  for (int i = 0; i < VCCs.Len(); i++) {
    TIntV* VCC_i = &VCCs[i];
    for (int j = i + 1; j < VCCs.Len(); j++) {
      TIntV* VCC_j = &VCCs[j];
      // skip overlap vcc
      if (get_common_set(*VCC_i, *VCC_j).size() > 0) continue;
      int t = 0, w = 0;
      std::unordered_map<int, bool> flag;
      for (TIntV::TIter uI = VCC_i->BegI(); uI != VCC_i->EndI(); uI++) {
        int u = *uI;
        if (flag.find(u) != flag.end() && flag[u] == true) continue;
        for (int d = 0; d < G->GetNI(u).GetDeg(); d++) {
          int v = G->GetNI(u).GetNbrNId(d);
          if (flag.find(v) != flag.end() && flag[v] == true) continue;
          if (VCC_j->IsIn(v)) {
            t += 1;
            flag[u] = true;
            w += 1;
            flag[v] = true;
          }
        }
      }
      gamma[make_pair(i, j)] = min(t, w);
    }
  }

  while (acost < b) {
    auto key = find_max_in_gamma(gamma);
    int i = key.first;
    int j = key.second;
    int gamma_i_j = gamma[key];
    cout << "i: " << i << " j: " << j << " gamma_i_j: " << gamma_i_j << endl;

    TIntV N_i, N_j;
    TIntV S_1, S_2, S_3, Psi;

    GetBoundary(VCCs[i], N_i);
    // cout << "N_i: " << N_i.Len() <<endl;
    // for(TIntV::TIter v = N_i.BegI(); v!=N_i.EndI(); v++) {
    //   cout <<*v<<" ";
    // }
    // cout << endl;
    GetBoundary(VCCs[j], N_j);
    // cout << "N_j: " << N_j.Len() <<endl;
    // for(TIntV::TIter v = N_j.BegI(); v!=N_j.EndI(); v++) {
    //   cout <<*v<<" ";
    // }
    // cout << endl;

    if (gamma_i_j == get_common_set(VCCs[i], N_j).size()) {
      S_1 = get_difference_set(VCCs[i], N_j);
      Psi = VCCs[j];
      S_2 = get_difference_set(VCCs[j], N_i);
      vector<int> tmp = get_common_set(VCCs[j], N_i);
      for (auto u : tmp) {
        S_3.Add(u);
      }
    } else {
      S_1 = get_difference_set(VCCs[j], N_i);
      Psi = VCCs[i];
      S_2 = get_difference_set(VCCs[i], N_j);
      vector<int> tmp = get_common_set(VCCs[i], N_j);
      for (auto u : tmp) {
        S_3.Add(u);
      }
    }

    if (acost + k - gamma_i_j > b) {
      gamma.erase(key);
      continue;
    }
    int blank = k - gamma_i_j;
    sort_by_deg(S_1);
    sort_by_deg(S_2);
    sort_by_deg(S_3);
    // cout << "S_1: " << S_1.Len() <<endl;
    // for(TIntV::TIter v = S_1.BegI(); v!=S_1.EndI(); v++) {
    //   cout <<*v<<" ";
    // }
    // cout << endl;
    // cout << "S_2: " << S_2.Len() <<endl;
    // for(TIntV::TIter v = S_2.BegI(); v!=S_2.EndI(); v++) {
    //   cout <<*v<<" ";
    // }
    // cout << endl;
    // cout << "S_3: " << S_3.Len() <<endl;
    // for(TIntV::TIter v = S_3.BegI(); v!=S_3.EndI(); v++) {
    //   cout <<*v<<" ";
    // }
    // cout << endl;
    int s2_idx = 0, s3_idx = 0;
    if (blank == 0) {
      cout << "gamma_i_j >= k" << endl;
      gamma[key] = -1;
    }
    for (int s1_idx = 0; s1_idx < S_1.Len() && blank > 0; s1_idx++) {
      int v = S_1[s1_idx];
      int u;
      if (!S_2.Empty()) {
        u = S_2[s2_idx++];
        // S_2.DelIfIn(u);
      } else {
        u = S_3[s3_idx++];
        // S_3.DelIfIn(u);
      }
      if (Inserted_Edge.find({v, u}) != Inserted_Edge.end() &&
          Inserted_Edge.find({u, v}) != Inserted_Edge.end())
        continue;
      Inserted_Edge.insert({v, u});
      blank -= 1;
      acost += 1;

      // S_2.DelIfIn(u);
      if (blank <= 0) {
        gamma[key] = -1;
        break;
      }
    }
  }
  cout << "acost: " << acost << " b: " << b << endl;
  debug_print("Inserted_Edge: " << Inserted_Edge.size() << endl);
  for (auto& edge : Inserted_Edge) {
    debug_print(edge.first << " " << edge.second << endl);
  }
  debug_print("Expanded_Vertex: " << Expanded_Vertex.Len() << endl);
  return;
}

void Master::CalConnectKVcc(TIntVIntV& VCCs, vector<int>& I, vector<int>& J,
                            vector<vector<int>>& T, vector<double>& R) {
  std::priority_queue<std::pair<std::pair<int, int>, double>,
                      std::vector<std::pair<std::pair<int, int>, double>>,
                      Compare>
      r;
  std::unordered_map<pair<int, int>, vector<int>, pair_hash> t;
  cout << "VCCs: " << VCCs.Len() << endl;
  // TODO: sort VCCs by size
  VCCs.Sort();

  for (int i = 0; i < VCCs.Len(); i++) {
    TIntV* VCC_i = &VCCs[i];
    for (int j = i + 1; j < VCCs.Len(); j++) {
      TIntV* VCC_j = &VCCs[j];

      int cost = k, gain = 0;
      int p_ij = get_common_set(*VCC_i, *VCC_j).size();
      cost -= p_ij;

      TIntV N_i, N_j;
      GetBoundary(VCCs[i], N_i);
      GetBoundary(VCCs[j], N_j);
      TIntV diff_ij = get_difference_set(VCCs[i], VCCs[j]);
      TIntV diff_ji = get_difference_set(VCCs[j], VCCs[i]);
      vector<int> S_L = get_common_set(diff_ij, N_j);
      vector<int> S_R = get_common_set(diff_ji, N_i);

      vector<int> M;
      std::vector<std::pair<int, int>> matches;
      // cout << "i: " << i << ", j: " << j << endl;
      PUNGraph G_LR = RemoveInternalEdges(G, S_L, S_R);
      if (G_LR->GetEdges() == 0) continue;
      if (S_L.size() != 0) {
        int MaxMatch = HopcroftKarp(G_LR, S_L, S_R, matches);
        for (auto match : matches) {
          t[make_pair(i, j)].push_back(match.first);
          t[make_pair(j, i)].push_back(match.second);
        }
        cost -= MaxMatch;
      }

      if (cost < k) {
        gain = 2 * VCCs[i].Len() * VCCs[j].Len() + p_ij * VCCs[j].Len();
        r.push({make_pair(i, j), (double)gain / cost});
      }
    }
  }

  bool flag = true;
  while (flag && !r.empty()) {
    auto top = r.top();
    r.pop();
    int i_star = top.first.first, j_star = top.first.second,
        r_ij_star = top.second;
    std::cout << "Key: (" << i_star << ", " << j_star
              << "), Value: " << r_ij_star << std::endl;
    if (r_ij_star == 0) {
      flag = false;
      break;
    }
    I.push_back(i_star);
    T.push_back(t[make_pair(i_star, j_star)]);
    T.push_back(t[make_pair(j_star, i_star)]);
    J.push_back(j_star);
    R.push_back(r_ij_star);
  }
}

void Master::CalMulVerices(TIntVIntV& VCCs, vector<int>& I,
                           vector<vector<int>>& MC, vector<double>& R) {
  std::priority_queue<std::pair<std::pair<int, int>, double>,
                      std::vector<std::pair<std::pair<int, int>, double>>,
                      Compare>
      r;

  vector<vector<vector<int>>>
      allCliques;  // 新增数组，用于保存每个 VCC_i 计算得到的 cliques
  for (int i = 0; i < VCCs.Len(); i++) {
    vector<vector<int>>
        currentCliques;  // 用于存储当前计算得到的所有 cliques 转换后的结果

    TIntV VCC_i = VCCs[i];
    // calcuate one-hop neighbor of Graph G as candidate set
    TIntV delta_S, delta_S_bar;
    delta_S = GetBoundary(VCC_i, delta_S_bar);
    TIntV G_sub = delta_S;
    G_sub.AddVMerged(delta_S_bar);
    PUNGraph G_Cand = TSnap::GetSubGraph(G, G_sub);

    // calculate the edge num from vertex in delta_S_bar to k-vcc,
    // divide it into S_(k), S_(k-1), S_(k-2), ...
    TIntVIntV S;
    for (int i = 0; i <= k; i++) {
      S.Add({});
    }
    TIntV res = {};
    TIntV nb_u1 = {}, nb_u2 = {};

    // in_neighs stores the neighbors of each candidate node belonging to k-vcc
    // out_neighs stores the neighbors of each candidate node belonging to
    // delta_S_bar
    TIntIntVH in_neighs, out_neighs;
    for (TIntV::TIter TI = delta_S_bar.BegI(); TI < delta_S_bar.EndI(); TI++) {
      nb_u1.Clr();
      nb_u2.Clr();
      TUNGraph::TNodeI v = G_Cand->GetNI(*TI);
      for (int i = 0; i < v.GetInDeg(); ++i) {
        if (delta_S.IsIn(v.GetInNId(i)))  // only consider neighbors in
          // k-vcc
          nb_u1.AddMerged(v.GetInNId(i));
        else
          nb_u2.AddMerged(v.GetInNId(i));
      }
      in_neighs.AddDat(*TI, nb_u1);
      out_neighs.AddDat(*TI, nb_u2);
      int index = nb_u1.Len();

      // impossible situation
      // if (index >= k) index = k;
      S[index].AddMerged(*TI);
    }

    // calculate maximal cliques MC in each S
    int num = 0;
    TIntV S_total;
    TIntVIntV cliques;
    TIntVTIntVH cliques_neighs_union;

    // priority_queue<pair<vector<int>, pair<int, int>>,
    //                vector<pair<vector<int>, pair<int, int>>>, decltype(comp)>
    //     clique_in_neighs(comp);

    int level = 0;
    int clq_idx = 0;
    for (int layer = k - 1; layer > 0; layer--) {
      // debug_print("S[" << i << "]: ");
      TIntV S_layer = S[layer];
      TIntV neigh_Union;
      TIntV S_layer_temp;

      S_total.Clr();
      S_total.AddVMerged(S[layer]);
      // to be considered how to deal with?
      // if (S_total.Len() < k - i + 1) continue;
      PUNGraph sub_G = TSnap::GetSubGraph(G_Cand, S_total);
      // PUNGraph sub_core = TSnap::GetKCore(sub_G, k - i);
      // to be considered how to deal with clique size smaller than k - i + 1
      TCliqueOverlap::GetMaxCliques(sub_G, 0, cliques);
      // for (int clique_idx = 0; clique_idx < cliques.Len(); clique_idx++) {

      for (int j = 0; j < cliques.Len(); j++) {
        TIntV cq = cliques[j];
        neigh_Union.Clr();
        int flag = 0;
        vector<int> vec_cq;
        for (TIntV::TIter TI = cq.BegI(); TI < cq.EndI(); TI++) {
          vec_cq.push_back(*TI);
          neigh_Union.AddV(in_neighs.GetDat(*TI));
        }
        currentCliques.push_back(vec_cq);
        neigh_Union.Merge();
        cliques_neighs_union.AddDat(cq, neigh_Union);
        // clique_in_neighs.push({vec_cq, pai r<int, int>(i + vec_cq.size(),
        // i)});
        vec_cq.clear();

        int cost = k - neigh_Union.Len();
        if (cost == 0) {
          // TODO: 这里有问题，存在cost=0的点
          // cout << "false" << endl;
          // // 输出 neigh_Union 的值
          // cout << "neigh_Union for clique (" << i << ", " << j << "): ";
          // for (TIntV::TIter TI = neigh_Union.BegI(); TI < neigh_Union.EndI();
          // TI++) {
          //     cout << *TI << " ";
          // }
          // cout << endl;
          cost = INT_MAX;
        }
        double gain = (2 * VCC_i.Len() + cq.Len()) * cq.Len();
        r.push({make_pair(i, clq_idx + j), (double)gain / cost});
      }
      clq_idx += cliques.Len();
    }
    allCliques.push_back(currentCliques);
  }

  bool flag = true;
  while (flag && !r.empty()) {
    auto top = r.top();
    r.pop();
    int i_star = top.first.first, j_star = top.first.second,
        r_ij_star = top.second;
    std::cout << "Key: (" << i_star << ", " << j_star
              << "), Value: " << r_ij_star << std::endl;
    if (r_ij_star == 0) {
      flag = false;
      break;
    }
    I.push_back(i_star);
    MC.push_back(allCliques[i_star][j_star]);
    R.push_back(r_ij_star);
  }
}

void Master::ExpSinVertices(
    TIntVIntV& VCCs, double& r,
    unordered_set<pair<int, int>, pair_hash>& Inserted_Edge) {
  int best_cost, best_v;
  TIntV best_vcc, best_in_neighs;
  VCCs.Sort();
  TIntV VCC_max = VCCs.Last();
  for (int i = VCCs.Len() - 1; i >= 0; --i) {
    if (VCCs[i].Len() < (2 * VCC_max.Len() - k + 2) / (2 * k - 2)) {
      break;
    }
    TIntV Vcc_i = VCCs[i];
    TIntV delta_S, delta_S_bar;
    delta_S = GetBoundary(Vcc_i, delta_S_bar);
    TIntV G_sub = delta_S;
    G_sub.AddVMerged(delta_S_bar);
    PUNGraph G_Cand = TSnap::GetSubGraph(G, G_sub);
    TIntVIntV S;
    for (int i = 0; i <= k; i++) {
      S.Add({});
    }
    TIntV res = {};
    TIntV nb_u1 = {}, nb_u2 = {};

    TIntIntVH in_neighs, out_neighs;
    for (TIntV::TIter TI = delta_S_bar.BegI(); TI < delta_S_bar.EndI(); TI++) {
      nb_u1.Clr();
      nb_u2.Clr();
      TUNGraph::TNodeI v = G_Cand->GetNI(*TI);
      for (int i = 0; i < v.GetInDeg(); ++i) {
        if (delta_S.IsIn(v.GetInNId(i)))  // only consider neighbors in
          // k-vcc
          nb_u1.AddMerged(v.GetInNId(i));
        else
          nb_u2.AddMerged(v.GetInNId(i));
      }
      in_neighs.AddDat(*TI, nb_u1);
      out_neighs.AddDat(*TI, nb_u2);
      int deg = nb_u1.Len();

      int cost = k - deg;
      double gain = 2 * Vcc_i.Len() + 1;
      if (r < double(gain) / cost) {
        r = double(gain) / cost;
        best_cost = cost;
        best_v = *TI;
        best_vcc = Vcc_i;
        best_in_neighs = in_neighs.GetDat(best_v);
      }
    }
  }

  TIntV S_cand;
  for (TIntV::TIter TI = best_vcc.BegI(); TI < best_vcc.EndI(); TI++) {
    if (!best_in_neighs.IsIn(*TI)) {
      S_cand.Add(*TI);
    }
  }
  sort_by_deg(S_cand);
  acost += best_cost;
  while (best_cost > 0) {
    int u = S_cand[0];
    if (Inserted_Edge.find({u, best_v}) != Inserted_Edge.end()) continue;
    Inserted_Edge.insert({u, best_v});
    best_cost -= 1;
    S_cand.DelIfIn(u);
  }
}