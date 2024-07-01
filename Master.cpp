#include "Master.h"

Master::Master(PUNGraph G, int k, int b) {
    this->G = G;
    acost = 0;
    nfs = 0;
    this->b = b;
    this->k = k;
}

void Master::Anchoring(string alg) {
    int round = 0;
    double group_anchor_time = 0.0, vertex_anchor_time = 0.0;
    TIntVIntV kvcc_array;
    TIntV kvcc, delta_S, delta_S_bar;
    TIntV Expanded_Vertex;
    Load_kvcc(kvcc_array);
    // select kvcc to expand
    kvcc = kvcc_array[0];
    while (acost < b) {
        cout << " -- Anchoring round: " << round++ << endl;
        double vertex_begin = (double)clock();
        double node_score = 0;

        // Compute by Multiple Vertex Anchoring
        vector<double> group;
        if (alg == string("t")) {
            group =
                GroupSelection_together(kvcc, delta_S, delta_S_bar, Expanded_Vertex);
        }
        else if (alg == string("m")) {
            group =
                GroupSelection_multi_vertex(kvcc, delta_S, delta_S_bar, Expanded_Vertex);
        }
        else if (alg == string("s")) {
            group = GroupSelection_single_vertex(kvcc, delta_S, delta_S_bar,
                Expanded_Vertex);
        }
    }
}

vector<double> Master::GroupSelection_together(TIntV& G_S, TIntV& delta_S,
    TIntV& delta_S_bar,
    TIntV& Expanded_Vertex) {


}
vector<double> Master::GroupSelection_multi_vertex(TIntV& G_S, TIntV& delta_S,
    TIntV& delta_S_bar,
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
            if (delta_S.IsIn(v.GetInNId(i))) // only consider neighbors in
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

    // use priority queue to select the clique with max size + level (|c| + i)
    auto comp = [](const std::pair<std::vector<int>, std::pair<int, int>>& a,
        const std::pair<std::vector<int>, std::pair<int, int>>& b) {
            const auto& a_first = a.second.first;
            const auto& b_first = b.second.first;
            const auto& a_second = a.second.second;
            const auto& b_second = b.second.second;

            return (a_first == b_first) ? (a_second > b_second)
                : (a_first < b_first);
        };

    priority_queue<pair<vector<int>, pair<int, int>>,
        vector<pair<vector<int>, pair<int, int>>>, decltype(comp)>
        clique_in_neighs(comp);

    int level = 0;
    for (int i = k - 1; i > 0; i--) {
        cout << "S[" << i << "]: ";
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
            clique_in_neighs.push(
                { vec_cq, pair<int, int>(i + vec_cq.size(), i) });
            vec_cq.clear();
        }

        // (clique_in_neighs.empty()) continue;
        if (clique_in_neighs.empty())
            continue;
        vector<int> c = clique_in_neighs.top().first;
        TIntV c_tmp;
        int flag = 0;
        for (auto& v : c) {
            c_tmp.Add(v);
        }
        clique_in_neighs.pop();
        TIntV c_neighs = cliques_neighs_union.GetDat(c_tmp);
        for (auto& v : c) {
            if (flag == 1)
                break;
            int need = max(1, k - i - static_cast<int>(c.size()) +
                1); // eahc v need how much edges to insert,
            // 至少需要插入1条边
            cout << v << " " << endl;
            for (TIntV::TIter TI = G_S.BegI(); TI < G_S.EndI(); TI++) {
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
                    cout << "(insert: " << *TI << " " << v << ") " << endl;
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
                    cout << "(insert: " << *TI << " " << v << ") " << endl;
                    // b--;
                    need--;
                    acost++;
                    c_neighs.Add(*TI);
                }
            }
        }
        if (flag == 0 || flag == 1 &&
            k - i - static_cast<int>(c.size()) + 1 <= 0 &&
            c_neighs.Len() >= k) {
            for (auto& v : c) {
                // v have k neighbors in kvcc
                Expanded_Vertex.Add(v);
                cout << "expanded: " << v << endl;
                // v is expanded, update its neighbors
                update_neighbour(S, in_neighs, out_neighs, v, Expanded_Vertex,
                    level);
            }
        }
        if (flag == 1) {
            cout << "acost: " << acost << endl;
            cout << "gain: " << Expanded_Vertex.Len() << endl;
            terminate();
        }
    }
    cliques.Clr();
    // how to select edge to insert
}

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

void Master::Load_kvcc(TIntVIntV& kvcc_array) {
    try {
        TFIn inFile("CA-GrQc_k=5.kvcc");
        kvcc_array.Load(inFile);
        cout << kvcc_array.Len() << endl;
    }
    catch (TPt<TExcept>) {
        cout << endl << "***kvcc result does not exist.***" << endl;
    }
    for (TIntVIntV::TIter TI = kvcc_array.BegI(); TI < kvcc_array.EndI();
        TI++) {
        PUNGraph G_kvcc = TSnap::GetSubGraph(G, *TI);
        cout << "kvcc_nodes: " << G_kvcc->GetNodes()
            << " kvcc_edges: " << G_kvcc->GetEdges() << endl;
    }
}

vector<double> Master::GroupSelection_single_vertex(TIntV& G_S, TIntV& delta_S,
    TIntV& delta_S_bar,
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
            if (delta_S.IsIn(v.GetInNId(i))) // only consider neighbors in
                // k-vcc
                nb_u1.AddMerged(v.GetInNId(i));
            else
                nb_u2.AddMerged(v.GetInNId(i));
        }
        in_neighs.AddDat(*TI, nb_u1);
        out_neighs.AddDat(*TI, nb_u2);
        int index = nb_u1.Len();

        // impossible situation
        if (index >= k)
            index = k;
        S[index].AddMerged(*TI);
    }
    int level;
    for (int i = k; i > 0; i--) {
        cout << "S[" << i << "]: ";
        level = i;
        for (int j = 0; j < S[i].Len(); j++) {
            int v = S[i][j];
            cout << v << " ";
        }
        cout << endl;
    }
    for (int i = k; i > 0;) {
        cout << "S[" << i << "]: ";
        level = i;
        for (int j = 0; j < S[i].Len(); j++) {
            int v = S[i][j];
            cout << v << " ";
            int need = k - i; // v need how much edges to insert
            for (TIntV::TIter TI = G_S.BegI(); TI < G_S.EndI(); TI++) {
                if (!in_neighs.GetDat(v).IsIn(*TI)) {
                    if (acost >= b) {
                        cout << "acost: " << acost << endl;
                        cout << "gain: " << Expanded_Vertex.Len() << endl;
                        terminate();
                    }
                    // insert_edges.Add({*TI, v});
                    cout << "(insert: " << *TI << " " << v << ") " << endl;
                    // b--;

                    need--;
                    acost++;
                    if (need == 0) {
                        // v have k neighbors in kvcc
                        Expanded_Vertex.Add(v);
                        cout << "expanded: " << v << endl;
                        // v is expanded, update its neighbors
                        update_neighbour(S, in_neighs, out_neighs, v,
                            Expanded_Vertex, level);
                        break;
                    }
                }
            }
            // cout<<S[i][j]<<" ";
        }
        // exist a vertex in S[level] that is not expanded
        if (level > i) {
            i = level;
        }
        else {
            i--;
        }

        cout << endl;
    }
    cout << "acost: " << acost << endl;
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
            1; //上面的循环已经更新过，所以这里为获得idx初始值应该-1
        /*cout << idx << endl;*/
        cout << "updated_idx: " << idx + 1 << " v:" << *NI << endl;
        // 回到上一层
        if (idx + 1 > level) {
            level = idx + 1;
        }
        if (idx < 1)
            continue;
        if (idx == k - 1) {
            S[idx].DelIfIn(*NI);
            // S[idx + 1].AddMerged(*NI);
            res.AddMerged(*NI);

            in_neighs.GetDat(*NI) = {};
            update_neighbour(S, in_neighs, out_neighs, *NI, res, level);
        }
        else if (idx < k - 1) {
            S[idx].DelIfIn(*NI);
            S[idx + 1].AddMerged(*NI);
        }
    }
}