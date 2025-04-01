#include "Master.h"

#define _DEBUG
#ifdef _DEBUG
#define debug_print(msg) std::cout << msg
#else
#define debug_print(msg)
#endif

Master::Master(PUNGraph G, int k, int b, bool disable_cck, bool disable_cmv) {
  this->G = G;
  acost = 0;
  nfs = 0;
  this->b = b;
  this->k = k;
  this->disable_cck = disable_cck;
  this->disable_cmv = disable_cmv;
}

// void Exact_Anchoring(std::string alg, std::string vcc_data) {
//   double t_begin = (double)clock();
//   int round = 0;
//   double group_anchor_time = 0.0, vertex_anchor_time = 0.0;
//   TIntVIntV kvcc_array;
//   TIntV kvcc, delta_S, delta_S_bar;
//   TIntV Expanded_Vertex;
//   Load_kvcc(kvcc_array, vcc_data);
//   std::unordered_set<std::pair<int, int>, pair_hash> total_Inserted_Edge;
//   double total_gain = 0.0;
//   auto start_time = std::chrono::high_resolution_clock::now();

//   // 1. 用long
//   //
//   long计算图中所有点的k顶点连通能力值，每个点的k顶点连通能力值=他所属于的最大的那个kvcc的大小
//   int num_vertices = 0;
//   for (const auto& kvcc : kvcc_array) {
//     for (int vertex : kvcc) {
//       num_vertices = std::max(num_vertices, vertex + 1);
//     }
//   }
//   std::vector<int> original_connectivity =
//       calculateKVertexConnectivity(kvcc_array, num_vertices);

//   // 2. 枚举图中顶点对应完全图的边 - 图中所有的边， 作为候选集的边
//   std::unordered_set<std::pair<int, int>, pair_hash> candidate_edges;
//   for (int i = 0; i < num_vertices; ++i) {
//     for (int j = i + 1; j < num_vertices; ++j) {
//       if (G->IsEdge(i, j) == false) candidate_edges.insert({i, j});
//     }
//   }

//   // 3. 从候选集边中依次枚举，每次选择b条边，进行插入，并且计算收益
//   std::vector<std::unordered_set<std::pair<int, int>, pair_hash>>
//       edge_combinations = generateEdgeCombinations(candidate_edges, b);
//   double max_gain = -DBL_MAX;
//   std::unordered_set<std::pair<int, int>, pair_hash> best_edges;
//   for (const auto& combination : edge_combinations) {
//     double gain = calculateGain(kvcc_array, combination,
//     original_connectivity); if (gain > max_gain) {
//       max_gain = gain;
//       best_edges = combination;
//     }
//   }

//   // 4. 统计收益最大的b条边
//   total_Inserted_Edge = best_edges;
//   total_gain = max_gain;

//   auto end_time = std::chrono::high_resolution_clock::now();
//   auto duration =
//       std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time)
//           .count();

//   std::cout << "最终插入边数量: " << total_Inserted_Edge.size() << std::endl;
//   std::cout << "总的 gain: " << total_gain << std::endl;
//   std::cout << "算法耗时: " << duration << " 秒" << std::endl;
// }

void Master::Anchoring(std::string alg, std::string vcc_data,
                       double threshold) {
  double t_begin = (double)clock();
  int round = 0;
  double group_anchor_time = 0.0, vertex_anchor_time = 0.0;
  TIntVIntV kvcc_array;
  TIntV kvcc, delta_S, delta_S_bar;
  TIntV Expanded_Vertex;
  Load_kvcc(kvcc_array, vcc_data);
  std::unordered_set<std::pair<int, int>, pair_hash> total_Inserted_Edge;
  double total_gain = 0.0;

  // 用于记录所有插入的边
  std::vector<std::pair<int, int>> allInsertedEdges;

  auto start_time = std::chrono::high_resolution_clock::now();

  // 初始化各函数总耗时
  double total_cal_connect_kvcc_time = 0.0;
  double total_cal_mul_verices_time = 0.0;
  double total_exp_sin_vertices_time = 0.0;
  double total_mer_connect_kvcc_time = 0.0;
  double total_exp_mul_vertices_time = 0.0;
  double total_time = 0.0;

  // 记录上一轮的 ra 和 rb 值
  double prev_ra = 0.0;
  double prev_rb = 0.0;

  // 记录上一轮 CalConnectKVcc 的结果
  int prev_Io = 0;
  int prev_J = 0;
  std::vector<std::vector<int>> prev_T;
  double prev_Ro = 0.0;
  int prev_Co = 0;
  double prev_gain_mer = 0.0;

  int initial_b = b;
  double quarter_b = initial_b * 0.25;
  double half_b = initial_b * 0.5;
  double three_quarters_b = initial_b * 0.75;

  double quarter_time = 0.0;
  double half_time = 0.0;
  double three_quarters_time = 0.0;
  double full_time = 0.0;

  double quarter_gain = 0.0;
  double half_gain = 0.0;
  double three_quarters_gain = 0.0;
  double full_gain = 0.0;

  std::cout << "开始循环，初始 b: " << b << std::endl;
  while (b > 0) {
    int Io, Ie, J, Co, Ce;
    double Ro = 0, Re = 0;
    std::vector<int> MC;
    std::vector<std::vector<int>> T;
    std::unordered_set<std::pair<int, int>, pair_hash> Inserted_Edge;
    double gain_sin, gain_mul, gain_mer;

    bool use_optimized_strategy = (alg == "approximate");
    bool can_execute_cal_connect_kvcc = true;

    if (use_optimized_strategy) {
      can_execute_cal_connect_kvcc =
          (round == 0) || (prev_ra > prev_rb * threshold);
    }

    if (can_execute_cal_connect_kvcc) {
      // 记录 CalConnectKVcc 开始时间
      auto start_cal_connect_kvcc = std::chrono::high_resolution_clock::now();
      if (!disable_cck) {
        CalConnectKVcc(kvcc_array, Io, J, T, Ro, Co, gain_mer);
      } else {
        Io = 0;
        J = 0;
        T.clear();
        Ro = 0.0;
        Co = 0;
        gain_mer = 0.0;
      }
      // 记录 CalConnectKVcc 结束时间
      auto end_cal_connect_kvcc = std::chrono::high_resolution_clock::now();
      auto duration_cal_connect_kvcc =
          std::chrono::duration_cast<std::chrono::seconds>(
              end_cal_connect_kvcc - start_cal_connect_kvcc)
              .count();
      total_cal_connect_kvcc_time += duration_cal_connect_kvcc;

      // 更新上一轮的结果
      prev_Io = Io;
      prev_J = J;
      prev_T = T;
      prev_Ro = Ro;
      prev_Co = Co;
      prev_gain_mer = gain_mer;
    } else {
      // 沿用上一轮的值
      Io = prev_Io;
      J = prev_J;
      T = prev_T;
      Ro = prev_Ro;
      Co = prev_Co;
      gain_mer = prev_gain_mer;
    }

    // std::cout << "T size after CalConnectKVcc: " << T.size() << std::endl;

    // 记录 CalMulVerices 开始时间
    auto start_cal_mul_verices = std::chrono::high_resolution_clock::now();
    if (!disable_cmv) {
      CalMulVerices(kvcc_array, Ie, MC, Re, Ce, gain_mul);
    } else {
      Ie = 0;
      MC.clear();
      Re = 0.0;
      Ce = 0;
      gain_mul = 0.0;
    }
    // CalMulVerices(kvcc_array, Ie, MC, Re, Ce, gain_mul);
    // 记录 CalMulVerices 结束时间
    auto end_cal_mul_verices = std::chrono::high_resolution_clock::now();
    auto duration_cal_mul_verices =
        std::chrono::duration_cast<std::chrono::seconds>(end_cal_mul_verices -
                                                         start_cal_mul_verices)
            .count();
    total_cal_mul_verices_time += duration_cal_mul_verices;

    double ra = Ro;
    double rb = Re;
    double rc = 0;
    int v = -1;
    int single_vcc_idx = -1;

    // 记录 ExpSinVertices 开始时间
    auto start_exp_sin_vertices = std::chrono::high_resolution_clock::now();
    std::unordered_set<std::pair<int, int>, pair_hash> Inserted_Edge_single;
    ExpSinVertices(kvcc_array, rc, Inserted_Edge_single, v, single_vcc_idx,
                   gain_sin);
    // 记录 ExpSinVertices 结束时间
    auto end_exp_sin_vertices = std::chrono::high_resolution_clock::now();
    auto duration_exp_sin_vertices =
        std::chrono::duration_cast<std::chrono::seconds>(end_exp_sin_vertices -
                                                         start_exp_sin_vertices)
            .count();
    total_exp_sin_vertices_time += duration_exp_sin_vertices;

    std::cout << "当前 ra: " << ra << ", rb: " << rb << ", rc: " << rc
              << std::endl;

    double current_gain = 0.0;
    if (ra == 0 && rb == 0 && rc == 0) {
      break;
    }
    if (ra >= rb && ra >= rc) {
      std::cout << "选择 ra 分支" << std::endl;
      if (!T.empty()) {
        int idx = Io;
        int jdx = J;
        std::vector<int> t_ij = T.front();
        T.assign(T.begin() + 1, T.end());
        std::vector<int> t_ji = T.front();
        T.assign(T.begin() + 1, T.end());
        int mer_cost = Co;

        // 记录 MerConnectKVcc 开始时间
        auto start_mer_connect_kvcc = std::chrono::high_resolution_clock::now();
        std::unordered_set<std::pair<int, int>, pair_hash> Inserted_Edge_mc;
        MerConnectKVcc(kvcc_array[idx], kvcc_array[jdx], t_ij, t_ji, ra,
                       mer_cost, Inserted_Edge_mc);
        // 记录 MerConnectKVcc 结束时间
        auto end_mer_connect_kvcc = std::chrono::high_resolution_clock::now();
        auto duration_mer_connect_kvcc =
            std::chrono::duration_cast<std::chrono::seconds>(
                end_mer_connect_kvcc - start_mer_connect_kvcc)
                .count();
        total_mer_connect_kvcc_time += duration_mer_connect_kvcc;

        current_gain = gain_mer;
        std::cout << "MerConnectKVcc 执行完毕，插入边数量: "
                  << Inserted_Edge_mc.size() << std::endl;
        if (b >= static_cast<int>(Inserted_Edge_mc.size())) {
          total_Inserted_Edge.insert(Inserted_Edge_mc.begin(),
                                     Inserted_Edge_mc.end());
          // 将插入的边记录到 allInsertedEdges 中
          for (const auto& edge : Inserted_Edge_mc) {
            allInsertedEdges.push_back(edge);
          }
          b -= Inserted_Edge_mc.size();
          InsertEdgesIntoGraph(Inserted_Edge_mc);
          // 将 idx 和 jdx 对应的 kvcc 合并
          TIntV merged_kvcc = kvcc_array[idx];
          for (TIntV::TIter it = kvcc_array[jdx].BegI();
               it < kvcc_array[jdx].EndI(); ++it) {
            if (!merged_kvcc.IsIn(*it)) {
              merged_kvcc.Add(*it);
            }
          }
          kvcc_array[idx] = merged_kvcc;
          // 移除 jdx 对应的 kvcc
          kvcc_array.Del(jdx);
          std::cout << "插入边成功，剩余 b: " << b << std::endl;
        } else {
          std::cout << "成本不足，程序结束" << std::endl;
          break;
        }
      } else {
        std::cout << "T 为空，跳过此分支" << std::endl;
      }
    } else if (rb >= ra && rb >= rc) {
      std::cout << "选择 rb 分支" << std::endl;
      if (!MC.empty()) {
        int idx = Ie;
        std::vector<int> mc_j = MC;
        std::unordered_set<std::pair<int, int>, pair_hash> Inserted_Edge_emv;
        TIntV mc_j_TIntV;
        for (int vertex : mc_j) {
          mc_j_TIntV.Add(vertex);
        }

        // 记录 ExpMulVertices 开始时间
        auto start_exp_mul_vertices = std::chrono::high_resolution_clock::now();
        ExpMulVertices(kvcc_array[idx], mc_j_TIntV, rb, Inserted_Edge_emv);
        // 记录 ExpMulVertices 结束时间
        auto end_exp_mul_vertices = std::chrono::high_resolution_clock::now();
        auto duration_exp_mul_vertices =
            std::chrono::duration_cast<std::chrono::seconds>(
                end_exp_mul_vertices - start_exp_mul_vertices)
                .count();
        total_exp_mul_vertices_time += duration_exp_mul_vertices;

        current_gain = gain_mul;
        std::cout << "ExpMulVertices 执行完毕，插入边数量: "
                  << Inserted_Edge_emv.size() << std::endl;
        if (b >= static_cast<int>(Inserted_Edge_emv.size())) {
          total_Inserted_Edge.insert(Inserted_Edge_emv.begin(),
                                     Inserted_Edge_emv.end());
          // 将插入的边记录到 allInsertedEdges 中
          for (const auto& edge : Inserted_Edge_emv) {
            allInsertedEdges.push_back(edge);
          }
          b -= Inserted_Edge_emv.size();
          InsertEdgesIntoGraph(Inserted_Edge_emv);
          // 将 mc_j_TIntV 中的点添加到 kvcc_array[idx]
          for (TIntV::TIter it = mc_j_TIntV.BegI(); it < mc_j_TIntV.EndI();
               ++it) {
            if (!kvcc_array[idx].IsIn(*it)) {
              kvcc_array[idx].Add(*it);
            }
          }
          std::cout << "插入边成功，剩余 b: " << b << std::endl;
        } else {
          std::cout << "成本不足，程序结束" << std::endl;
          break;
        }
      } else {
        std::cout << "MC 为空，跳过此分支" << std::endl;
      }
    } else {
      current_gain = gain_sin;
      std::cout << "选择 rc 分支" << std::endl;
      if (b >= static_cast<int>(Inserted_Edge_single.size())) {
        total_Inserted_Edge.insert(Inserted_Edge_single.begin(),
                                   Inserted_Edge_single.end());
        // 将插入的边记录到 allInsertedEdges 中
        for (const auto& edge : Inserted_Edge_single) {
          allInsertedEdges.push_back(edge);
        }
        b -= Inserted_Edge_single.size();
        InsertEdgesIntoGraph(Inserted_Edge_single);
        // 将顶点 v 加到下标为 single_vcc_idx 的 kvcc 中
        if (!kvcc_array[single_vcc_idx].IsIn(v)) {
          kvcc_array[single_vcc_idx].Add(v);
        }
        std::cout << "插入边成功，剩余 b: " << b << std::endl;
      } else {
        std::cout << "成本不足，程序结束" << std::endl;
        break;
      }
    }
    total_gain += current_gain;
    round++;
    std::cout << "第 " << round << " 轮循环结束，当前总 gain: " << total_gain
              << std::endl;

    // 更新上一轮的 ra 和 rb 值
    prev_ra = ra;
    prev_rb = rb;

    // 检查 b 的消耗情况
    if (initial_b - b >= quarter_b && quarter_time == 0) {
      auto current_time = std::chrono::high_resolution_clock::now();
      quarter_time = std::chrono::duration_cast<std::chrono::seconds>(
                         current_time - start_time)
                         .count();
      quarter_gain = total_gain;
    }
    if (initial_b - b >= half_b && half_time == 0) {
      auto current_time = std::chrono::high_resolution_clock::now();
      half_time = std::chrono::duration_cast<std::chrono::seconds>(
                      current_time - start_time)
                      .count();
      half_gain = total_gain;
    }
    if (initial_b - b >= three_quarters_b && three_quarters_time == 0) {
      auto current_time = std::chrono::high_resolution_clock::now();
      three_quarters_time = std::chrono::duration_cast<std::chrono::seconds>(
                                current_time - start_time)
                                .count();
      three_quarters_gain = total_gain;
    }
    if (b == 0 && full_time == 0) {
      auto current_time = std::chrono::high_resolution_clock::now();
      full_time = std::chrono::duration_cast<std::chrono::seconds>(
                      current_time - start_time)
                      .count();
      full_gain = total_gain;
    }
  }

  auto current_time = std::chrono::high_resolution_clock::now();
  total_time = std::chrono::duration_cast<std::chrono::seconds>(current_time -
                                                                start_time)
                   .count();

  std::cout << "循环结束，最终插入边数量: " << total_Inserted_Edge.size()
            << std::endl;
  std::cout << "最终 b: " << b << std::endl;
  std::cout << "总的 gain: " << total_gain << std::endl;

  // 输出各函数总耗时
  std::cout << "CalConnectKVcc 总耗时: " << total_cal_connect_kvcc_time << " 秒"
            << std::endl;
  std::cout << "CalMulVerices 总耗时: " << total_cal_mul_verices_time << " 秒"
            << std::endl;
  std::cout << "ExpSinVertices 总耗时: " << total_exp_sin_vertices_time << " 秒"
            << std::endl;
  std::cout << "MerConnectKVcc 总耗时: " << total_mer_connect_kvcc_time << " 秒"
            << std::endl;
  std::cout << "ExpMulVertices 总耗时: " << total_exp_mul_vertices_time << " 秒"
            << std::endl;

  // 输出 b 消耗 25%、50%、75%、100% 时的时间和 gain
  std::cout << "b 消耗 25% 时，时间: " << quarter_time
            << " 秒，gain: " << quarter_gain << std::endl;
  std::cout << "b 消耗 50% 时，时间: " << half_time
            << " 秒，gain: " << half_gain << std::endl;
  std::cout << "b 消耗 75% 时，时间: " << three_quarters_time
            << " 秒，gain: " << three_quarters_gain << std::endl;
  std::cout << "b 消耗 100% 时，时间: " << full_time
            << " 秒，gain: " << full_gain << std::endl;
  std::cout << "总耗时，时间: " << total_time << " 秒，gain: " << total_gain
            << std::endl;

  // 按行输出所有插入的边
  std::cout << "所有插入的边：" << std::endl;
  for (const auto& edge : allInsertedEdges) {
    std::cout << edge.first << " " << edge.second << std::endl;
  }
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
//   // kvcc = kvcc_array[1];
//   // unordered_set<pair<int, int>, pair_hash> Inserted_Edge;
//   vector<int> Io, Ie, J;
//   vector<double> Ro, Re;
//   vector<vector<int>> T, MC;
//   unordered_set<pair<int, int>, pair_hash> Inserted_Edge;

//   CalConnectKVcc(kvcc_array, Io, J, T, Ro);
//   CalMulVerices(kvcc_array, Ie, MC, Re);

//   double ra = *Ro.begin();
//   double rb = *Re.begin();
//   double rc;
//   unordered_set<pair<int, int>, pair_hash> Inserted_Edge_single;
//   unordered_set<pair<int, int>, pair_hash> Inserted_Edge_multi;
//   ExpSinVertex(kvcc_array, rc, Inserted_Edge_single);

//   if (ra >= rb && ra >= rc) {
//     alg = "m";
//   }

//   // while (acost < b) {
//   //   if (!Ro.empty() && !Re.empty()) {
//   //     int ra = *Ro.begin();
//   //     int rb = *Re.begin();
//   //     ExpSinVertex(kvcc_array, r, Inserted_Edge);
//   //     // std::cout << "ra: " << ra << ", rb: " << rb << std::endl;
//   //   } else {
//   //     break;
//   //   }
//   // }

//   // while (acost < b) {
//   //   cout << " -- Anchoring round: " << round++ << endl;
//   //   double vertex_begin = (double)clock();
//   //   double node_score = 0;

//   //   // Compute by Multiple Vertex Anchoring
//   //   // vector<double> group;

//   // }

//   // cout << "acost: " << acost << endl;
//   // cout << "gain: " << Expanded_Vertex.Len() << endl;
//   // cout << "Expanded_Vertex:";
//   // for (TIntV::TIter NI = Expanded_Vertex.BegI(); NI <
//   Expanded_Vertex.EndI();
//   //      NI++) {
//   //   cout << *NI << " ";
//   // }
//   // cout << endl;
//   // double t_end = (double)clock();
//   // cout << "the anchoring time is:" << (t_end - t_begin) / CLOCKS_PER_SEC
//   <<
//   // "s."
//   //      << endl;
// }

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

void Master::CalConnectKVcc(TIntVIntV& VCCs, int& i_star, int& j_star,
                            std::vector<std::vector<int>>& t_star,
                            double& r_star, int& cost_star, double& best_gain) {
  std::priority_queue<std::tuple<int, int, int, double>,
                      std::vector<std::tuple<int, int, int, double>>, Compare>
      r;
  std::unordered_map<std::pair<int, int>, std::vector<int>, pair_hash> t;
  std::cout << "VCCs: " << VCCs.Len() << std::endl;
  // TODO: sort VCCs by size
  VCCs.Sort();
  r_star = 0;

// 并行化外层循环
#pragma omp parallel for shared(VCCs, r_star, i_star, j_star, t_star, \
                                cost_star, best_gain, t)
  for (int i = 0; i < VCCs.Len(); i++) {
    TIntV* VCC_i = &VCCs[i];
    for (int j = i + 1; j < VCCs.Len(); j++) {
      TIntV* VCC_j = &VCCs[j];

      int cost = k;
      int p_ij = get_common_set(*VCC_i, *VCC_j).size();

      if (2 * VCCs[i].Len() * VCCs[j].Len() + p_ij * VCCs[j].Len() < r_star)
        continue;
      cost -= p_ij;

      TIntV N_i, N_j;
      GetBoundary(VCCs[i], N_i);
      GetBoundary(VCCs[j], N_j);
      TIntV diff_ij = get_difference_set(VCCs[i], VCCs[j]);
      TIntV diff_ji = get_difference_set(VCCs[j], VCCs[i]);
      std::vector<int> S_L = get_common_set(diff_ij, N_j);
      std::vector<int> S_R = get_common_set(diff_ji, N_i);

      std::vector<int> M;
      std::vector<std::pair<int, int>> matches;
      // std::cout << "i: " << i << ", j: " << j << std::endl;
      PUNGraph G_LR = RemoveInternalEdges(G, S_L, S_R);
      if (G_LR->GetEdges() == 0 || S_L.size() == 0) continue;

      int MaxMatch = HopcroftKarp(G_LR, S_L, S_R, matches);
      for (auto match : matches) {
#pragma omp critical
        {
          t[std::make_pair(i, j)].push_back(match.first);
          t[std::make_pair(j, i)].push_back(match.second);
        }
      }
      cost -= MaxMatch;

      if (cost < k) {
        double gain = 2 * VCCs[i].Len() * VCCs[j].Len() + p_ij * VCCs[j].Len();
        double r_ij = (double)gain / cost;
#pragma omp critical
        {
          if (r_star < r_ij && b >= cost) {
            i_star = i;
            j_star = j;
            t_star.clear();
            t_star.push_back(t[std::make_pair(i, j)]);
            t_star.push_back(t[std::make_pair(j, i)]);
            r_star = r_ij;
            cost_star = cost;
            best_gain = gain;
          }
        }
      }
    }
  }
}

void Master::CalMulVerices(TIntVIntV& VCCs, int& i_star,
                           std::vector<int>& mc_star, double& r_star,
                           int& cost_star, double& best_gain) {
  std::priority_queue<std::tuple<int, int, int, double>,
                      std::vector<std::tuple<int, int, int, double>>, Compare>
      r;

  std::vector<std::vector<std::vector<int>>>
      allCliques;  // 新增数组，用于保存每个 VCC_i 计算得到的 cliques

// 并行化外层循环
#pragma omp parallel for shared(VCCs, allCliques, i_star, mc_star, r_star, \
                                cost_star, best_gain)
  for (int i = 0; i < VCCs.Len(); i++) {
    std::vector<std::vector<int>>
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
    double local_r_star = 0;  // 每个线程的局部 r_star
    int local_i_star = 0;
    std::vector<int> local_mc_star;
    int local_cost_star = 0;
    double local_best_gain = 0;

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
      if (index >= k) index = k;
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
      // debug_print("S[" << layer << "]: ");
      // for (auto it = S[layer].BegI(); it != S[layer].EndI(); ++it)
      // {
      //   cout<<*it<<endl;
      // }
      // cout << endl;

      TIntV S_layer = S[layer];
      TIntV neigh_Union;
      TIntV S_layer_temp;

      S_total.Clr();
      S_total.AddVMerged(S[layer]);
      // to be considered how to deal with?
      if (S_total.Len() < k - i + 1) continue;
      PUNGraph sub_G = TSnap::GetSubGraph(G_Cand, S_total);
      // PUNGraph sub_core = TSnap::GetKCore(sub_G, k - i);
      // to be considered how to deal with clique size smaller than k - i + 1
      TCliqueOverlap::GetMaxCliques(sub_G, 0, cliques);
      // for (int clique_idx = 0; clique_idx < cliques.Len(); clique_idx++) {

      for (int j = 0; j < cliques.Len(); j++) {
        TIntV cq = cliques[j];
        neigh_Union.Clr();
        int flag = 0;
        std::vector<int> vec_cq;
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
          cost = INT_MAX;
        }
        double gain = (2 * VCC_i.Len() + cq.Len()) * cq.Len();
        double r_ij = (double)gain / cost;
        if (local_r_star < r_ij && b >= cost) {
          local_i_star = i;
          local_r_star = r_ij;
          local_mc_star.clear();
          local_best_gain = gain;
          // 将 cq 中的元素逐个添加到 mc_star 中
          for (TIntV::TIter TI = cq.BegI(); TI < cq.EndI(); ++TI) {
            local_mc_star.push_back(*TI);
          }
          local_cost_star = cost;
        }
      }
      clq_idx += cliques.Len();
    }

// 临界区更新全局变量
#pragma omp critical
    {
      if (r_star < local_r_star) {
        r_star = local_r_star;
        i_star = local_i_star;
        mc_star = local_mc_star;
        cost_star = local_cost_star;
        best_gain = local_best_gain;
      }
      allCliques.push_back(currentCliques);
    }
  }
}

void Master::ExpSinVertices(
    TIntVIntV& VCCs, double& r,
    std::unordered_set<std::pair<int, int>, pair_hash>& Inserted_Edge, int& v,
    int& vcc_idx, double& best_gain) {
  int best_cost = 0, best_v = 0, best_i = 0;
  TIntV best_vcc, best_in_neighs;
  VCCs.Sort();
  TIntV VCC_max = VCCs.Last();

  // Parallelize the outer loop over VCCs
  #pragma omp parallel
  {
    // Local variables for each thread to avoid race conditions
    double local_r = r;
    int local_best_cost = best_cost;
    int local_best_v = best_v;
    int local_best_i = best_i;
    TIntV local_best_vcc;
    TIntV local_best_in_neighs;
    double local_best_gain = best_gain;

    #pragma omp for nowait
    for (int i = VCCs.Len() - 1; i >= 0; --i) {
      if (VCCs[i].Len() < (2 * VCC_max.Len() - k + 2) / (2 * k - 2)) {
        continue;
      }
      TIntV Vcc_i = VCCs[i];
      TIntV delta_S, delta_S_bar;
      delta_S = GetBoundary(Vcc_i, delta_S_bar);
      TIntV G_sub = delta_S;
      G_sub.AddVMerged(delta_S_bar);
      PUNGraph G_Cand = TSnap::GetSubGraph(G, G_sub);

      // Inner loop over delta_S_bar can also be parallelized
      #pragma omp parallel for
      for (int j = 0; j < delta_S_bar.Len(); ++j) {
        TInt u = delta_S_bar[j];
        TIntV nb_u1, nb_u2;
        TUNGraph::TNodeI v_node = G_Cand->GetNI(u);
        for (int k = 0; k < v_node.GetInDeg(); ++k) {
          if (delta_S.IsIn(v_node.GetInNId(k))) {
            nb_u1.AddMerged(v_node.GetInNId(k));
          } else {
            nb_u2.AddMerged(v_node.GetInNId(k));
          }
        }
        int deg = nb_u1.Len();
        int cost = k - deg;
        double gain = 2 * Vcc_i.Len() + 1;

        // Thread-local comparison to avoid locks
        if (local_r < gain / cost && cost <= b) {
          local_r = gain / cost;
          local_best_cost = cost;
          local_best_v = u;
          local_best_vcc = Vcc_i;
          local_best_i = i;
          local_best_in_neighs = nb_u1; // Assuming nb_u1 is the in_neighs
          local_best_gain = gain;
        }
      }
    }

    // Critical section to update global best values
    #pragma omp critical
    {
      if (local_r > r) {
        r = local_r;
        best_cost = local_best_cost;
        best_v = local_best_v;
        best_vcc = local_best_vcc;
        best_i = local_best_i;
        best_in_neighs = local_best_in_neighs;
        best_gain = local_best_gain;
      }
    }
  }

  // Non-parallel part: Process the best candidate
  TIntV S_cand;
  for (TIntV::TIter TI = best_vcc.BegI(); TI < best_vcc.EndI(); TI++) {
    if (!best_in_neighs.IsIn(*TI)) {
      S_cand.Add(*TI);
    }
  }

  std::cout << "S_cand size: " << S_cand.Len() << std::endl;
  sort_by_deg(S_cand);

  if (S_cand.Empty()) {
    std::cout << "S_cand is empty! Cannot access S_cand[0]." << std::endl;
  } else {
    while (best_cost > 0) {
      int u = S_cand[0];
      if (Inserted_Edge.find({u, best_v}) != Inserted_Edge.end()) continue;
      Inserted_Edge.insert({u, best_v});
      best_cost -= 1;
      S_cand.DelIfIn(u);
    }
  }
  v = best_v;
  vcc_idx = best_i;
}

// 获取顶点在指定集合中的度数
int Master::get_degree_in_set(int vertex, const TIntV& set) {
  int degree = 0;
  TUNGraph::TNodeI node = G->GetNI(vertex);
  for (int i = 0; i < node.GetOutDeg(); ++i) {
    int neighbor = node.GetNbrNId(i);
    if (set.IsIn(neighbor)) {
      degree++;
    }
  }
  return degree;
}

void Master::MerConnectKVcc(
    TIntV& VCC_i, TIntV& VCC_j, std::vector<int> t_ij, std::vector<int> t_ji,
    double r_ij, int cost,
    std::unordered_set<std::pair<int, int>, pair_hash>& Inserted_Edge) {
  TIntV VCC_i_temp, VCC_j_temp;
  VCC_i_temp = VCC_i;
  for (int element : t_ji) {
    VCC_i_temp.Add(element);
  }

  VCC_j_temp = VCC_j;
  for (int element : t_ij) {
    VCC_j_temp.Add(element);
  }

  TIntV S_L = get_difference_set(VCC_i, VCC_j_temp);
  TIntV S_R = get_difference_set(VCC_j, VCC_i_temp);

  // 打印 VCC_i 和 VCC_j 的长度以及 r_ij 的值
  std::cout << "MerConnectKVcc: VCC_i size: " << VCC_i.Len()
            << ", VCC_j size: " << VCC_j.Len() << ", r_ij: " << r_ij
            << std::endl;

  // int cost = 2 * (VCC_i.Len() + VCC_j.Len()) / r_ij;

  std::cout << "MerConnectKVcc: S_L size: " << S_L.Len()
            << ", S_R size: " << S_R.Len() << ", cost: " << cost << std::endl;

  while (cost > 0) {
    if (S_L.Empty() || S_R.Empty()) {
      break;
    }

    // 从 S_L 中选取顶点 v，使得 v 在 VCC_i 中的度数最小
    int v = -1;
    int min_degree = std::numeric_limits<int>::max();
    for (int i = 0; i < S_L.Len(); ++i) {
      int degree = get_degree_in_set(S_L[i], VCC_i);
      if (degree < min_degree) {
        min_degree = degree;
        v = S_L[i];
      }
    }

    // 从 S_R 中选取顶点 u，使得 u 在 VCC_j 中的度数最小
    int u = -1;
    min_degree = std::numeric_limits<int>::max();
    for (int i = 0; i < S_R.Len(); ++i) {
      int degree = get_degree_in_set(S_R[i], VCC_j);
      if (degree < min_degree) {
        min_degree = degree;
        u = S_R[i];
      }
    }

    // 将边 (v, u) 加入到 Inserted_Edge 中
    Inserted_Edge.insert({v, u});
    std::cout << "MerConnectKVcc: Inserted edge (" << v << ", " << u << ")"
              << std::endl;

    // cost 减 1
    cost--;

    // 从 S_L 中移除顶点 v，从 S_R 中移除顶点 u
    S_L.DelIfIn(v);
    S_R.DelIfIn(u);
  }

  // 步骤 5: 返回插入的边集合
  // return Inserted_Edge;
}

TIntV Master::getNeighborhood(const TIntV& set, const TIntV& targetSet) {
  TIntV neighborhood;
  for (int i = 0; i < set.Len(); ++i) {
    TUNGraph::TNodeI node = G->GetNI(set[i]);
    for (int j = 0; j < node.GetOutDeg(); ++j) {
      int neighbor = node.GetNbrNId(j);
      if (targetSet.IsIn(neighbor) && !neighborhood.IsIn(neighbor)) {
        neighborhood.Add(neighbor);
      }
    }
  }
  return neighborhood;
}

TIntV Master::getIntersection(const TIntV& set1, const TIntV& set2) {
  TIntV result;
  for (int i = 0; i < set1.Len(); ++i) {
    if (set2.IsIn(set1[i])) {
      result.Add(set1[i]);
    }
  }
  return result;
}

void Master::ExpMulVertices(
    TIntV& VCC_i, TIntV& mc_j, double r_ij,
    unordered_set<pair<int, int>, pair_hash>& Inserted_Edge) {
  // 步骤 1: 计算 S_L
  TIntV neighborhood_mc_in_psi = getNeighborhood(mc_j, VCC_i);
  TIntV intersection_psi_neigh = getIntersection(VCC_i, neighborhood_mc_in_psi);
  TIntV S_L = get_difference_set(VCC_i, intersection_psi_neigh);

  // 步骤 2: 计算 S_R
  // TIntV neighborhood_psi_in_mc = getNeighborhood(VCC_i, mc_j);
  // TIntV intersection_mc_neigh = getIntersection(mc_j,
  // neighborhood_psi_in_mc); TIntV S_R = get_difference_set(mc_j,
  // intersection_mc_neigh);
  TIntV S_R = mc_j;
  // 步骤 3: 计算初始代价 cost
  int cost = (2 * VCC_i.Len() * mc_j.Len() + (mc_j.Len() * mc_j.Len())) / r_ij;

  std::cout << "ExpMulVertices: S_L size: " << S_L.Len()
            << ", S_R size: " << S_R.Len() << ", cost: " << cost << std::endl;
  // 步骤 4: 循环处理
  while (cost > 0) {
    if (S_L.Empty() || S_R.Empty()) {
      break;
    }

    // 在 S_L 中找到使 d_psi_k(v) 最小的顶点 v
    int v = -1;
    int min_degree_psi = std::numeric_limits<int>::max();
    for (int i = 0; i < S_L.Len(); ++i) {
      int degree = get_degree_in_set(S_L[i], VCC_i);
      if (degree < min_degree_psi) {
        min_degree_psi = degree;
        v = S_L[i];
      }
    }

    // 在 S_R 中找到使 d_mc_j(u) 最小的顶点 u
    int u = -1;
    int min_degree = std::numeric_limits<int>::max();
    for (int i = 0; i < S_R.Len(); ++i) {
      int degree = get_degree_in_set(S_R[i], mc_j);
      if (degree < min_degree) {
        min_degree = degree;
        u = S_R[i];
      }
    }

    // 将边 (v, u) 添加到边集合 Inserted_Edge 中
    Inserted_Edge.insert({v, u});
    std::cout << "ExpMulVertices: Inserted edge (" << v << ", " << u << ")"
              << std::endl;

    // 更新邻域信息 Neighbours(u, Q)（这里简单跳过，可根据实际需求实现）
    // update_neighbour
    // 将 cost 减 1
    cost--;

    // 从 S_L 中移除顶点 v，从 S_R 中移除顶点 u
    S_L.DelIfIn(v);
    S_R.DelIfIn(u);
  }

  // 步骤 5: 返回插入的边集合 Inserted_Edge
  // return Inserted_Edge;
}

// // 辅助函数：计算图中所有点的k顶点连通能力值
// std::vector<int> Master::calculateKVertexConnectivity(const TIntVIntV&
// kvcc_array,
//                                               int num_vertices) {
//   std::vector<int> k_connectivity(num_vertices, 0);
//   for (const auto& kvcc : kvcc_array) {
//     int size = kvcc.size();
//     for (int vertex : kvcc) {
//       k_connectivity[vertex] = std::max(k_connectivity[vertex], size);
//     }
//   }
//   return k_connectivity;
// }

// // 辅助函数：计算插入边后的收益
// double Master::calculateGain(
//     const TIntVIntV& original_kvcc_array,
//     const std::unordered_set<std::pair<int, int>, pair_hash>& inserted_edges,
//     const std::vector<int>& original_connectivity) {
//   TIntVIntV new_kvcc_array = original_kvcc_array;
//   for (const auto& edge : inserted_edges) {
//     if(!G->IsEdge(edge.first, edge.second)) G->AddEdge(edge.first,
//     edge.second);
//   }
//   int num_vertices = original_connectivity.size();
//   std::vector<int> new_connectivity =
//       calculateKVertexConnectivity(new_kvcc_array, num_vertices);
//   double gain = 0;
//   for (int i = 0; i < num_vertices; ++i) {
//     gain += new_connectivity[i] - original_connectivity[i];
//   }
//   for (const auto& edge : inserted_edges) {
//     if(G->IsEdge(edge.first, edge.second)) G->DelEdge(edge.first,
//     edge.second);
//   }
//   return gain;
// }

// // 辅助函数：生成所有可能的b条边的组合
// std::vector<std::unordered_set<std::pair<int, int>, pair_hash>>
// Master::generateEdgeCombinations(
//     const std::unordered_set<std::pair<int, int>, pair_hash>&
//     candidate_edges, int b) {
//   std::vector<std::unordered_set<std::pair<int, int>, pair_hash>>
//   combinations; std::vector<bool> mask(candidate_edges.size());
//   std::fill(mask.begin(), mask.begin() + b, true);

//   std::vector<std::pair<int, int>> edge_list(candidate_edges.begin(),
//                                              candidate_edges.end());
//   do {
//     std::unordered_set<std::pair<int, int>, pair_hash> combination;
//     for (int i = 0; i < candidate_edges.size(); ++i) {
//       if (mask[i]) {
//         combination.insert(edge_list[i]);
//       }
//     }
//     combinations.push_back(combination);
//   } while (std::prev_permutation(mask.begin(), mask.end()));

//   return combinations;
// }