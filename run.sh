#!/bin/bash
# exp dataset

# dataset="CA-CondMat" # k = 3 5 7 8 9 10 11 13 15 20 23 25 34
# dataset="CA-GrQc"       # k = 3 5 13 
# dataset="bauru5727" # k = 3 4 5 6
# dataset="bio-pdb1HYS" # k = 30 60

# dataset="cit-patent" # k = 13 15 17 19 26 38 51 64
# dataset="coAuthorsCitesee" # k = 7 9 11 13 15 17 21 26 34 52 69 86
# dataset="DBLP" # k = 7 9 11 13 15 20 23 34 45 46 68 69 86 90 92 192 258 344
# dataset="socfb-konect" # k = 6 7 8 9 10 11 13 15 
# dataset="web-google-dir" # k = 13 15 16 17 18 19 21 23 26 38 51 64


# dataset="ca-HepPh"    # k = 3 4 5 6 7
# dataset="ca-MathSciNet" # k = 5 7 9 10 11 13 14 15 19 20 24
# dataset="email-enron-large" # k = 4 5 6 7 9
# dataset="fb-pages-company"  # k = 
# dataset="sc-shipsec1"     # k = 
# dataset="sc-shipsec5" # k = 6 16 18 19 20 21 22 23 24 25 27 29 50 75 100 125
# dataset="soc-epinions"    #
# dataset="Stanford" # k = 15
# dataset="web-arabic-2005" # k = 4 5 6 7 8 9 10 11 12 13 14 17 20 30 40 60 
# dataset="web-it-2004"     # k = 11 13 15 17 19 22 43
# dataset="web-sk-2005"     # k = 4 5 6 8 10 16 24 32 48 49 54 
# dataset="web-uk-2005"     # k = 100 150 200 250 299 300
# dataset="web-webbase-2001-all"    # k = 1100 1200 1300 1400

declare -A dataset_k_map=(
    ["CA-CondMat"]="3 5"
    ["CA-GrQc"]="3 5"
    ["bauru5727"]="3 4"
    ["bio-pdb1HYS"]="30 60"
    ["cit-patent"]="13 15"
    ["coAuthorsCitesee"]="7 9"
    ["DBLP"]="7 9"
    ["socfb-konect"]="6 7"
    ["web-google-dir"]="13 15"
    ["ca-HepPh"]="3 4"
    ["ca-MathSciNet"]="5 7"
    ["email-enron-large"]="4 5"
    ["sc-shipsec5"]="6 16"
    ["Stanford"]="15"
    ["web-arabic-2005"]="4 5"
    ["web-it-2004"]="11 13"
    ["web-sk-2005"]="4 5"
    ["web-uk-2005"]="100 150"
    ["web-webbase-2001-all"]="1100 1200"
)

# declare -A dataset_k_map=(
#     ["CA-CondMat"]="3 5 7 8 9 10 11 13 15 20 23 25 34"
#     ["CA-GrQc"]="3 5 13"
#     ["bauru5727"]="3 4 5 6"
#     ["bio-pdb1HYS"]="30 60"
#     ["cit-patent"]="13 15 17 19 26 38 51 64"
#     ["coAuthorsCitesee"]="7 9 11 13 15 17 21 26 34 52 69 86"
#     ["DBLP"]="7 9 11 13 15 20 23 34 45 46 68 69 86 90 92 192 258 344"
#     ["socfb-konect"]="6 7 8 9 10 11 13 15"
#     ["web-google-dir"]="13 15 16 17 18 19 21 23 26 38 51 64"
#     ["ca-HepPh"]="3 4 5 6 7"
#     ["ca-MathSciNet"]="5 7 9 10 11 13 14 15 19 20 24"
#     ["email-enron-large"]="4 5 6 7 9"
#     ["sc-shipsec5"]="6 16 18 19 20 21 22 23 24 25 27 29 50 75 100 125"
#     ["Stanford"]="15"
#     ["web-arabic-2005"]="4 5 6 7 8 9 10 11 12 13 14 17 20 30 40 60"
#     ["web-it-2004"]="11 13 15 17 19 22 43"
#     ["web-sk-2005"]="4 5 6 8 10 16 24 32 48 49 54"
#     ["web-uk-2005"]="100 150 200 250 299 300"
#     ["web-webbase-2001-all"]="1100 1200 1300 1400"
# )

# 算法列表
alg_list=("exact" "approximate")
b_list=(100)
pids=()
# 遍历所有数据集
for dataset in "${!dataset_k_map[@]}"
do
    # 获取对应的 k 值
    k_list=(${dataset_k_map[$dataset]})

    # 如果 k_list 为空，跳过该数据集
    if [ ${#k_list[@]} -eq 0 ]; then
        echo "Skipping dataset $dataset: no k values specified."
        continue
    fi

    # 输出目录
    output_dir="./output_0219/$dataset"
    if [ ! -d "$output_dir" ]; then
        mkdir -p "$output_dir"
    fi

    echo "Running experiments for dataset: $dataset"
    echo "k_list: ${k_list[@]}"

    # 运行实验
    for alg in ${alg_list[@]}
    do
        for k in ${k_list[@]}
        do
            for b in ${b_list[@]}
            do
                output_file="./output_0219/$dataset/${dataset}_followers_b=${b}_k=${k}_alg=${alg}.txt"
                echo "Running: dataset=$dataset, k=$k, b=$b, alg=$alg"
                (/usr/bin/time -v nohup ./main -d /home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/dataset/useful/${dataset} -c /home/lhy/Snap-For-KVCC/examples/testgraph/community/$dataset -k $k -b $b -a $alg) >& "$output_file" &
                pids+=($!)
                if [ ${#pids[@]} -ge 16 ]; then
                    # 等待第一个完成
                    wait ${pids[0]}
                    # 删除已完成任务的 PID
                    pids=("${pids[@]:1}")
                fi
            done
        done
    done
done

for pid in "${pids[@]}"; do
    wait $pid
done

echo "All experiments have been started in the background."



# alg_list=("exact" "approximate")
# # alg_list=("s" "m" "t" "ma" "mo")
# b_list=(50 100)
# k_list=(4 5 6 7 8 9 10)

# output_dir="./output_0219/$dataset"
# if [ ! -d "$output_dir" ]; then
#     mkdir -p "$output_dir"
# fi

# for alg in ${alg_list[@]}
# do
#     for k in ${k_list[@]}
#     do
#         for b in ${b_list[@]}
#         do
#             (/usr/bin/time -v nohup ./main -d /home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/dataset/useful/${dataset} -c /home/lhy/Snap-For-KVCC/examples/testgraph/community/$dataset -k $k -b $b -a $alg) >&  ./output_1123_new/$dataset/${dataset}_followers_b=${b}_k=${k}_alg=${alg}.txt &
#         done
#     done
# done