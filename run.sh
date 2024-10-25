#!/bin/bash
# exp dataset
# dataset="twitter_copen" # n=8580 f:0
# dataset="pkustk02"
# dataset="Email-Enron"
# dataset="Facebook"
# dataset="Gowalla" #n=63731 f:0
# dataset="DBLP" #n=317080 f:0
# dataset="Amazon"
# dataset="soc-youtube"
# dataset="Google" #n=875713 f:0
# dataset="soc-lastfm"
# dataset="pokec" #n=1632803 f:0
# dataset="soc-flixster"
# dataset="LiveJournal" #n=4847571 f:1

# useful dataset
# dataset="Berkstan" # n=685230 f:0
# dataset="BlogCatalog" #n=88784 f:0
# dataset="orkut" #n=3072626 f:0
# dataset="CA-CondMat" 

# dataset=sample_graph
# dataset="scc_reality" # n:6809 f:0
# dataset="fb_messages" #n=1899 f:0
# dataset="infect_dublin" #n=10972 f:0
# dataset="fb-pages-food"
# dataset="socfb-Bingham82"
# dataset="ia-infect-hyper"
# dataset="soc-dolphins"
# dataset="Journals"
# dataset="ck104"
# dataset="polbooks"
# dataset="fb-pages-sport"
# dataset="fb-pages-public-figure"
# dataset="air03"

# pre_list=(0 1)  # 1 pre 0 without pre
# mode_list=(0 1 2)  # 0 vertex 1 group 2 vertex+group
# b_list=(10 50)
alg_list=("m" "t" "s")
b_list=(10)
k_list=(5)

output_dir="./output/$dataset"
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

for b in ${b_list[@]}
do
    for k in ${k_list[@]}
    do
        for alg in ${alg_list[@]}
        do
            (/usr/bin/time -v nohup ./main -a $alg -d /home/public/lxw/datasets/${dataset}/${dataset}_new.txt -k $k -b $b) >&  ./output/$dataset/${dataset}_followers_b=${b}_k=${k}_alg=${alg}.txt
            # echo ./output/$dataset/${dataset}_LhCDScvx_k=${k}_h=${h}_t=${t}.txt
        done
    done
done