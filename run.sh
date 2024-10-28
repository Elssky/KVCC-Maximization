#!/bin/bash
# exp dataset
# dataset="CA-GrQc"       # k = 3
# dataset="ca-HepPh"    # k = 4
# dataset="ca-MathSciNet" # k = 5
# dataset="email-enron-large" # k = 4, 5, 6, 7, 9
# dataset="fb-pages-company"  # k = 
# dataset="sc-shipsec1"     # k = 
# dataset="sc-shipsec5"     # k = 6, 16, 18, 19, 20, 21, 22, 23, 24, 25, 27, 29, 50, 75, 100, 125
# dataset="soc-epinions"    #
# dataset="Stanford" # k = 15
# dataset="web-arabic-2005" # k = 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 17, 20, 30, 40, 60 
# dataset="web-it-2004"     # k = 11, 13, 15, 17, 19, 22, 43
# dataset="web-sk-2005"     # k = 4, 5, 6, 8, 10, 16, 24, 32, 48, 49, 54, 
# dataset="web-uk-2005"     # k = 100, 150, 200, 250, 299, 300
dataset="web-webbase-2001-all"    # k = 1100, 1200, 1300, 1400,

alg_list=("m" "t" "s")
b_list=(10)
k_list=(1100 1200 1300)

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
            (/usr/bin/time -v nohup ./main -d /home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/dataset/useful/${dataset} -c /home/lhy/Snap-For-KVCC/examples/testgraph/community/$dataset -k $k -b $b -a $alg) >&  ./output/$dataset/${dataset}_followers_b=${b}_k=${k}_alg=${alg}.txt
        done
    done
done