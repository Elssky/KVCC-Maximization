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
dataset="web-arabic-2005" # k = 4 5 6 7 8 9 10 11 12 13 14 17 20 30 40 60 
# dataset="web-it-2004"     # k = 11 13 15 17 19 22 43
# dataset="web-sk-2005"     # k = 4 5 6 8 10 16 24 32 48 49 54 
# dataset="web-uk-2005"     # k = 100 150 200 250 299 300
# dataset="web-webbase-2001-all"    # k = 1100 1200 1300 1400

alg_list=("s" "m" "t")
# alg_list=("s" "m" "t" "ma" "mo")
b_list=(50 100)
k_list=(4 5 6 7 8 9 10)

output_dir="./output_1123_new/$dataset"
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

for alg in ${alg_list[@]}
do
    for k in ${k_list[@]}
    do
        for b in ${b_list[@]}
        do
            (/usr/bin/time -v nohup ./main -d /home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/dataset/useful/${dataset} -c /home/lhy/Snap-For-KVCC/examples/testgraph/community/$dataset -k $k -b $b -a $alg) >&  ./output_1123_new/$dataset/${dataset}_followers_b=${b}_k=${k}_alg=${alg}.txt &
        done
    done
done