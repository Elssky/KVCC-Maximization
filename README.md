# KVCC-Maximization
## Usage
```shell
# for exact real-time update
./main -d /home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/dataset/useful/web-arabic-2005 -c /home/lhy/Snap-For-KVCC/examples/testgraph/community/web-arabic-2005 -k 5 -b 50 -a exact

# for apprroximate real-time update
./main -d /home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/dataset/useful/web-arabic-2005 -c /home/lhy/Snap-For-KVCC/examples/testgraph/community/web-arabic-2005 -k 5 -b 50 -a apprroximate

```

Use `#define _DEBUG` to activate `debug_print` function

(/usr/bin/time -v nohup  ./main -d /home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/dataset/useful/web-arabic-2005 -c /home/lhy/Snap-For-KVCC/examples/testgraph/community/web-arabic-2005 -k 5 -b 50) >& ./web-arabic-2005_k=5.txt &

(/usr/bin/time -v nohup  ./main -d /home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/dataset/useful/CA-CondMat -c /home/lhy/Snap-For-KVCC/examples/testgraph/community/CA-CondMat -a exact -k 7 -b 50) >& ./CA-conadmat_k=7_b=50_alg=exact.txt &

(/usr/bin/time -v nohup  ./main -d /home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/dataset/useful/web-uk-2005 -c /home/lhy/Snap-For-KVCC/examples/testgraph/community/web-uk-2005 -a approximate -k 100 -b 100) >& ./output_0219/web-uk-2005/web-uk-2005_b=100_k=100_alg=approximate.txt &

(/usr/bin/time -v nohup  ./main -d /home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/dataset/useful/web-uk-2005 -c /home/lhy/Snap-For-KVCC/examples/testgraph/community/web-uk-2005 -a exact -k 100 -b 100) >& ./output_0219/web-uk-2005/web-uk-2005_b=100_k=100_alg=exact.txt &

(/usr/bin/time -v nohup  ./main -d /home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/dataset/useful/CA-CondMat -c /home/lhy/Snap-For-KVCC/examples/testgraph/community/CA-CondMat -a approximate -k 11 -b 100 -th 0.9 -u all) >& ./CA-conadmat_k=11_b=100_alg=approximate_th=0.9_unuse=all.txt &

(/usr/bin/time -v nohup  ./main -d /home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/dataset/useful/CA-CondMat -c /home/lhy/Snap-For-KVCC/examples/testgraph/community/CA-CondMat -a approximate -k 11 -b 100 -th 0.9 -u cck) >& ./CA-conadmat_k=11_b=100_alg=approximate_th=0.9_unuse=cck.txt &

(/usr/bin/time -v nohup  ./main -d /home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/dataset/useful/CA-CondMat -c /home/lhy/Snap-For-KVCC/examples/testgraph/community/CA-CondMat -a approximate -k 11 -b 100 -th 0.9 -t 8) >& ./CA-conadmat_k=11_b=100_alg=approximate_th=0.9_t=8.txt &

(/usr/bin/time -v nohup  ./main -d /home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/dataset/useful/road-chesapeake -c /home/lhy/Snap-For-KVCC/examples/testgraph/community/road-chesapeake -a exact -k 5 -b 100) >& ./road-chesapeake_k=5_b=100_alg=exact_.txt &

./main -d /home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/dataset/useful/road -c /home/lhy/Snap-For-KVCC/examples/testgraph/community/road -a exact -k 3 -b 20

./kvcm -d /home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/dataset/useful/dolphins -c /home/lhy/Snap-For-KVCC/examples/testgraph/community/dolphins -k 3 -b 20