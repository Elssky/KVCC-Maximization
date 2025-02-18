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

(/usr/bin/time -v nohup  ./main -d /home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/dataset/useful/CA-CondMat -c /home/lhy/Snap-For-KVCC/examples/testgraph/community/CA-CondMat -a approximate -k 7 -b 50) >& ./CA-conadmat_k=7_b=50_alg=approximate.txt &