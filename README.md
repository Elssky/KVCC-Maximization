# KVCC-Maximization
## Usage
```shell
# for single vertex maximization
./main -d /home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/dataset/useful/web-arabic-2005 -c /home/lhy/Snap-For-KVCC/examples/testgraph/community/web-arabic-2005 -k 5 -b 50 -a s

# for multi vertex maximization
./main -d CA-GrQc -k 5 -b 10 -a m

# for together maximization
./main -d CA-GrQc -k 5 -b 10 -a t
```

Use `#define _DEBUG` to activate `debug_print` function