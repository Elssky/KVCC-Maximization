#include <stdio.h>
#include <iostream>
#include "getopt.h"
#include "Master.h"
using namespace std;

int main(int argc, char* argv[]) {
	int k, threads = 1;
	string alg, seed, mergeMethod, expandMethod;
	TStr dataset;
    int o;
	const char* optstring = "d:a:k:s:m:e:t:"; 
	while ((o = getopt(argc, argv, optstring)) != -1) {
		switch (o) {
		case 'd':
			dataset = optarg;
			printf("opt is d, oprarg is: %s\n", optarg);
			break;
		case 'a':
			alg = optarg;
			printf("opt is a, oprarg is: %s\n", optarg);
			break;
		case 'k':
			k = atoi(optarg);
			printf("opt is k, oprarg is: %s\n", optarg);
			break;
		case 's':
			seed = optarg;
			printf("opt is s, oprarg is: %s\n", optarg);
			break;
		case 'm':
			mergeMethod = optarg;
			printf("opt is m, oprarg is: %s\n", optarg);
			break;
		case 'e':
			expandMethod = optarg;
			printf("opt is e, oprarg is: %s\n", optarg);
			break;
		case 't':
			threads = atoi(optarg);
			printf("opt is t, oprarg is: %s\n", optarg);
			break;
		case '?':
			printf("Correct Usage:\n");
			return 0;
		}
	}

}