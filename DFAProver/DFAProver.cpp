#include <iostream>
#include <fstream>
#include "nfa.h"
using namespace std;


// Ez( z = x + y) -> all pairs (x,y) in N^2 [0,0][1,1] 
//
//

void test2(nfa* a, int n) {
	int count = a->dim;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			int input[2] = { i, j };
			printf("(%d,%d) -> %d\n",i,j,nfa_check(a, input));
		}
	}
}

void test2(nfa* a, nfa* b, int n) {
	int count = a->dim;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			int input[2] = { i, j };
			printf("(%d,%d) -> a1: %d a2: %d\n", i, j, nfa_check(a, input), nfa_check(b, input));
		}
	}
}

int main()
{
	nfa* a = nfa_init(1, 1, node_get(0), node_get(0));
	nfa_add(a, 0, 0, 0);

	nfa* b = nfa_init(1, 1, node_get(1), node_get(1));
	nfa_add(a, 1, 1, 1);

	nfa* c = nfa_union(a, b);
	nfa_to_dot(c, "file1.dot");
}
