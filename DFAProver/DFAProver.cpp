#include <iostream>
#include <fstream>
#include "cli.h"

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
	nfa* a = nfa_linear_equals(3);
	nfa* b = nfa_to_dfa(a);
	int k = 0, m = 0;
	for (int i = 0; i < 1000; i++) {
		int c[1] = { i };
		k = nfa_check(a, c);
		m = nfa_check(b, c);
		if (k != m) cout << -1;
	}

}
