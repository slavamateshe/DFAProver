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
	nfa* a = nfa_read("sum.txt"); // x + y = z
	nfa* b = nfa_projection(a,2);
	node* start = node_get(0);
	nfa* zero = nfa_init(2, 1, start, start);
	nfa_add(zero, 0, 0, 0);
	nfa* correct_add = left_quot(b, zero);
	//test2(correct_add, 10);

	node* start_odd = node_get(0);
	node* final_odd = node_get(1);
	nfa* odd = nfa_init(2, 2, start_odd, final_odd);
	nfa_add(odd, 0, 3, 1);
	for(int i = 0; i < 4; i++) nfa_add(odd, 1, i, 1);
	//test2(odd, 10); 
	nfa* test_odd = left_quot(odd, odd); //Ex(x in odd and y = xz)
	test2(odd, test_odd, 10);
	//nfa_to_dot(correct_add, "corr_add.dot");
	//nfa_to_dot(odd, "odd.dot");
	//nfa_to_dot(test_odd, "test_odd.dot");

}
