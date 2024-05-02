#include <iostream>
#include <fstream>
#include "nfa.h"

using namespace std;

void test_mult_const(nfa* a, int t) {
	int input[2] = { 0, 0 };
	for (int i = 0; i < t; ++i) {
		for (int j = 0; j < t; ++j) {
			input[0] = i, input[1] = j;
			if (nfa_check(a, input))
				printf("(%d, %d) -> accept!!!!!!!!\n", i, j);
			else 
				printf("(%d, %d) -> reject\n", i, j);
		}
	}
}

int main()
{
	nfa* t = nfa_linear_equals(5); // y = 5*x (15,3) i.e. (5*x,x)
	cout << t->n;
	nfa* b = nfa_del_unrechable(t);
	cout << b->n << endl;
	test_mult_const(t, 6);

}