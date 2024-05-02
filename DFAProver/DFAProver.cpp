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
	nfa* a = nfa_read("example.txt");
	nfa* b = nfa_minimize(a);
	for (int i = 0; i < 1000; i++) {
		int input[1] = { i };
		cout << i << " " << nfa_check(a, input) << " " << nfa_check(b, input) << endl;
	}
}