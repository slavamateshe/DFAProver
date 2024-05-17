#include <iostream>
#include <fstream>
#include "cli.h"

using namespace std;


void test() {
	nfa* a = nfa_linear_equals(2); 
	cout << " " << a->n << "!!!!!!!!!!!!!!" << endl;
	nfa* b = nfa_to_dfa(a);
	nfa* b = nfa_to_dfa(a);
	for (int i = 0; i < 20; i++) {
		for (int j = 0; j < 20; j++) {
			int input[2] = { i, j };
			cout << i << " " << j << ": ";  
			if (nfa_check(a, input)) {
				cout << "true " << endl;
			}
		}
	}
}

int main()
{  
	test();
}