#include <iostream>
#include <fstream>
#include "cli.h"

using namespace std;


void test() {
	nfa* a = nfa_linear_equals(5);
	//nfa* b = nfa_to_dfa(a);
	for (int i = 0; i < 20; i++) {
		for (int j = 0; j < 20; j++) {
			int input[2] = { i, j };
			cout << i << " " << j << ": ";  
			if (nfa_check(a, input)) {
				cout << "true " << endl;
			}
			else {
				cout << "false " << endl;
			}

			/*if (nfa_check(b, input)) {
				cout << "true " << endl;
			}
			else {
				cout << "false " << endl;
			}*/
		}
	}
}

int main()
{
	//test();
	//cli();

	nfa* a = power_of2(2);
	nfa_to_dot(a, "test.dot");
	for (int i = 0; i < 20; ++i)
	{
		for (int j = 0; j < 20; ++j) {
			int t[2] = { i, j };
			cout << i << " " << j << " " << nfa_check(a, t) << endl;
		}
	}
}