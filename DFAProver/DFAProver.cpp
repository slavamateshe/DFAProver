#include <iostream>
#include <fstream>
#include "cli.h"

using namespace std;


int main()
{  
	nfa* a = nfa_linear_equals(2);
	nfa* b = nfa_to_dfa(a);
	nfa* c = nfa_minimize(b);
	for (int i = 0; i < 20; i++) {
		for (int j = 0; j < 20; j++) {
			int input[2] = { i, j };
			cout << i << " " << j << ": ";
			if (nfa_check(a, input)) {
				cout << "true ";
			}
			else {
				cout << "false ";
			}
			if (nfa_check(b, input)) {
				cout << "true ";
			}
			else {
				cout << "false ";
			}
			if (nfa_check(c, input)) {
				cout << "true" << endl;
			}
			else {
				cout << "false" << endl;
			}
		}
	}
	///cli();
}