#include <iostream>
#include <fstream>
#include "cli.h"

using namespace std;


void test() {
	nfa* a = nfa_read("2div.txt");
	nfa* b = nfa_read("3div.txt");
	nfa* c = nfa_intersect(a, b);
	nfa* d = nfa_complement2(c);
	//nfa_to_dot(b, "div3.dot");
	for (int i = 0; i < 100; i++) {
		int t[1] = { i };
		cout << i << ": ";
		cout << nfa_check(a, t) << " ";
		cout << nfa_check(b, t) << " ";
		cout << nfa_check(c, t) << " ";
		cout << nfa_check(d, t) << " " << endl;
	}
}

int main()
{  
	nfa* a = nfa_linear_equals(1);
	//nfa* b = nfa_to_dfa(a);
	//nfa* c = nfa_minimize(b);
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {
			int input[2] = { i, j };
			cout << i << " " << j << ": ";
			if (nfa_check(a, input)) {
				cout << "true " << endl;
			}
			else {
				cout << "false " << endl;
			}
			/*if (nfa_check(b, input)) {
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
			}*/
		}
	}
	///cli();
}