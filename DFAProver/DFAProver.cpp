#include <iostream>
#include <fstream>
#include "cli.h"

using namespace std;


void test() {
	nfa* a = nfa_linear_equals(7);
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
	cli();

	/*nfa* a = nfa_read("2div.txt");
	nfa* a1 = nfa_to_dfa(a);
	nfa* b = nfa_complement2(a1);
	nfa* c = nfa_intersect(a1, b);

	for (int i = 0; i < 20; ++i)
	{
		int g[1] = { i };
		cout << i << " " << nfa_check(a1, g) << " " << nfa_check(b, g) << " " << nfa_check(c, g) << endl;
	}*/
}