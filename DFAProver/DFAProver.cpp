#include <iostream>
#include <fstream>
#include "nfa.h"

using namespace std;


int main()
{
	nfa* t = nfa_linear_equals(4);

	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {
			int b[5] = { i, j,};
			if (nfa_check(t, b)) {
				cout << i << " " << " " << j << endl;
			}
		}
	}
}