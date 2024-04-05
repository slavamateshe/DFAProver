#include <iostream>
#include <fstream>
#include "nfa.h"

using namespace std;


int main()
{
	nfa* x_pl_y = nfa_read("sum.txt");
	nfa* x_eq_y = nfa_read("equals.txt");
	nfa* z = nfa_extend(x_eq_y, 0);
	
	nfa* a = nfa_intersect(x_pl_y, z);
	nfa* s = nfa_projection(a, 2);
	nfa* p = nfa_projection(s, 2);
	nfa_write(p, "2divc.txt");
	for (int k = 0; k < 100; k++) {
		int* b = (int*)malloc(sizeof(int));
		b[0] = k;
		if (nfa_check(p, b)) {
			cout << k << endl;
		}
	}
	
}