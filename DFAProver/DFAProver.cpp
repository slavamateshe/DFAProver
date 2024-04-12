#include <iostream>
#include <fstream>
#include "nfa.h"

using namespace std;


int main()
{
	nfa* a = nfa_read("equals.txt");
	nfa* b = nfa_read("sum.txt");
	nfa* c = nfa_intersect(a, b);
	nfa* d = nfa_intersect(c, b);
	nfa* e = nfa_projection(d, 1);
	nfa* f = nfa_projection(e, 1);
	nfa* g = nfa_projection(f, 1);
	int* r = (int*)malloc(sizeof(int));
	for (int i = 0; i < 10; i++) {
		r[0] = i;
		if (nfa_check(g, r)) cout << i << endl;
	}
}