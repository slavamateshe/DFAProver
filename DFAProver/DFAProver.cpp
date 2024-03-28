
#include <iostream>
#include <fstream>
#include "nfa.h"

using namespace std;


int main()
{
    nfa* a = nfa_read("2div.txt");
	nfa* b = nfa_read("3div.txt");


	//nfa_add(a, 0, 1, 0);
	nfa* c = nfa_union(a, b);
	nfa_to_dot(c, "2_2div.txt");
	for (int i = 0; i < 100; i++) {
		if (nfa_check(c, i)) cout << i << endl;
	}
}