#include <iostream>
#include <fstream>
#include "nfa.h"

using namespace std;


int main()
{
    nfa* a = nfa_read("2div.txt");
	nfa* b = nfa_read("3div.txt");

	nfa* c = nfa_intersect(a, b);
	cout << sizeof(c->start);

	for (int i = 0; i < 100; i++) {
		if (nfa_check(c, i)) {
			cout << i << endl;
		}
	}
}