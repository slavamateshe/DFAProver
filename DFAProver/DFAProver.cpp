#include <iostream>
#include <fstream>
#include "nfa.h"

using namespace std;


int main()
{
	nfa* a = nfa_read("2div.txt");
	nfa* b = nfa_read("3div.txt");

	nfa* c = nfa_union(a, b);
	for (node* n = c->start; n; n = n->next) {
		cout << n->q << " ";
	}
	cout << endl;
	for (node* n = c->end; n; n = n->next) {
		cout << n->q << " ";
	}
	cout << endl;
	
	for (int i = 0; i < 100; i++) {
		if (nfa_check(c, i)) cout << i << endl;
	}
}