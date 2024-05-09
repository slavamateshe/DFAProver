#include <iostream>
#include <fstream>
#include "cli.h"

using namespace std;


int main()
{  
	nfa* a = nfa_read("2div.txt");
	nfa* b = nfa_minimize(nfa_to_dfa(a));
	for (int i = 0; i < 100; i++) {
		int bb[1] = { i };
		if (!(nfa_check(b, bb))) cout << i << endl;
	}
}