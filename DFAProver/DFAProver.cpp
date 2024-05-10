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
	test();
}