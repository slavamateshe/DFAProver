#include <iostream>
#include <fstream>
#include "nfa.h"

using namespace std;


int main()
{
	nfa* a = nfa_read("2div.txt");
	nfa* a_c = nfa_complement(a);
	nfa* all = nfa_union(a, a_c);
	

	for (int i = 0; i < 100; i++) {
		int input[1] = { i };
		if (nfa_check(all, input)) {
			cout << i << endl;
		}
	}
	//node* n0 = node_get(0);
	//nfa* sum = nfa_init(3, 2, n0, n0);
	//nfa_add(sum, 0, 0, 0);
	//nfa_add(sum, 0, 3, 0);
	//nfa_add(sum, 0, 5, 0);

	//nfa_add(sum, 0, 6, 1);
	//nfa_add(sum, 1, 1, 0);

	//nfa_add(sum, 1, 2, 1);
	//nfa_add(sum, 1, 4, 1);
	//nfa_add(sum, 1, 7, 1);
	//int input[3] = { 14,10,24 }; // 14 == 0111 ; 10 == 0101; 25 == 00011
	//for (int i = 0; i < 30; i++) {
	//	input[2] = i;
	//	cout << nfa_check(sum, input);
}