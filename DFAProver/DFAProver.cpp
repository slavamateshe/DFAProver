#include <iostream>
#include <fstream>
#include "nfa.h"

using namespace std;


int main()
{
	nfa* a = nfa_read("sum.txt");
	cout << a->g->adj_list[0].symbols[5].head->q << endl;
	nfa* b = nfa_projection(a, 1);
	nfa_write(b, "sum_2.txt");
	/*nfa_add(a, 0, 1, 0);
	nfa_write(a, "3div_prime.txt");
	a = nfa_read("3div_prime.txt");
	nfa_to_dot(a, "3div.dot");*/
}