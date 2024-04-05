#include <iostream>
#include <fstream>
#include "nfa.h"

using namespace std;


int main()
{
	nfa* a = nfa_read("sum.txt");
	nfa* b = nfa_swap(a, 0, 1);
	nfa_to_dot(b, "sum_to_dot.dot");
}