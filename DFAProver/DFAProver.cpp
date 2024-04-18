#include <iostream>
#include <fstream>
#include "nfa.h"

using namespace std;


int main()
{
	nfa* a = nfa_read("а1.txt");
	nfa* b = nfa_read("a2_2.txt");
	nfa* c = right_quot(a, b);
	nfa_to_dot(c, "sum.ptoj");
}