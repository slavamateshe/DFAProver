#include <iostream>
#include <fstream>
#include "nfa.h"

using namespace std;


int main()
{
	nfa* b = nfa_read("3div.txt");
	cout << nfa_is_dfa(b);
}