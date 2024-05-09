#include <iostream>
#include <fstream>
#include "cli.h"

using namespace std;


int main()
{  
	nfa* a = nfa_linear_equals(52);
	cout << a->n;
	cli();
}