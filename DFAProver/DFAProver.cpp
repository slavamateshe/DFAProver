#include <iostream>
#include <fstream>
#include "cli.h"
#include <filesystem>
#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
using namespace std;

void test() {
	nfa* a = nfa_linear_equals(7); 
	cout << " " << a->n << "!!!!!!!!!!!!!!" << endl;
	nfa* b = nfa_projection(a, 0);
	nfa* c = nfa_by_word("101");
	for (int i = 0; i < 100; i++) {
		int input[1] = {i};
		if (nfa_check(b, input)) 
			cout << i << " ";
	}
}

int main(void)
{
	test();
	//cli();
}