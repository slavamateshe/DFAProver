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
	for (int i = 0; i < 100; i++) {
		int input[1] = {i};
		if (nfa_check(b, input)) 
			cout << i << " ";
	}
}

void test_for(nfa* b, int n) {
	for (int i = 0; i < n; i++) {
		int input[1] = { i };
		if (nfa_check(b, input))
			cout << i << " ";
	}
}

void test2()
{
	const char* s = "( (01) | (10))*";//"( (1101) | (11111) | ((101) | (10))? )";
	nfa* a = nfa_from_regex(s);
	nfa* b = nfa_cut_leading_zeros(a);
	test_for(a, 100);
	cout << endl;
	test_for(b, 100);
	nfa_to_dot(a, "test.dot");
}

int main(void)
{
	test2();
	//cli();
}