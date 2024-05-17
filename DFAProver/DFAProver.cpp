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
	nfa* a = nfa_linear_equals(4); 
	cout << " " << a->n << "!!!!!!!!!!!!!!" << endl;
	for (int i = 0; i < 20; i++) {
		for (int j = 0; j < 20; j++) {
			int input[2] = { i, j };
			cout << i << " " << j << ": ";  
			if (nfa_check(a, input)) {
				cout << "true " << endl;
			}
			else {
				cout << "false" << endl;
			}
		}
	}
}

int main(void)
{
	cli();
}