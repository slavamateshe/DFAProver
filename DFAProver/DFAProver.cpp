#include <iostream>
#include <fstream>
#include "nfa.h"

using namespace std;


int main()
{
	nfa* a = nfa_read("2div.txt");
	nfa* b = nfa_extend(a, 0);

	for (int i = 0; i < 100; i++) {
		for (int j = 0; j < 10; j++) {
			int* str = (int*)malloc(2 * sizeof(int));
			str[0] = i;
			str[1] = j;
			if (nfa_check(b, str)) {
				cout << i << " " << j << endl;
			}
		}
	}

}