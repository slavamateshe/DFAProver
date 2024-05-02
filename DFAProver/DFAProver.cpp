#include <iostream>
#include <fstream>
#include "nfa.h"

using namespace std;

void test_mult_const(nfa* a, int t) {
	int input[2] = { 0, 0 };
	for (int i = 0; i < t; ++i) {
		for (int j = 0; j < t; ++j) {
			input[0] = i, input[1] = j;
			if (nfa_check(a, input))
				printf("(%d, %d) -> accept!!!!!!!!\n", i, j);
			else 
				printf("(%d, %d) -> reject\n", i, j);
		}
	}
}

int main()
{
	int buffsize = 128;
	char* input = (char*)malloc(buffsize * sizeof(char));
	fgets(input, buffsize, stdin);
	stack* rpn = infix_to_rpn(input);
	while (rpn->size != 0) {
		cout << stack_top(rpn) << " ";
		stack_pop(rpn);
	}
}