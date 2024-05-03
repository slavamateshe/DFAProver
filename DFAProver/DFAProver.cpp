#include <iostream>
#include <fstream>
#include "cli.h"

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

void test_rpn() {
	/// (($div2(x) & $div3(x)) | ~$div2(x))
	
	int buffsize = 128;
	char* input = (char*)malloc(buffsize * sizeof(char));
	fgets(input, buffsize, stdin);
	stack* rpn = infix_to_rpn(input);
	nfa* a = rpn_to_nfa(rpn);
	if (!a) {
		return;
	}
	for (int i = 0; i < 20; i++) {
		int input[1] = { i };
		if (nfa_check(a, input)) {
			cout << i << " accept!!!!" << endl;
		}
		else {
			cout << i << " reject" << endl;
		}
	}
	cout << a->n << endl;

	nfa* b = nfa_minimize(a);
	for (int i = 0; i < 20; i++) {
		int input[1] = { i };
		if (nfa_check(b, input)) {
			cout << i << " accept!!!!" << endl;
		}
		else {
			cout << i << " reject" << endl;
		}
	}
	cout << b->n;
}

//void text_prover() {
//	int buffsize = 128;
//	char* input = (char*)malloc(buffsize * sizeof(char));
//	fgets(input, buffsize, stdin);
//	if (substr())
//
//}

int main()
{  
	test_rpn();
	
}