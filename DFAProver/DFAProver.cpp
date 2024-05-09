#include <iostream>
#include <fstream>
#include "cli.h"

using namespace std;


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
}

int main()
{  
	cli();
}