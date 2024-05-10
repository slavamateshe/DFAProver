#include <iostream>
#include <fstream>
#include "cli.h"

using namespace std;


int main()
{  
<<<<<<< HEAD
	nfa* a = nfa_read("2div.txt");
	nfa* b = nfa_minimize(nfa_to_dfa(a));
	for (int i = 0; i < 100; i++) {
		int bb[1] = { i };
		if (!(nfa_check(b, bb))) cout << i << endl;
	}
=======
	nfa* a = nfa_linear_equals(52);
	cout << a->n;
	cli();
>>>>>>> 745e94ccd2f072c56cc96fd08bf0c167a311eebd
}