#include <iostream>
#include <fstream>
#include "nfa.h"

using namespace std;


/*nfa* intersect_nfa(nfa* n1, nfa* n2) {
	node* new_start = NULL;
	node* new_end = NULL;
	int c1, c2, q1, q2;


}*/

int main()
{
	nfa* a = nfa_read("3div.txt");
	//nfa_to_dot(a, "3div.dot");
	//int x;
	//cin >> x;
	for(int i=0; i < 100; i++)
	 if (nfa_check(a, i)) cout << i << endl;
}