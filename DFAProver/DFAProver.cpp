
#include <iostream>
#include <fstream>
#include "nfa.h"

using namespace std;


int main()
{
    nfa* a = nfa_read("2div.txt");

    for (int i = 0; i < 100; i++) {
        if (nfa_check(a, i)) cout << i << endl;
    }
}


