#include <iostream>
#include <fstream>
#include "nfa.h"

using namespace std;


int nfa_check(nfa* NFA, int str) {
	node* current = NFA->start;
	for (; str; str >>= 1) {
		node* n = current;
		node* qnew = NULL;
		while (n) {
			node* q = NFA->g->adj_list[n->q].symbols[str & 1].head;
			while (q) {
				qnew = add_to_list(qnew, q);
				q = q->next;
			}
			n = n->next;
		}
		current = qnew;
	}

	while (current) {
		node* n = NFA->end;
		while (n) {
			if (current->q == n->q) {
				return 1;
			}
			n = n->next;
		}
		current = current->next;
	}
	return 0;
}

nfa* union_nfa(nfa* n1, nfa * n2) {
	node* new_start = NULL;
	node* new_end = NULL;
	int c1, c2, q1, q2;

	for (node* start1 = n1->start; start1; start1 = start1->next) {
		for (node* start2 = n2->start; start2; start2 = start2->next) {
			c1 = start1->q * n1->n + start2->q;
			new_start = add_to_list(new_start, node_get(c1));
		}
	}

	for (node* end1 = n1->end; end1; end1 = end1->next) {
		for (node* end2 = n2->end; end2; end2 = end2->next) {
			c2 = end1->q * n1->n + end2->q;
			new_end = add_to_list(new_end, node_get(c2));
		}
	}

	nfa* new_n = nfa_init(max(n1->dim, n2->dim), n1->n * n2->n, new_start, new_end);

	for (int symb = 0; symb < pow(2, min(n1->dim, n2->dim)); symb++) {
		for (int i = 0; i < n1->n; i++) {
			for (int j = 0; j < n2->n; j++) {
				q1 = i * n1->n + j;
				for (node* nd1 = n1->g->adj_list[i].symbols[symb].head; nd1; nd1 = nd1->next) {
					for (node* nd2 = n2->g->adj_list[j].symbols[symb].head; nd2; nd2 = nd2->next) {
						q2 = nd1->q * n1->n + nd2->q;
						nfa_add(new_n, q1, symb, q2);
					}
				}

			}

		}
	}
	return new_n;
}

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