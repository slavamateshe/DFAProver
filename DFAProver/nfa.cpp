#include <iostream>
#include <fstream>
#include "nfa.h"

using namespace std;

nfa* nfa_init(int dim, int n, node* start, node* end) {
	nfa* NFA = (nfa*)malloc(sizeof(nfa));
	NFA->dim = dim;
	NFA->n = n;
	NFA->g = graph_init(n, dim);
	NFA->start = start;
	NFA->end = end;
	return NFA;
}

void nfa_free(nfa* NFA) {
	graph_free(NFA->g, NFA->dim);
	free(NFA->start);
	free(NFA->end);
}

void nfa_add(nfa* NFA, int q, int symb, int q_new) {
	graph_add_arc(NFA->g, q, symb, q_new);
}

void nfa_remove(nfa* NFA, int q, int symb, int q_del) {
	graph_del_arc(NFA->g, q, symb, q_del);
}

nfa* nfa_read(const char* s) {
	int dim, n, start_len, end_len, q, symb, q_new, edges_len;

	fstream f(s);
	f >> dim;
	f >> n;

	f >> start_len;
	int x;

	node* start = NULL;
	for (int i = 0; i < start_len; i++) {
		f >> x;
		start = list_add(start, node_get(x));
	}

	f >> end_len;
	node* end = NULL;
	for (int i = 0; i < start_len; i++) {
		f >> x;
		end = list_add(end, node_get(x));
	}

	nfa* NFA = nfa_init(dim, n, start, end);
	f >> edges_len;
	for (int i = 0; i < edges_len; i++) {
		f >> q;
		f >> symb;
		f >> q_new;
		nfa_add(NFA, q, symb, q_new);
	}
	f.close();
	return NFA;
}

void nfa_to_dot(nfa* NFA, const char* s) {
	ofstream f(s);

	f << "digraph{" << endl;

	for (int i = 0; i < NFA->n; i++) {
		for (int j = 0; j < (1 << NFA->dim); j++) {
			node* n = NFA->g->adj_list[i].symbols[j].head;
			if (n) {
				for (; n; n = n->next)
					f << i << " -> " << n->q << " [label=" << j << "]" << endl;
			}
		}
	}

	f << "}";
	f.close();
}

int nfa_check(nfa* NFA, int str) {
	node* current = NFA->start;
	for (; str; str >>= 1) {
		node* n = current;
		node* qnew = NULL;
		while (n) {
			node* q = NFA->g->adj_list[n->q].symbols[str & 1].head;
			while (q) {
				qnew = list_add(qnew, q);
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