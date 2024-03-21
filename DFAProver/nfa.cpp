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
	for (int i = 0; i < end_len; i++) {
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

nfa* nfa_cartesian(nfa* n1, nfa* n2) {
	node* new_start = NULL;
	int c1, c2, q1, q2;

	for (node* start1 = n1->start; start1; start1 = start1->next) {
		for (node* start2 = n2->start; start2; start2 = start2->next) {
			c1 = start1->q * n1->n + start2->q;
			new_start = list_add(new_start, node_get(c1));
		}
	}

	nfa* new_n = nfa_init(max(n1->dim, n2->dim), n1->n * n2->n, new_start, NULL);

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

nfa* nfa_intersect(nfa* n1, nfa* n2) {
	nfa* new_n = nfa_cartesian(n1, n2);
	node* new_end = NULL;
	int state_num = 0;

	for (node* end1 = n1->end; end1; end1 = end1->next) {
		for (node* end2 = n2->end; end2; end2 = end2->next) {
			state_num = end1->q * n1->n + end2->q;
			new_end = list_add(new_end, node_get(state_num));
		}
	}

	new_n->end = list_add(new_n->end, new_end);
	return new_n;
}

nfa* nfa_union(nfa* n1, nfa* n2) {
	nfa* new_n = nfa_cartesian(n1, n2);
	node* new_end = NULL;
	int state_num = 0;

	for (int i = 0; i < n1->n; i++) {
		for (node* end2 = n2->end; end2; end2 = end2->next) {
			state_num = i * n1->n + end2->q;
			new_end = list_add(new_end, node_get(state_num));
		}
	}


	for (int i = 0; i < n2->n; i++) {
		for (node* end1 = n1->end; end1; end1 = end1->next) {
			state_num = i * n2->n + end1->q;
			new_end = list_add(new_end, node_get(state_num));
		}
	}

	new_n->end = list_add(new_n->end, new_end);

	return new_n;
}
