#include <iostream>
#include <fstream>
#include "nfa.h"

using namespace std;

nfa* nfa_init(int dim, int n, node* start, node* end) {
	nfa* NFA = (nfa*)malloc(sizeof(nfa));
	if (NFA) {
		NFA->dim = dim;
		NFA->n = n;
		NFA->g = graph_init(n, dim);
		NFA->start = start;
		NFA->end = end;
	}
	return NFA;
}

void nfa_free(nfa* NFA) {
	graph_free(NFA->g, NFA->dim);
	list_free(NFA->start);
	list_free(NFA->end);
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
	if (!f.is_open()) {
		cout << "file not exist" << endl;
		return NULL;
	}

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

void nfa_write(nfa* NFA, const char* s) {
	int k1 = 0, k2 = 0;
	ofstream f(s);

	f << NFA->dim << endl;
	f << NFA->n << endl;
	for (node* curr = NFA->start; curr; curr = curr->next, k1++);
	f << k1 << endl; 
	for (node* curr = NFA->start; curr; curr = curr->next)
		f << curr->q << endl;
	for (node* curr = NFA->end; curr; curr = curr->next, k2++);
	f << k2 << endl;
	for (node* curr = NFA->end; curr;  curr = curr->next)
		f << curr->q << endl;

	for (int i = 0; i < NFA->n; i++) {
		for (int j = 0; j < (1 << NFA->dim); j++) {
			node* n = NFA->g->adj_list[i].symbols[j].head;
			if (n) {
				for (; n; n = n->next)
					f << i << " " << j << " " << n->q << endl;
			}
		}
	}
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

int nfa_check(nfa* NFA, int* str) {
	node* current = NFA->start;
	for (int k = 0; ; k++) {
		int symb = 0;
		for (int i = 0; i < NFA->dim; i++) {
			symb += ((str[i] >> k) & 1) << i;
		}
		if (symb == 0) {
			break;
		}

		node* n = current;
		node* qnew = NULL;
		while (n) {
			node* q = NFA->g->adj_list[n->q].symbols[symb].head;
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
			c1 = start2->q * n1->n + start1->q;
			new_start = list_add(new_start, node_get(c1));
		}
	}

	nfa* new_n = nfa_init(max(n1->dim, n2->dim), n1->n * n2->n, new_start, NULL);

	for (int symb = 0; symb < pow(2, min(n1->dim, n2->dim)); symb++) {
		for (int i = 0; i < n1->n; i++) {
			for (int j = 0; j < n2->n; j++) {
				q1 = j * n1->n + i;
				for (node* nd1 = n1->g->adj_list[i].symbols[symb].head; nd1; nd1 = nd1->next) {
					for (node* nd2 = n2->g->adj_list[j].symbols[symb].head; nd2; nd2 = nd2->next) {
						q2 = nd2->q * n1->n + nd1->q;
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
			state_num = end2->q * n1->n + end1->q;
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

	for (node* end2 = n2->end; end2; end2 = end2->next) {
		for (int i = 0; i < n1->n; i++) {
			state_num = end2->q * n1->n + i;
			new_end = list_add(new_end, node_get(state_num));
		}
	}

	for (node* end1 = n1->end; end1; end1 = end1->next) {
		for (int i = 0; i < n2->n; i++) {
			state_num = i * n1->n + end1->q;
			new_end = list_add(new_end, node_get(state_num));
		}
	}

	new_n->end = new_end;
	return new_n;
}

int nfa_is_dfa(nfa* n) {
	node* nd_1;
	node* nd_2;

	for (int i = 0; i < n->n; i++) {
		for (int j = 0; j < pow(2, n->dim); j++) {

			nd_1 = n->g->adj_list[i].symbols[j].head;
			if (nd_1) {
				nd_2 = nd_1->next;
				if (nd_2 && nd_2->next) return 0;
			}
		}
	}
	return 1;
}

nfa* nfa_projection(nfa* a, int n) {
	int new_symb;
	node* start = NULL;

	for (node* x = a->start; x; x = x->next) {
		start = list_add(start, node_get(x->q));
	}
	node* end = NULL;
	for (node* x = a->end; x; x = x->next) {
		end = list_add(end, node_get(x->q));
	}
	nfa* b = nfa_init(a->dim - 1, a->n, start, end);

	for (int i = 0; i < a->n; i++) {
		for (int symb = 0; symb < (1 << a->dim); symb++) {
			for (node* nd = a->g->adj_list[i].symbols[symb].head; nd; nd = nd->next) {
				new_symb = ((symb >> n) << (n - 1)) + (((1 << n) - 1) & symb);
				nfa_add(b, i, symb, nd->q);
			}
		}
	}
	return b;
}


nfa* nfa_extend(nfa* a, int n) {
	node* start = NULL;
	for (node* x = a->start; x; x = x->next) {
		start = list_add(start, node_get(x->q));
	}
	node* end = NULL;
	for (node* x = a->end; x; x = x->next) {
		end = list_add(end, node_get(x->q));
	}
	nfa* b = nfa_init(a->dim + 1, a->n, start, end);

	for (int i = 0; i < a->n; i++) {
		for (int symb = 0; symb < (1 << a->dim); symb++) {
			for (node* nd = a->g->adj_list[i].symbols[symb].head; nd; nd = nd->next) {
				int right = symb & ((1 << n) - 1);
				int left = ((symb - right) << (n + 1));
				nfa_add(b, i, left + 0 + right, nd->q);
				nfa_add(b, i, left + (1 << n) + right, nd->q);
			}
		}
	}
	return b;
}