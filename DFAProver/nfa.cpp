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

nfa* nfa_copy(nfa* a) {
	node* new_start = NULL;
	node* new_end = NULL;
	for (node* x = a->start; x; x = x->next) {
		new_start = list_add(new_start, node_get(x->q));
	}
	for (node* x = a->end; x; x = x->next) {
		new_end = list_add(new_end, node_get(x->q));
	}
	nfa* new_nfa = nfa_init(a->dim, a->n, new_start, new_end);

	for (int i = 0; i < a->n; i++) {
		for (int symb = 0; symb < (1 << a->dim); symb++) {
			for (node* nd = a->g->adj_list[i].symbols[symb].head; nd; nd = nd->next) {
				nfa_add(new_nfa, i, symb, nd->q);
			}
		}
	}
	return new_nfa;
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
	int k1 = 0, k2 = 0, k = 0;
	ofstream f(s);

	f << NFA->dim << endl;
	f << NFA->n << endl;
	for (node* curr = NFA->start; curr; curr = curr->next, k1++);
	f << k1 << endl;
	for (node* curr = NFA->start; curr; curr = curr->next)
		f << curr->q << endl;
	for (node* curr = NFA->end; curr; curr = curr->next, k2++);
	f << k2 << endl;
	for (node* curr = NFA->end; curr; curr = curr->next)
		f << curr->q << endl;

	for (int i = 0; i < NFA->n; i++) {
		for (int j = 0; j < (1 << NFA->dim); j++) {
			node* n = NFA->g->adj_list[i].symbols[j].head;
			if (n) k++;
		}
	}
	f << k << endl;

	for (int i = 0; i < NFA->n; i++) {
		for (int j = 0; j < (1 << NFA->dim); j++) {
			node* n = NFA->g->adj_list[i].symbols[j].head;
			if (n) {
				for (; n; n = n->next)
					f << i << " " << j << " " << n->q << endl;
			}
		}
	}
	f.close();
}

void nfa_to_dot(nfa* NFA, const char* s) {
	ofstream f(s);

	f << "digraph{" << endl;

	f << "node [shape = doublecircle] ";

	for (node* curr = NFA->end; curr; curr = curr->next)
		f << curr->q << " ";

	f << endl;

	f << "node [shape = circle]" << endl;

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

int nfa_DFS(node* n, int* checked, nfa* NFA) {
	checked[n->q] = 1;
	for (node* t = NFA->end; t; t = t->next) {
		if (n->q == t->q) return 1;
	}
	for (node* t = NFA->g->adj_list[n->q].symbols[0].head; t; t = t->next) {
		if (checked[t->q] == 0 && nfa_DFS(t, checked, NFA)) {
			return 1;
		}
	}
	return 0;
}

int nfa_check(nfa* NFA, int* str) {
	//when ||str|| == 0 the function doesnt work correctly
	//we have to find all the reachable states from initial states
	//and then compare them with the final states

	if (NFA->dim == 0) {
		int k = 0;
		int* checked = (int*)calloc(NFA->n, sizeof(int));
		for (node* n = NFA->start; n; n = n->next) {
			if (checked[n->q] == 0) {
				if (nfa_DFS(n, checked, NFA)) {
					return 1;
				}
			}
			checked = (int*)calloc(NFA->n, sizeof(int));
		}
		return 0;
	}

	node* current = NFA->start;
	for (int k = 0; ; k++) {
		int symb = 0;
		bool flag = false;
		for (int i = 0; i < NFA->dim; i++) {
			symb += ((str[i] >> k) & 1) << (NFA->dim - i - 1);
			if ((str[i] >> k) > 0) {
				flag = true;
			}
		}
		if (!flag) {
			break;
		}

		node* qnew = NULL;
		for (node* n = current; n; n = n->next) {
			node* q = NFA->g->adj_list[n->q].symbols[symb].head;
			while (q) {
				qnew = list_add(qnew, node_get(q->q));
				q = q->next;
			}
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

nfa* nfa_del_unrechable(nfa* a) {
	node* reachable = node_get(a->start->q);
	node* new_states = node_get(a->start->q);
	while (new_states != NULL) {
		node* temp = NULL;
		for (node* n = new_states; n != NULL; n = n->next) {
			for (int i = 0; i < (1 << (a->dim)); i++) {
				for (node* t = a->g->adj_list[n->q].symbols[i].head; t != NULL; t = t->next) {
					temp = list_add(temp, node_get(t->q));
				}
			}
		}
		new_states = NULL;
		for (node* n = temp; n != NULL; n = n->next) {
			if (!node_in_list(n, reachable)) {
				new_states = list_add(new_states, node_get(n->q));
			}
		}
		for (node* n = new_states; n != NULL; n = n->next) {
			reachable = list_add(reachable, node_get(n->q));
		}
	}
	int* new_q = (int*)malloc(a->n * sizeof(int));
	int k = 0;
	for (int i = 0; i < a->n; i++) {
		if (node_in_list(node_get(i), reachable)) {
			new_q[i] = k;
			k++;
		}
		else {
			new_q[i] = -1;
		}
	}

	node* new_start = NULL;
	for (node* n = a->start; n != NULL; n = n->next) {
		if (new_q[n->q] != -1) {
			new_start = list_add(new_start, node_get(new_q[n->q]));
		}
	}
	node* new_end = NULL;
	for (node* n = a->end; n != NULL; n = n->next) {
		if (new_q[n->q] != -1) {
			new_end = list_add(new_end, node_get(new_q[n->q]));
		}
	}
	nfa* b = nfa_init(a->dim, k, new_start, new_end);
	for (int i = 0; i < a->n; i++) {
		for (int symb = 0; symb < (1 << a->dim); symb++) {
			for (node* n = a->g->adj_list[i].symbols[symb].head; n != NULL; n = n->next) {
				if (new_q[i] != -1 && new_q[n->q] != -1) {
					nfa_add(b, new_q[i], symb, new_q[n->q]);
				}
			}
		}
	}
	return b;
}

bool equal_states(node** partition, nfa* a, int q1, int q2, int t) {
	for (int symb = 0; symb < (1 << a->dim); symb++) {
		for (int i = 0; i < t; i++) {
			if (node_in_list(a->g->adj_list[q1].symbols[symb].head, partition[i]) ^
				node_in_list(a->g->adj_list[q2].symbols[symb].head, partition[i])) {
				return false;
			}
		}
	}
	return true;
}

node** new_partition(node** partition, nfa* a, int* t, bool* changed) {
	node** n_partition = (node**)malloc(0 * sizeof(node*));
	int t_start = *t;
	int p = 0;
	for (int i = 0; i < *t; i++) {
		node* s = partition[i];
		node** subpart = (node**)malloc(0 * sizeof(node*));
		int k = 0;
		for (node* n = s; n; n = n->next) {
			bool eq_found = false;
			for (int j = 0; j < k; j++) {
				if ((n->q != subpart[j]->q) && equal_states(partition, a, n->q, subpart[j]->q, *t)) {
					subpart[j] = list_add(subpart[j], node_get(n->q));
					eq_found = true;
					break;
				}
			}
			if (!eq_found) {
				subpart = (node**)realloc(subpart, (k + 1) * sizeof(node*));
				subpart[k] = node_get(n->q);
				k++;
			}
		}
		n_partition = (node**)realloc(n_partition, (p + k) * sizeof(node*));
		for (int j = p; j < p + k; j++) {
			n_partition[j] = subpart[j - p];
		}
		p += k;
	}
	*t = p;
	if (t_start == *t) {
		*changed = false;
	}
	return n_partition;
}

nfa* nfa_minimize(nfa* a) {
	a = nfa_del_unrechable(a);
	node** partition = (node**)malloc(2 * sizeof(node*));
	partition[0] = a->end;
	partition[1] = NULL;
	int t = 2;
	for (int i = 0; i < a->n; i++) {
		if (!node_in_list(node_get(i), a->end)) {
			partition[1] = list_add(partition[1], node_get(i));
		}
	}

	bool changed = true;
	while (true) {
		node** n_partition = new_partition(partition, a, &t, &changed);
		if (!changed) {
			break;
		}
		partition = n_partition;
	}

	node* start = NULL;
	node* end = NULL;
	for (int i = 0; i < t; i++) {
		for (node* n = a->start; n; n = n->next) {
			if (node_in_list(n, partition[i])) {
				start = list_add(start, node_get(i));
			}
		}
		for (node* n = a->end; n; n = n->next) {
			if (node_in_list(n, partition[i])) {
				end = list_add(end, node_get(i));
			}
		}
	}

	nfa* m = nfa_init(a->dim, t, start, end);
	for (int i = 0; i < t; i++) {
		for (node* n = partition[i]; n; n = n->next) {
			for (int symb = 0; symb < (1 << a->dim); symb++) {
				for (int j = 0; j < t; j++) {
					if (a->g->adj_list[n->q].symbols[symb].head && node_in_list(a->g->adj_list[n->q].symbols[symb].head, partition[j])) {
						nfa_add(m, i, symb, j);
					}
				}
			}
		}
	}
	return m;
}

nfa* miss_trans(nfa* a)
{
	nfa* b = nfa_init(a->dim, a->n + 1, NULL, NULL);
	node* curr = NULL;

	for (int i = 0; i < a->n; ++i)
	{
		for (int symb = 0; symb < (1 << a->dim); ++symb)
		{
			curr = a->g->adj_list[i].symbols[symb].head;
			if (curr) {
				for (; curr; curr = curr->next)
					nfa_add(b, i, symb, curr->q);
			}
			else nfa_add(b, i, symb, b->n - 1);
		}
	}

	for (int symb = 0; symb < (1 << a->dim); ++symb)
		nfa_add(b, b->n - 1, symb, b->n - 1);

	node* start = NULL;
	node* end = NULL;

	for (node* nd = a->start; nd; nd = nd->next)
		start = list_add(start, node_get(nd->q));

	for (node* nd = a->end; nd; nd = nd->next)
		end = list_add(end, node_get(nd->q));

	b->start = start;
	b->end = end;
	return b; // nfa_minimize(nfa_to_dfa(b)); idk is it must have
}

nfa* nfa_to_dfa(nfa* a) {
	if (nfa_is_dfa(a))
		return nfa_copy(a);

	int start_num = 0;

	for (node* curr = a->start; curr; curr = curr->next)
		start_num += (1 << (curr->q));

	node* start = node_get(start_num);
	nfa* result = nfa_init(a->dim, (1 << a->n), start, NULL);
	node* end = NULL;

	for (node* curr = a->end; curr; curr = curr->next) {
		for (int num = 1; num < result->n; ++num) {
			if (num & (1 << (curr->q)))
				end = list_add(end, node_get(num));
		}
	}

	result->end = end;

	int pos, trans, q;

	for (int i = 1; i < result->n; i++) {
		for (int symb = 0; symb < (1 << a->dim); symb++) {
			trans = 0;
			for (int k = 0; (i >> k); k++) {
				if (!((1 << k) & i)) continue;
				for (node* curr = a->g->adj_list[k].symbols[symb].head; curr; curr = curr->next) {
					pos = (1 << (curr->q));
					trans |= pos;
				}

			}
			if (trans) nfa_add(result, i, symb, trans);
		}
	}

	return result;
}

nfa* nfa_cartesian(nfa* n1, nfa* n2) {
	int c1, c2, q1, q2;

	nfa* new_n = nfa_init(n1->dim, n1->n * n2->n, NULL, NULL);

	for (int symb = 0; symb < (1 << n1->dim); symb++) {
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

nfa* nfa_intersect(nfa* a, nfa* b) {

	nfa* n1 = miss_trans(a);
	nfa* n2 = miss_trans(b);
	nfa* new_n = nfa_cartesian(n1, n2);

	int state_num = 0;
	node* new_start = NULL;
	for (node* start1 = n1->start; start1; start1 = start1->next) {
		for (node* start2 = n2->start; start2; start2 = start2->next) {
			state_num = start2->q * n1->n + start1->q;
			new_start = list_add(new_start, node_get(state_num));
		}
	}
	new_n->start = new_start;

	node* new_end = NULL;
	for (node* end1 = n1->end; end1; end1 = end1->next) {
		for (node* end2 = n2->end; end2; end2 = end2->next) {
			state_num = end2->q * n1->n + end1->q;
			new_end = list_add(new_end, node_get(state_num));
		}
	}
	new_n->end = new_end;
	nfa_free(n1);
	nfa_free(n2);
	return new_n;//nfa_minimize(nfa_to_dfa(new_n));
}

nfa* nfa_union(nfa* a, nfa* b) {

	//first we need to check if there are states 
	//from which there are no transitions acros all symbols
	//add_miss_trans add missing transitions to given NFA
	//same for intersection

	nfa* n1 = miss_trans(a);
	nfa* n2 = miss_trans(b);
	nfa* new_n = nfa_cartesian(n1, n2);

	node* new_start = NULL;
	int state_num = 0;
	for (node* start1 = n1->start; start1; start1 = start1->next) {
		for (node* start2 = n2->start; start2; start2 = start2->next) {
			state_num = start2->q * (n1->n) + start1->q;
			new_start = list_add(new_start, node_get(state_num));
		}
	}

	new_n->start = new_start;

	node* new_end = NULL;
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

	nfa_free(n1);
	nfa_free(n2);

	return new_n;
	return nfa_minimize(new_n);// nfa_minimize(nfa_to_dfa(new_n));
}

nfa* nfa_union_minimal(nfa* a, nfa* b) {
	nfa* result = nfa_union(a, b);
	nfa* temp = nfa_to_dfa(result);
	nfa_free(result);
	result = nfa_minimize(temp), nfa_free(temp);
	return result;
}

nfa* nfa_intersect_minimal(nfa* a, nfa* b) {
	nfa* result = nfa_intersect(a, b);
	nfa* temp = nfa_to_dfa(result);
	nfa_free(result);
	result = nfa_minimize(temp), nfa_free(temp);
	return result;
}

nfa* nfa_complement(nfa* a) { //my version (it works)
	nfa* c = miss_trans(a);
	nfa* b = nfa_to_dfa(c);
	nfa_free(c);
	c = nfa_minimize(b);
	int fl;
	node* new_end = NULL;

	for (int i = 0; i < c->n; ++i)
	{
		fl = 1;
		for (node* end = c->end; end; end = end->next)
		{
			if (i == end->q) {
				fl = 0;
				break;
			}

		}
		if (fl) new_end = list_add(new_end, node_get(i));
	}

	b->end = new_end;
	list_free(c->end);
	c->end = new_end;

	return c;

}

int nfa_is_dfa(nfa* n) {
	node* nd_1;
	node* nd_2;

	for (int i = 0; i < n->n; i++) {
		for (int j = 0; j < pow(2, n->dim); j++) {

			nd_1 = n->g->adj_list[i].symbols[j].head;
			if (nd_1) {
				nd_2 = nd_1->next;
				if (nd_2) return 0;
			}
		}
	}
	return 1;
}

nfa* nfa_projection(nfa* a, int n) {
	n = a->dim - n - 1;
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
				new_symb = ((symb >> (n + 1)) << n) + (((1 << n) - 1) & symb);
				nfa_add(b, i, new_symb, nd->q);
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
				int left = ((symb - right) << 1);
				nfa_add(b, i, left + 0 + right, nd->q);
				nfa_add(b, i, left + (1 << n) + right, nd->q);
			}
		}
	}
	return b;
}

nfa* nfa_swap(nfa* n, int i, int j) {
	node* start = NULL;
	for (node* x = n->start; x; x = x->next) {
		start = list_add(start, node_get(x->q));
	}
	node* end = NULL;
	for (node* x = n->end; x; x = x->next) {
		end = list_add(end, node_get(x->q));
	}
	nfa* new_n = nfa_init(n->dim, n->n, start, end);

	int p = max(i, j), q = min(i, j);
	int t1 = (1 << p), t2 = (1 << q), new_symb;

	for (int k = 0; k < n->n; k++) {
		for (int symb = 0; symb < (1 << n->dim); symb++) {
			new_symb = symb - (t1 & symb) - (t2 & symb) + ((t1 & symb) >> (p - q)) + ((t2 & symb) << (p - q));
			for (node* current = n->g->adj_list[k].symbols[symb].head; current; current = current->next) {
				nfa_add(new_n, k, new_symb, current->q);
			}
		}
	}
	return new_n;
}


/// <summary>
/// Composes an automaton for a sum of the right-hand-sides
/// </summary>
/// <param name="a">an automaton for y = a_1*x</param>
/// <param name="b">an automaton for y = b_1*x</param>
/// <returns>an automaton for y = (a_1+b_1)*x</returns>
nfa* nfa_sum_equals(nfa* a, nfa* b) {
	if (!a) {
		return b;
	}
	nfa* temp = nfa_to_dfa(a);
	nfa* u = nfa_minimize(temp);
	u = nfa_extend(u, 1);
	u = nfa_extend(u, 3);
	u = nfa_extend(u, 3);

	nfa_free(temp);
	temp = nfa_to_dfa(b);
	nfa* v = nfa_minimize(temp);
	v = nfa_extend(v, 1);
	v = nfa_extend(v, 1);
	v = nfa_extend(v, 4);

	nfa* w = nfa_read("automata_lib\\sum.txt");
	w = nfa_extend(w, 0);
	w = nfa_extend(w, 0);
	w = nfa_swap(w, 2, 4);

	nfa* eq = nfa_read("automata_lib\\equals.txt");
	eq = nfa_extend(eq, 0);
	eq = nfa_extend(eq, 2);
	eq = nfa_extend(eq, 2);

	u = nfa_intersect_minimal(u, v);
	u = nfa_intersect_minimal(u, w);
	u = nfa_intersect_minimal(u, eq);

	u = nfa_projection(u, 2); // if i try to add nfa_to_dfa without minimization program fails
	u = nfa_projection(u, 2); // if i try to add nfa_to_dfa with minimization - wrong answer
	u = nfa_projection(u, 1); // both situations tested for linear_equals(2)

	nfa_free(temp), nfa_free(v), nfa_free(w), nfa_free(eq);
	temp = nfa_to_dfa(u), nfa_free(u);
	u = nfa_minimize(temp), nfa_free(temp);
	return u;
}

nfa* nfa_cut_leading_zeros(nfa* a) {
	node* start = node_get(0);
	nfa* zeros = nfa_init(a->dim, 1, start, start);
	nfa_add(zeros, 0, 0, 0);
	return nfa_left_quot(a, zeros);
}

/// <summary>
/// automaton for y = a * x
/// </summary>
/// <param name="a">coefficinet</param>
/// <returns></returns>
nfa* nfa_linear_equals(int a) {
	int k = 0;
	for (; (a >> k) > 0; k++);
	nfa** deg2 = (nfa**)malloc(k * sizeof(nfa*));
	deg2[0] = nfa_read("automata_lib\\equals.txt"); // x = y
	for (int i = 1; i < k; i++) {
		deg2[i] = nfa_sum_equals(deg2[i - 1], deg2[i - 1]); // x = (2^k)*y
	}

	nfa* ans = NULL;
	bool fl = false;
	for (int i = 0; (a >> i) > 0; i++) {
		if (((a >> i) & 1) == 1) {
			if (fl) {
				ans = nfa_sum_equals(ans, deg2[i]);
			}
			else {
				ans = deg2[i];
				fl = true;
			}
		}
	}
	ans = nfa_swap(ans, 1, 0);
	return ans;
}

void nfa_dfs(nfa* a, int q, int n, int* vis) {
	if (q >= a->n) return;

	int state_num;

	for (int symb = 0; symb < (1 << a->dim); symb++) {
		for (node* curr = a->g->adj_list[q].symbols[symb].head; curr; curr = curr->next) {
			state_num = curr->q * n + q;
			if (!vis[curr->q]) {
				vis[curr->q] = 1;
				nfa_dfs(a, curr->q, n, vis);
			}
		}
	}
}

nfa* nfa_left_quot(nfa* a, nfa* b) {
	int* vis = NULL;
	node* initial_start = a->start;
	nfa* lq = nfa_copy(a);
	nfa* n = NULL;
	int fl;

	list_free(lq->end);

	lq->end = NULL;

	for (int i = 0; i < a->n; i++) {
		fl = 0;

		a->start = node_get(i);
		n = nfa_intersect(a, b);

		for (node* curr = n->end; curr; curr = curr->next)
		{
			if (i == curr->q) fl = 1;
		}

		if (fl) {
			lq->end = list_add(lq->end, node_get(i));
			continue;
		}

		// if L(n) is non-empty -> add i to final states of lq.
		// I.e. just check non-emptiness: construct a set of reachable states R(i) in n from its strating state.
		// if R(i) intersect F(n) != 0 -> add i to final states of lq.

		vis = (int*)calloc(a->n * b->n, sizeof(int));
		nfa_dfs(n, i, a->n, vis);

		for (node* curr = n->end; curr; curr = curr->next) {
			if (vis[curr->q]) {
				lq->end = list_add(lq->end, node_get(i));
				break;
			}
		}

		free(vis);
		nfa_free(n);
		free(a->start);
	}
	a->start = initial_start;
	return lq;//nfa_minimize(nfa_to_dfa(lq));
}

nfa* right_quot(nfa* a, nfa* b) {
	int* vis = NULL;
	int state_num, s;
	node* initial_end = a->end;
	nfa* rq = nfa_copy(a);
	nfa* n = NULL;

	free(rq->start);

	rq->start = NULL;

	for (int i = 0; i < a->n; i++) {
		a->end = node_get(i);
		n = nfa_intersect(a, b);
		s = 0;


		for (node* curr = n->start; curr; curr = curr->next) {
			vis = (int*)calloc(a->n * b->n, sizeof(int));
			nfa_dfs(n, curr->q, a->n, vis);

			for (int j = 0; j < b->n; j++) {
				state_num = j * a->n + i;
				s += vis[state_num];
			}
			if (s) {
				rq->start = list_add(rq->start, node_get(i));
				free(vis);
				break;
			}
			free(vis);
		}
		nfa_free(n);
		free(a->end);
	}
	a->end = initial_end;
	return rq;//nfa_minimize(nfa_to_dfa(rq));
}

nfa* nfa_by_word(const char* word, int size)
{
	int symb;
	nfa* a = nfa_init(1, size + 1, node_get(0), node_get(size));

	for (int i = 0; i < size; ++i)
	{
		symb = (int)(word[size - i - 1] - '0');
		nfa_add(a, i, symb, i + 1);
	}
	return a;
}

void nfa_closure(nfa* a) //a+
{
	node* curr = NULL;

	for (node* end = a->end; end; end = end->next)
	{
		for (int symb = 0; symb < (1 << a->dim); ++symb)
		{
			for (node* start = a->start; start; start = start->next)
			{
				for (node* curr = a->g->adj_list[start->q].symbols[symb].head; curr; curr = curr->next)
					nfa_add(a, end->q, symb, curr->q);
			}
		}
	}
}


//* - 0 or more
//? - 0 or 1

nfa* add_empty_word(nfa* a) //a?
{
	nfa* b = nfa_init(0, 1, node_get(0), node_get(0)); //"empty" automata that accepts only empty word
	nfa* c = nfa_extend(b, 0);
	return nfa_union(a, c);
}

nfa* get_regex(const char* input, char** p)
{
	int len = strlen(input), size = 0, k = 0;
	char* str = NULL;
	nfa* curr = NULL;
	nfa** nfas = NULL;
	stack* s = stack_init();

	char* symb1 = (char*)malloc(1); //(
	char* symb2 = (char*)malloc(1); //|

	*symb1 = '(';
	*symb2 = '|';

	char* symb = (char*)input;

	for (; *symb != '\0'; ++symb)
	{

		switch (*symb)
		{
		case '(':
			stack_push(s, symb1);
			break;

		case ')':
			if (stack_is_empty(s))
			{
				*p = symb;

				for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
				free(nfas); free(symb1);
				free(symb2); stack_free(s);
				return NULL;
			}

			if (s->top->str[0] == '(')
			{
				if (!curr)
				{
					curr = nfa_by_word(str, size);
					free(str); str = NULL;
					size = 0;
				}
				else
				{
					*p = symb - size - 1;

					for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
					free(nfas); free(symb1);
					free(symb2); stack_free(s);
					return NULL;
				}
			}
			else
			{
				if (!curr) {
					*p = symb;

					for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
					free(nfas); free(symb1);
					free(symb2); stack_free(s);
					return NULL;
				}
				curr = nfa_union(curr, nfas[k - 1]);
				nfa_to_dot(curr, "test.dot");
				nfas = (nfa**)realloc(nfas, (k - 1) * sizeof(nfa*));
				k--;

				while (s->top->str[0] != '(' && !stack_is_empty(s))
					stack_pop(s);
				if (stack_is_empty(s))
				{
					*p = symb;

					for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
					free(nfas); free(symb1);
					free(symb2); stack_free(s);
					return NULL;
				}
			}
			stack_pop(s);
			break;

		case '|':
			if (s->top->str[0] == '|' && !curr)
			{
				*p = symb;
				for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
				free(nfas); free(symb1);
				free(symb2); stack_free(s);
				return NULL;
			}
			nfas = (nfa**)realloc(nfas, (k + 1) * sizeof(nfa*));
			nfas[k] = curr;
			curr = NULL;
			k++;
			stack_push(s, symb2);
			break;

		case '+':
			if (*(symb - 1) != ')')
			{
				*p = symb;
				for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
				free(nfas); free(symb1);
				free(symb2); stack_free(s);
				return NULL;
			}
			nfa_closure(curr);
			break;

		case '?':
			if (*(symb - 1) != ')')
			{
				*p = symb;
				for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
				free(nfas); free(symb1);
				free(symb2); stack_free(s);
				return NULL;
			}
			curr = add_empty_word(curr);
			break;

		case '*':
			if (*(symb - 1) != ')')
			{
				*p = symb;
				for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
				free(nfas); free(symb1);
				free(symb2); stack_free(s);
				return NULL;
			}
			nfa_closure(curr);
			curr = add_empty_word(curr);
			break;

		case ' ':
			break;

		default:
			if (!str)
				str = (char*)malloc(1);
			else
				str = (char*)realloc(str, size + 1);
			str[size] = *symb;
			size++;
			break;
		}
	}

	if (!stack_is_empty(s))
	{
		*p = symb;
		for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
		free(nfas); free(symb1);
		free(symb2); stack_free(s);
		return NULL;
	}

	for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
	free(nfas); free(symb1);
	free(symb2); stack_free(s);

	return curr;
}

nfa* nfa_from_regex(const char* s)
{
	char* p = NULL;
	nfa* a = get_regex(s, &p);
	if (p)
	{
		cout << "Invalid input" << endl;
		cout << s << endl;
		for (char* symb = (char*)s; symb != p; ++symb) cout << '~';
		cout << '^';
		return NULL;
	}
	return a;
}