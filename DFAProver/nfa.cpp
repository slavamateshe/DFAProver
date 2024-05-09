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
	for (node* curr = NFA->end; curr;  curr = curr->next)
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

int nfa_check(nfa* NFA, int* str) {
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
	node* nd = NULL;
	for (int symb = 0; symb < (1 << a->dim); symb++) {
		for (int i = 0; i < t; i++) {
			nd = a->g->adj_list[q1].symbols[symb].head;
			if (nd && (node_in_list(nd, partition[i]) ^
				node_in_list(nd, partition[i]))) {
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
	node* nd = NULL;
	for (int i = 0; i < t; i++) {
		for (node* n = partition[i]; n; n = n->next) {
			for (int symb = 0; symb < (1 << a->dim); symb++) {
				for (int j = 0; j < t; j++) {
					nd = a->g->adj_list[n->q].symbols[symb].head;
					if (nd && node_in_list(nd, partition[j])) {
						nfa_add(m, i, symb, j);
					}
				}
			}
		}
	}
	return m;
}

nfa* nfa_to_dfa(nfa* a) {
	int start_num = 0;

	for (node* curr = a->start; curr; curr = curr->next)
		start_num += (1 << (curr->q));

	node* start = node_get(start_num);
	nfa* result = nfa_init(a->dim, 1 << a->n, start, NULL);
	node* end = NULL;

	for (node* curr = a->end; curr; curr = curr->next) {
		for (int num = 0; num < result->n; ++num) {
			if (num & (1 << (curr->q)))
				end = list_add(end, node_get(num));
		}
	}

	result->end = end;

	int pos, trans, t;

	for (int i = 0; i < (1 << a->n); i++) {
		for (int symb = 0; symb < (1 << a->dim); symb++) {
			trans = 0;
			for (int j = 0; j < a->n; j++) {
				t = (1 << j & i);
				if (!t) continue;
				for (node* curr = a->g->adj_list[j].symbols[symb].head; curr; curr = curr->next) {
					pos = (1 << (curr->q));
					if (!(pos & trans)) trans += pos;
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

nfa* nfa_intersect(nfa* n1, nfa* n2) {
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

	return nfa_minimize(nfa_to_dfa(new_n));
}

nfa* nfa_union(nfa* n1, nfa* n2) {
	nfa* new_n = nfa_cartesian(n1, n2);

	node* new_start = NULL;
	int state_num = 0;
	for (node* start1 = n1->start; start1; start1 = start1->next) {
		for (node* start2 = n2->start; start2; start2 = start2->next) {
			state_num = start2->q * n1->n + start1->q;
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

	return nfa_minimize(nfa_to_dfa(new_n));
}

nfa* nfa_complement(nfa* a) {
	nfa* b = nfa_init(a->dim, a->n, a->start, NULL);
	b->g = a->g;
	node* end = NULL;
	for (int i = 0; i < b->n; i++) {
		bool fl = true;
		for (node* nd = a->end; nd; nd = nd->next) {
			if (i == nd->q) {
				fl = false;
			}
		}
		if (fl) {
			end = list_add(end, node_get(i));
		}
	}
	b->end = end;
	return nfa_minimize(nfa_to_dfa(b));
}

int nfa_is_dfa(nfa* n) {
	node* nd_1;
	node* nd_2;

	for (int i = 0; i < n->n; i++) {
		for (int j = 0; j < pow(2, n->dim); j++) {

			nd_1 = n->g->adj_list[i].symbols[j].head;
			if (nd_1) {
				nd_2 = nd_1->next;
				if (nd_2) {
					cout << i << " " << j << endl;
					return 0;
				}
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
	return nfa_minimize(nfa_to_dfa(b));
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
	nfa* c = nfa_to_dfa(b);
	return nfa_minimize(c);
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
	return nfa_minimize(nfa_to_dfa(new_n));
}

/// <summary>
/// Composes an automaton for a sum of the right-hand-sides
/// </summary>
/// <param name="a">an automaton for y = a_1*x</param>
/// <param name="b">an automaton for y = b_1*x</param>
/// <returns>an automaton for y = (a_1+b_1)*x</returns>
nfa* nfa_sum_equals(nfa* a, nfa* b) {
	nfa* u = a;
	u = nfa_extend(u, 0);
	u = nfa_extend(u, 1);
	u = nfa_extend(u, 4);
	u = nfa_swap(u, 2, 3);

	nfa* v = b;
	v = nfa_extend(v, 0);
	v = nfa_extend(v, 2);
	v = nfa_extend(v, 4);
	v = nfa_swap(v, 1, 3);

	nfa* w = nfa_read("sum.txt");
	w = nfa_extend(w, 3);
	w = nfa_extend(w, 4);

	nfa* eq = nfa_read("equals.txt");
	eq = nfa_extend(eq, 1);
	eq = nfa_extend(eq, 1);
	eq = nfa_extend(eq, 1);

	u = nfa_intersect(u, v);
	u = nfa_intersect(u, w);
	u = nfa_intersect(u, eq);

	u = nfa_projection(u, 0);
	u = nfa_projection(u, 0);
	u = nfa_projection(u, 0);

	return u;
}

nfa* nfa_cut_leading_zeros(nfa *a) {
	node* start = node_get(0);
	nfa* zeros = nfa_init(a->dim, 1, start, start);
	nfa_add(zeros, 0, 0, 0);
	return left_quot(a, zeros);
	//left_quot
}

/// <summary>
/// returns an automaton for (y = 2*x /\ a(x)) if a is a unary automaton
/// </summary>
/// <param name="a"></param>
/// <returns></returns>
nfa* nfa_double(nfa* a) {
	nfa* result = NULL;
	if (a->dim == 1) {
		node* start = node_get(0);
		nfa* mul2 = nfa_init(2, 2, start, start);
		nfa_add(mul2, 0, 0, 0);
		nfa_add(mul2, 0, 1, 1);
		nfa_add(mul2, 1, 2, 1);
		nfa_add(mul2, 1, 3, 0);
		nfa* a1 = nfa_extend(a, 0);
		nfa* conj = nfa_intersect(mul2, a1);
		nfa* semi_result = nfa_projection(conj, 0);
		nfa_free(a1);
		nfa_free(mul2);
	}
	return result;
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
	deg2[0] = nfa_read("equals.txt"); // x = y
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
	return ans;
}

void nfa_dfs(nfa* a, int q, int n, int* vis) {
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

nfa* left_quot(nfa* a, nfa* b) {
	int* vis = NULL;
	node* initial_start = a->start;
	nfa* lq = nfa_copy(a);
	nfa* n = NULL;

	free(lq->end);

	lq->end = NULL;

	for (int i = 0; i < a->n; i++) {
		a->start = node_get(i);
		n = nfa_intersect(a, b);

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
	return nfa_minimize(nfa_to_dfa(lq));
}

nfa* nfa_right_quot(nfa* a, nfa* b) {
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
	return nfa_minimize(nfa_to_dfa(rq));
}