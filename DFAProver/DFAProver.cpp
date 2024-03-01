#include <iostream>
#include <fstream>

using namespace std;

typedef struct node {
    int q;
    struct node* next;
} node;

typedef struct list_of_nodes {
    node* head;
} list_of_nodes;

typedef struct list_of_lists {
	list_of_nodes* symbols;
};

typedef struct graph {
    int count;
    list_of_lists* adj_list;
} graph;

node* node_get(int q) {
	node* edge = (node*)malloc(sizeof(node));
	edge->q = q;
	edge->next = NULL;
	return edge;
}

void add_to_list(node* start, node* a) {
	if (!start) {
		start = a;
		return;
	}
	node* n = start;
	while (n) {
		n = n->next;
		if (n->q == a->q) {
			return;
		}
	}
	n->next = a;
}

void del_from_list(node* start, node* a) {
	node* n = start;
	while (n) {
		if (!n->next) {
			return;
		}
		if (n->next == a) {
			if (n->next->next) {
				n->next = n->next->next;
			}
			else {
				n->next = NULL;
			}
		}
	}
}

void free_list(node* n) {
	for (node* t = n->next; t; t = t->next) {
		free(n);
		n = t;
	}
	free(n);
}

graph* graph_init(int n, int dim) {
	graph* g = (graph*)malloc(sizeof(graph));
	g->count = n;
	g->adj_list = (list_of_lists*)malloc(n * sizeof(list_of_lists));
	for (int i = 0; i < n; i++) {
		g->adj_list[i].symbols = (list_of_nodes*)malloc(pow(2, n) * sizeof(list_of_lists));
		for (int j = 0; j < pow(2, dim); j++) {
			g->adj_list[i].symbols[j].head = NULL;
		}
	}
	return g;
}

void add_arc(graph* g, int q, int symb, int q_new) {
	list_of_nodes* l = g->adj_list[q].symbols;
	node* n = l[symb].head;

	if (n) {
		while (n->next) {
			if (n->q == q_new) return;
			n = n->next;
		}
		if (n->q != q_new) n->next = node_get(q_new);
	}
	else {
		g->adj_list[q].symbols[symb].head = node_get(q_new);
		g->adj_list[q].symbols[symb].head->next = NULL;
	}
}

void del_arc(graph* g, int q, int symb, int q_del) {
	node* n = g->adj_list[q].symbols[symb].head;
	if (n) {
		if (n->q == q_del) {
			if (n->next)
				g->adj_list[q].symbols[symb].head = n->next;
			else
				g->adj_list[q].symbols[symb].head = NULL;
			free(n);
		}
		else {
			while (n->next && n->next->q != q_del) {
				n = n->next;
			}
			if (n->next) {
				node* c = n->next;
				n->next = c->next;
				free(c);
			}
		}
	}
	return;
}


void graph_free(graph* g, int dim) {
	for (int i = 0; i < g->count; i++) {
		for (int j = 0; j < pow(2, dim); j++) {
			node* n = g->adj_list[i].symbols[j].head;
			if (n) {
				for (node* t = n->next; t; t = t->next) {
					free(n);
					n = t;
				}
				free(n);
			}
		}
		free(g->adj_list[i].symbols);
	}
	free(g->adj_list);
	free(g);
}


typedef struct nfa {
    int dim;
    int n;
    graph* g;
    node* start;
    node* end;
} nfa;


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
	add_arc(NFA->g, q, symb, q_new);
}

void nfa_remove(nfa* NFA, int q, int symb, int q_del) {
	del_arc(NFA->g, q, symb, q_del);
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
		add_to_list(start, node_get(x));
	}

	f >> end_len;
	node* end = NULL;
	for (int i = 0; i < start_len; i++) {
		f >> x;
		add_to_list(end, node_get(x));
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
		for (int j = 0; j < pow(2, NFA->dim); j++) {
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
	for(; str; str >>= 1) {
		node* n = current;
		node* qnew = NULL;
		while (n) {
			node* q = NFA->g->adj_list[current->q].symbols[str & 1].head;
			while (q) {
				add_to_list(qnew, q);
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


int main()
{
	nfa* a = nfa_read("3div.txt");
	nfa_to_dot(a, "3div.dot");
	int x;
	cin >> x;
	cout << "Work in progress...";//nfa_check(a, x);
}