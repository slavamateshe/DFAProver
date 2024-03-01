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
    int* start;
    int* end;
} nfa;


nfa* nfa_init(int dim, int n, int* start, int* end) {
	nfa* NFA = (nfa*)malloc(sizeof(nfa));
	NFA->dim = dim;
	NFA->n = n;
	NFA->g = graph_init(n, dim);
	NFA->start = start;
	NFA->end = end;
	return NFA;
}

void* nfa_free(nfa* NFA) {
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
	int* start = (int*)malloc(sizeof(int) * start_len);
	for (int i = 0; i < start_len; i++) f >> start[i];

	f >> end_len;
	int* end = (int*)malloc(sizeof(int) * end_len);
	for (int i = 0; i < end_len; i++) f >> end[i];
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
	int* q = NFA->start;
	int* used = (int*)malloc(sizeof(int) * NFA->n);

	for (int i = 0; str > 0; i++, str /= 2) {
		for (int j = 0; j < NFA->n; j++) {
			used[j] = 0;
		}

		int* qn = (int*)malloc(0);
		int c = 0;
		for (int j = 0; j < sizeof(q) / sizeof(int); j++) {
			node* n = NFA->g->adj_list[q[j]].symbols[str % 2].head;
			while (n) {
				if (used[n->q] == 0) {
					used[n->q] = 1;
					c++;
					qn = (int*)realloc(qn, c);
					qn[c - 1] = n->q;
				}
				n->next;
			}
		}
		q = qn;
	}

	for (int i = 0; i < sizeof(NFA->end) / sizeof(int); i++) {
		if (used[NFA->end[i]]) {
			return 1;
		}
	}
	return 0;
}


int main()
{
	nfa* a = nfa_read("graph.txt");
	nfa_to_dot(a, "graph2.dot");
}


