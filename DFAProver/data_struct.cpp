#include <iostream>
#include <fstream>
#include "data_struct.h"

using namespace std;

node* node_get(int q) {
	node* edge = (node*)malloc(sizeof(node));
	edge->q = q;
	edge->next = NULL;
	return edge;
}

node_str* node_str_get(char* str) {
	node_str* edge = (node_str*)malloc(sizeof(node_str));
	edge->str = str;
	edge->next = NULL;
	return edge;
}

node* list_add(node* start, node* a) {
	if (!start) {
		start = node_get(a->q);
		return start;
	}
	node* n = start;
	node* t = start;
	int is_first = true;
	while (n) {
		if (n->q == a->q) {
			return start;
		}
		if (!is_first) {
			t = t->next;
		}
		n = n->next;
		is_first = false;
	}
	t->next = a;
	return start;
}

void list_del(node* start, node* a) {
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

bool node_in_list(node* a, node* list) {
	if (!a) {
		return false;
	}
	for (node* n = list; n != NULL; n = n->next) {
		if (a->q == n->q) {
			return true;
		}
	}
	return false;
}

void list_free(node* n) {
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
		g->adj_list[i].symbols = (list_of_nodes*)malloc(pow(2, dim) * sizeof(list_of_lists));
		for (int j = 0; j < pow(2, dim); j++) {
			g->adj_list[i].symbols[j].head = NULL;
		}
	}
	return g;
}


void graph_add_arc(graph* g, int q, int symb, int q_new) {
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

void graph_del_arc(graph* g, int q, int symb, int q_del) {
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
		for (int j = 0; j < (1 << dim); j++) {
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

stack* stack_init() {
	stack* s = (stack*)malloc(sizeof(stack));
	s->top = NULL;
	s->size = 0;
	return s;
}

void stack_push(stack* s, char* x) {
	node_str* n = node_str_get(x);
	if (!s->top) {
		s->top = n;
	}
	else
	{
		n->next = s->top;
		s->top = n;
	}
	s->size++;
}

void stack_pop(stack* s) {
	if (s->top) {
		if (s->top->next) {
			node_str* n = s->top->next;
			s->top = n;
		}
		else {
			s->top = NULL;
		}
		s->size--;
	}
}

char* stack_top(stack* s) {
	return s->top->str;
}