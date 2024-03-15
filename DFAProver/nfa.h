#pragma once

#include "data_struct.h"

typedef struct nfa {
    int dim;
    int n;
    graph* g;
    node* start;
    node* end;
} nfa;

nfa* nfa_init(int dim, int n, node* start, node* end);

void nfa_free(nfa* NFA);
void nfa_add(nfa* NFA, int q, int symb, int q_new);
void nfa_remove(nfa* NFA, int q, int symb, int q_del);
nfa* nfa_read(const char* s);

void nfa_to_dot(nfa* NFA, const char* s);
int nfa_check(nfa* NFA, int str);

nfa* union_nfa(nfa* n1, nfa* n2);