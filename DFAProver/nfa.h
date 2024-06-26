#pragma once

#include "data_struct.h"

typedef struct nfa {
    int dim;
    long long int n;
    graph* g;
    node* start;
    node* end;
} nfa;

nfa* nfa_init(int dim, int n, node* start, node* end);

void nfa_free(nfa* NFA);
void nfa_add(nfa* NFA, int q, int symb, int q_new);
void nfa_remove(nfa* NFA, int q, int symb, int q_del);
nfa* nfa_read(const char* s);
void nfa_write(nfa* NFA, const char* s);
void nfa_to_dot(nfa* NFA, const char* s);
int nfa_check(nfa* NFA, int* str);

nfa* nfa_cartesian(nfa* n1, nfa* n2);

nfa* nfa_intersect(nfa* n1, nfa* n2);

nfa* nfa_union(nfa* n1, nfa* n2);

int nfa_is_dfa(nfa* n);

nfa* nfa_extend(nfa* a, int n);
nfa* nfa_extend(nfa* a, int n);
nfa* nfa_projection(nfa* a, int n);

nfa* nfa_complement(nfa* a);

nfa* nfa_swap(nfa* n, int i, int j);

nfa* nfa_double(nfa* a);

nfa* power_of2(int k);

nfa* nfa_linear_equals(int a);
nfa* nfa_linear_equals(int *a, int n);
nfa* nfa_sum_equals(nfa* a, nfa* b);
nfa* nfa_del_unrechable(nfa* a);

nfa* nfa_sum_equals(nfa* a, nfa* b);

void nfa_dfs(nfa* a, int q, int n, int* vis);

nfa* nfa_minimize(nfa* x);
nfa* nfa_min(nfa* a);

stack* infix_to_rpn(char* input);
nfa* rpn_to_nfa(stack* rpn, nfa*** nfas, char*** names, int k);
char* substr(char* string, int start, int end);
nfa* nfa_right_quot(nfa* a, nfa* b);
nfa* nfa_left_quot(nfa* a, nfa* b);
nfa* nfa_cut_leading_zeros(nfa* a);

nfa* nfa_to_dfa(nfa* a);
nfa* nfa_union_minimal(nfa* a, nfa* b);
nfa* nfa_from_regex(const char* s);