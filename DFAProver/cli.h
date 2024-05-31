#pragma once
#include "nfa.h"

void cli();

nfa* nfa_by_word(const char* word, int size);
void nfa_closure(nfa* a);
nfa* get_regex(const char* input, char **p);
void nfa_closure(nfa* a);
nfa* add_empty_word(nfa* a);
nfa* nfa_from_regex(const char* s);