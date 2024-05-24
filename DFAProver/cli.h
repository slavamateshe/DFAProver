#pragma once
#include "nfa.h"

void cli();

nfa* nfa_by_word(const char* word);
void nfa_closure(nfa* a);
nfa* hr_lang(const char* input);