#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <io.h>
#include <time.h>
#include "cli.h"

using namespace std;

char* substr(char* string, int start, int end) {
	char* sub = (char*)malloc((end - start + 2) * sizeof(char));
	for (int i = start; i <= end; i++) {
		sub[i - start] = string[i];
	}
	sub[end - start + 1] = '\0';
	return sub;
}

bool str_eq(char* a, const char* b) {
	if (strlen(a) != strlen(b)) {
		return false;
	}
	for (int i = 0; i < strlen(a); i++) {
		if (a[i] != b[i]) {
			return false;
		}
	}
	return true;
}

int prec(char x) {
	if (x == '~') return 3;
	if (x == 'E') return 3;
	if (x == 'A') return 3;
	if (x == '&') return 2;
	if (x == '|') return 1;
	return -1;
}


// nfa_get_lin(int *a, int n) -> automaton for y = a[0]+a[1]x_1+...+a[n-1]x_{n-1}
// 
// Ex ($div2(x+1) & $div2(x)) -> False
// $div2(3*x1+5*x2+1) -> Ey (div2(y) /\ y = 3*x1+5*x2+1)
stack* infix_to_rpn(char* input) {
	stack* operators = stack_init();
	stack* rpn = stack_init();
	char *p = input;
	while (*p != '\n' && *p != '\0') {
		char* token = (char*)malloc(2 * sizeof(char));
		token[0] = *p;
		token[1] = '\0';
		if (*p == ' ' || *p == '\n') {
			p++;
			continue;
		}
		if (*p == '$') {
			int t = 0; // parse name
			while (*(p+t) != '(') ++t;
			p++;
			strncpy(token, p, t - 1);
			token[t - 1] = '\0';
			p += (t + 1);
			stack_push(rpn, token);
			while (*p != ')') p++;
			p++;
			// parse linear term until ')'
			// start: ad hoc #1: just read (x)
			// !!!!!!!!!!!!!!!!! nfa_linear_equals// at least for a*x+b
		}
		else if (prec(token[0]) >= 0) {
			if (token[0] == 'E' || token[0] == 'A') {
				token = (char*)malloc(4 * sizeof(char));
				token[0] = *p;
				token[1] = *(p + 1);
				token[2] = *(p + 2);
				token[3] = '\0';
				p += 2;
			}
			while (operators->size != 0 && (stack_top(operators)[0] != '(') && (prec(stack_top(operators)[0]) > prec(token[0]))) {
				stack_push(rpn, stack_top(operators));
				stack_pop(operators);
			}
			stack_push(operators, token);
			p++;
		}
		else if (token[0] == '(') {
			stack_push(operators, token);
			p++;
		}
		else if (token[0] == ')') {
			while (stack_top(operators)[0] != '(') {
				stack_push(rpn, stack_top(operators));
				stack_pop(operators);
			}
			stack_pop(operators);
			p++;
		}
		else {
			while (rpn->size != 0) {
				stack_pop(rpn);
			}
			while (operators->size != 0) {
				stack_pop(operators);
			}
			stack_push(rpn, (char*)"error");
			break;
		}
	}
	while (operators->size != 0) {
		stack_push(rpn, stack_top(operators));
		stack_pop(operators);
	}
	stack* ans = stack_init();
	while (rpn->size != 0) {
		stack_push(ans, stack_top(rpn));
		stack_pop(rpn);
	}
	return ans;
}

char* handle_name(char* a) {
	int c = strlen(a);
	char* b = (char*)malloc((c + 5) * sizeof(char));
	strncpy(b, a, 4);
	b[c] = '.', b[c + 1] = 't', b[c + 2] = 'x', b[c + 3] = 't', b[c + 4] = '\0';
	return b;
}

bool oprt(char* x) {
	if (x[0] == '~' || x[0] == '|' || x[0] == '&' || x[0] == 'E' || x[0] == 'A') return true;
	return false;
}

nfa* rpn_to_nfa(stack* rpn) {
	nfa** st = (nfa**)malloc(0);
	int c = 0;
	if (rpn->size == 1 && rpn->top->str == "error") {
		cout << "error" << endl;
		return NULL;
	}
	while (rpn->size != 0) {
		char* token = stack_top(rpn);
		if (!oprt(token)) {
			c++;
			st = (nfa**)realloc(st, c * sizeof(nfa*));
			char* x = handle_name(stack_top(rpn));
			char* y = (char*)malloc((strlen("automata_lib\\") + strlen(x) + 2) * sizeof(char));
			y[strlen("automata_lib\\")] = '\0';
			strncpy(y, "automata_lib\\", strlen("automata_lib\\"));
			int l1 = strlen(y);
			for (int i = l1; i < l1 + strlen(x); i++) {
				y[i] = x[i - l1];
			}
			y[l1 + strlen(x)] = '\0';

			st[c - 1] = nfa_read(y);
		}
		else if (token[0] == '&') {
			nfa* a = nfa_intersect(st[c - 1], st[c - 2]);
			st = (nfa**)realloc(st, (c - 1) * sizeof(nfa*));
			st[c - 2] = a;
			c--;
		}
		else if (token[0] == '|') {
			nfa* a = nfa_union(st[c - 1], st[c - 2]);
			st = (nfa**)realloc(st, (c - 1) * sizeof(nfa*));
			st[c - 2] = a;
			c--;
		}
		else if (token[0] == '~') {
			nfa* a = nfa_complement(st[c - 1]);
			st[c - 1] = a;
		}
		else if (token[0] == 'E') {
			nfa* a = st[c - 1];
			int coord = token[2] - '0';
			a = nfa_projection(a, coord);
			st[c - 1] = a;
		}
		else if (token[0] == 'A') {
			nfa* a = st[c - 1];
			int coord = token[1] - '0';
			a = nfa_complement(a);
			a = nfa_projection(a, coord);
			a = nfa_complement(a);
			st[c - 1] = a;
		}
		stack_pop(rpn);
	}
	return st[0];
}

int* parse_argument(char* input) {
	if (strlen(input) == 0) {
		int* a = (int*)malloc(0);
		return a;
	}
	int k = 0;
	for (int i = 0; i < strlen(input); i++) {
		if (input[i] == ',') k++;
	}
	int* a = (int*)malloc((k + 1) * sizeof(int));

	int p = 0;
	for (int i = 0; i < k; i++) {
		int x = 0;
		for (int j = p; j < strlen(input) && input[j] != ','; j++, p++) {
			x = x * 10 + (input[j] - '0');
		}
		a[i] = x;
		p++;
	}
	int x = 0;
	for (int j = p; j < strlen(input); j++) {
		x = x * 10 + (input[j] - '0');
	}
	a[k] = x;

	return a;
}

void parse_input(char* input, nfa*** nfas, char*** names, int* k) {
	int p = 0;
	int t = p;
	for (; input[t] != ' ' && input[t] != '\n'; t++);
	char* command = substr(input, p, t - 1);
	if (str_eq(command, "def")) {
		for (p = t + 1, t++; input[t] != ' ' && input[t] != '\n'; t++);
		char* name = substr(input, p, t - 1);
		int redef = -1;
		for (int i = 0; i < *k; i++) {
			if (str_eq(name, (*names)[i])) {
				redef = i;
				break;
			}
		}
		if (redef == -1) {
			*names = (char**)realloc(*names, (*k + 1) * sizeof(char*));
			(*names)[*k] = substr(input, p, t - 1); // def test "~($div2(x) & $div3(x))"
			for (p = t + 1, t += 2; input[t] != '\"'; t++);
			*nfas = (nfa**)realloc(*nfas, (*k + 1) * sizeof(nfa*));
			(*nfas)[*k] = (rpn_to_nfa(infix_to_rpn(substr(input, p + 1, t - 1))));
			(*k)++;
		}
		else {
			for (p = t + 1, t += 2; input[t] != '\"'; t++);
			(*nfas)[redef] = (rpn_to_nfa(infix_to_rpn(substr(input, p + 1, t - 1))));
		}
	}
	else if (str_eq(command, "eval")) {
		p = t + 3;
		t = p;
		int id = -1;
		for (t = p; input[t] != '('; t++);
		for (int i = 0; i < *k; i++) {
			if (str_eq((*names)[i], substr(input, p, t - 1))) {
				id = i;
				break;
			}
		}
		for (p = t + 1; input[t] != ')'; t++);
		int* argument = parse_argument(substr(input, p, t - 1)); 
		if (nfa_check((*nfas)[id], argument)) {
			cout << "True" << endl;
		}
		else {
			cout << "False" << endl;
		}
	}
	else if (str_eq(command, "eval_for")) {
		p = t + 3;
		t = p;
		int id = -1;
		for (t = p; input[t] != '('; t++);
		for (int i = 0; i < *k; i++) {
			if (str_eq((*names)[i], substr(input, p, t - 1))) {
				id = i;
				break;
			}
		}
		for (p = t + 1; input[t] != ')'; t++);
		int* argument = parse_argument(substr(input, p, t - 1));
		for (int i = 0; i < 20; i++) {
			int input[1] = { i };
			cout << i << ": ";
			if (nfa_check((*nfas)[id], input)) {
				cout << "True" << endl;
			}
			else {
				cout << "False" << endl;
			}
		}
	}
	else if (str_eq(command, "help")) {
		cout << "In progress..." << endl;
	}
	else {
		cout << "Inavalid input. See \"help\" command" << endl;
	}
	return;
}

// def name "$div2(x) & ~$div3(x)"
// eval "$name(33)"
// def leq "Ex0 ($sum(x0,x1,x2))"
// eval "$leq(2,3)"
// def test "Ex0 ($div3(x0) & $div2(x0))"
// eval "$test()"
void cli() {
	struct _finddata_t c_file;
	intptr_t hFile;

	nfa** nfas = (nfa**)malloc(0 * sizeof(nfa));
	char** names = (char**)malloc(0 * sizeof(char));
	hFile = _findfirst("automata_lib\\*.txt", &c_file);
	int i = 0;
	do {
		names = (char**)realloc(names, (i + 1) * sizeof(char*));
		nfas = (nfa**)realloc(nfas, (i + 1) * sizeof(nfa*));
		int j = 0;
		for (; j < 260 && c_file.name[j] != '.'; j++);
		names[i] = substr(c_file.name, 0, j - 1);

		char* n = substr((char*)"automata_lib\\", 0, 13);
		n = (char*)realloc(n, 13 + (j + 4) * (sizeof(char)));
		for (int h = 13; h < 13 + j; h++) {
			n[h] = c_file.name[h - 13];
		}
		n[13 + j] = '.', n[14 + j] = 't', n[15 + j] = '.x', n[16 + j] = 't', n[17 + j] = '\0';
		nfas[i] = nfa_read(n);
		i++;
	} while (_findnext(hFile, &c_file) == 0);
	_findclose(hFile);

	int k = i;
	while (true) {
		cout << "IkbalProver: ";
		int buffsize = 128;
		char* input = (char*)malloc(buffsize * sizeof(char));
		fgets(input, buffsize, stdin);
		parse_input(input, &nfas, &names, &k);
	}
}