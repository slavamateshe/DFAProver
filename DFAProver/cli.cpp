#include <iostream>
#include <fstream>
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
	int p = 0;
	while (p < strlen(input)) {
		char x = input[p];
		char* token = (char*)malloc(2 * sizeof(char));
		token[0] = x;
		token[1] = '\0';
		p++;
		if (x == ' ' || x == '\n') {
			continue;
		}
		if (x == '$') {
			int t = p; // parse name
			for (; input[t] != '('; t++);
			token = substr(input, p, t - 1);
			p = t + 1;
			stack_push(rpn, token);
			// parse linear term until ')'
			// start: ad hoc #1: just read (x)
			while (input[t] != ')') ++t;
			p = t + 1;
			//   end: ad hoc #1
		}
		else if (prec(token[0]) >= 0) {
			while (operators->size != 0 && (stack_top(operators)[0] != '(') && (prec(stack_top(operators)[0]) > prec(token[0]))) {
				stack_push(rpn, stack_top(operators));
				stack_pop(operators);
			}
			stack_push(operators, token);
		}
		else if (token[0] == '(') {
			stack_push(operators, token);
		}
		else if (token[0] == ')') {
			while (stack_top(operators)[0] != '(') {
				stack_push(rpn, stack_top(operators));
				stack_pop(operators);
			}
			stack_pop(operators);
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

char* handle_name(char* input) {
	int c = strlen(input);
	char* a = (char*)realloc(input, (c + 5) * sizeof(char));
	a[c] = '.', a[c + 1] = 't', a[c + 2] = 'x', a[c + 3] = 't', a[c + 4] = '\0';
	return a;
}

bool oprt(char* x) {
	if (x[0] == '~' || x[0] == '|' || x[0] == '&') return true;
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
			st[c - 1] = nfa_read(handle_name(stack_top(rpn)));
			stack_pop(rpn);
		}
		else if (token[0] == '&') {
			nfa* a = nfa_intersect(st[c - 1], st[c - 2]);
			st = (nfa**)realloc(st, (c - 1) * sizeof(nfa*));
			st[c - 2] = a;
			stack_pop(rpn);
			c--;
		}
		else if (token[0] == '|') {
			nfa* a = nfa_union(st[c - 1], st[c - 2]);
			st = (nfa**)realloc(st, (c - 1) * sizeof(nfa*));
			st[c - 2] = a;
			stack_pop(rpn);
			c--;
		}
		else if (token[0] == '~') {
			nfa* a = nfa_complement(st[c - 1]);
			st[c - 1] = a;
			stack_pop(rpn);
		}
	}
	return st[0];
}

int* parse_argument(char* input) {
	int* a = (int*)malloc(1 * sizeof(int));
	int x = 0;
	for (int i = 0; i < strlen(input); i++) {
		x = x * 10 + (input[i] - '0');
	}
	a[0] = x;
	return a;
}

void parse_input(char* input, nfa*** nfas, char*** names, int* k) {
	int p = 0;
	int t = p;
	for (; input[t] != ' ' && input[t] != '\n'; t++);
	if (str_eq(substr(input, p, t - 1), "def")) {
		for (p = t + 1, t++; input[t] != ' ' && input[t] != '\n'; t++);
		*names = (char**)realloc(*names, (*k + 1) * sizeof(char*));
		*names[*k] = substr(input, p, t - 1);
		for (p = t + 1, t += 2; input[t] != '\"'; t++);
		*nfas = (nfa**)realloc(*nfas, (*k + 1) * sizeof(nfa*));
		*nfas[*k] = rpn_to_nfa(infix_to_rpn(substr(input, p + 1, t - 1)));
		(*k)++;
	}
	else if (str_eq(substr(input, p, t - 1), "eval")) {
		p = t + 3;
		t = p;
		int id = -1;
		for (t = p; input[t] != '('; t++);
		for (int i = 0; i < *k; i++) {
			if (str_eq(*names[i], substr(input, p, t - 1))) {
				id = i;
				break;
			}
		}
		for (p = t + 1; input[t] != ')'; t++);
		int* argument = parse_argument(substr(input, p + 1, t - 1)); 
		if (nfa_check(*nfas[id], argument)) {
			cout << "True" << endl;
		}
		else {
			cout << "False" << endl;
		}
	}
	return;
}

// def name "$div2(x) & ~$div3(x)"
// eval "$name(33)"
void cli() {
	nfa** nfas = (nfa**)malloc(0 * sizeof(nfa));
	char** names = (char**)malloc(0 * sizeof(char));
	int k = 0;
	while (true) {
		cout << "IkbalProver: ";
		int buffsize = 128;
		char* input = (char*)malloc(buffsize * sizeof(char));
		fgets(input, buffsize, stdin);
		parse_input(input, &nfas, &names, &k);
	}
}