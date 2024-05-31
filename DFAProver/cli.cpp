#include "cli.h"

using namespace std;

const char* automata_folder = "automata_lib\\";

char* substr(char* string, int start, int end) {
	char* sub = (char*)malloc((end - start + 2) * sizeof(char));
	for (int i = start; i <= end; i++) {
		sub[i - start] = string[i];
	}
	sub[end - start + 1] = '\0';
	return sub;
}

char* int_to_str(int x) {
	if (x == 0) {
		char* a = (char*)malloc(2 * sizeof(char));
		return a;
	}
	int k = 0;
	int t = x;
	while (x != 0) {
		k++;
		x /= 10;
	}
	char* a = (char*)malloc((k + 1) * sizeof(char));
	int c = 0;
	while (t != 0) {
		a[k - c - 1] = char('0' + (t % 10));
		t /= 10;
		c++;
	}
	a[k] = '\0';
	return a;
}

char* str_sum(char* a, char* b) {
	int c = strlen(a) + strlen(b);
	char* s = (char*)malloc((c + 1) * sizeof(char));
	for (int i = 0; i < strlen(a); i++) {
		s[i] = a[i];
	}
	for (int i = 0; i < strlen(b); i++) {
		s[strlen(a) + i] = b[i];
	}
	s[c] = '\0';
	return s;
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
		// def theorem_magic "Ax0 Ax1 (~$div2(x0) | ~$div2(x1) | Ex2 ($div2()))"
		// def sum_even "Ex2 ($div2(x2) & $sum(x0,x1,x2))"
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
			stack_push(rpn, int_to_str(p - input));
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

nfa* rpn_to_nfa(stack* rpn, nfa*** nfas, char*** names, int k) {
	nfa** st = (nfa**)malloc(0);
	int c = 0;
	if (rpn->top->str == "error") {
		stack_pop(rpn);
		char* numb = stack_top(rpn);
		cout << "IkbalProver: unknown symbol at position ";  
		printf("%s\n", numb);
		return NULL;
	}
	while (rpn->size != 0) {
		char* token = stack_top(rpn);
		if (!oprt(token)) {
			c++;
			st = (nfa**)realloc(st, c * sizeof(nfa*));
			nfa* a = NULL;
			for (int i = 0; i < k; i++) {
				if (str_eq((*names)[i], token)) {
					a = (*nfas)[i];
				}
			}
			if (!a) {
				cout << "IkbalProver: no automata with name ";
				printf("%s\n", token);
				return NULL;
			}
			st[c - 1] = a;
		}
		else if (token[0] == '&') {
			if (c <= 1) {
				cout << "IkbalProver: something wrong with formula" << endl;
				return NULL;
			}
			nfa* a = nfa_intersect(st[c - 2], st[c - 1]);
			st = (nfa**)realloc(st, (c - 1) * sizeof(nfa*));
			st[c - 2] = a;
			c--;
		}
		else if (token[0] == '|') {
			if (c <= 1) {
				cout << "IkbalProver: something wrong with formula" << endl;
				return NULL;
			}
			nfa* a = nfa_union(st[c - 1], st[c - 2]);
			st = (nfa**)realloc(st, (c - 1) * sizeof(nfa*));
			st[c - 2] = a;
			c--;
		}
		else if (token[0] == '~') {
			if (c <= 0) {
				cout << "IkbalProver: something wrong with formula" << endl;
				return NULL;
			}
			nfa* a = nfa_complement(st[c - 1]);
			st[c - 1] = a;
		}
		else if (token[0] == 'E') {
			if (c <= 0) {
				cout << "IkbalProver: something wrong with formula" << endl;
				return NULL;
			}
			nfa* a = st[c - 1];
			int coord = token[2] - '0';
			a = nfa_projection(a, coord);
			st[c - 1] = a;
		}
		else if (token[0] == 'A') {
			if (c <= 0) {
				cout << "IkbalProver: something wrong with formula" << endl;
				return NULL;
			}
			nfa* a = st[c - 1];
			int coord = token[2] - '0';
			a = nfa_complement(a);
			a = nfa_projection(a, coord);
			nfa_to_dot(a, "123.txt");
			a = nfa_complement(a);
			st[c - 1] = a;
		}
		stack_pop(rpn);
	}
	if (c != 1) {
		cout << "IkbalProver: something wrong with formula" << endl;
		return NULL;
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

int parse_input(char* input, nfa*** nfas, char*** names, int* k) {
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
			(*names)[*k] = substr(input, p, t - 1);
			for (p = t + 1, t += 2; input[t] != '\"'; t++);
			*nfas = (nfa**)realloc(*nfas, (*k + 1) * sizeof(nfa*));
			nfa* a = rpn_to_nfa(infix_to_rpn(substr(input, p + 1, t - 1)), nfas, names, *k);
			if (!a) {
				*names = (char**)realloc(*names, (*k) * sizeof(char*));
				*nfas = (nfa**)realloc(*nfas, (*k) * sizeof(nfa*));
			}
			else {
				(*nfas)[*k] = a;
				(*k)++;
			}
		}
		else {
			for (p = t + 1, t += 2; input[t] != '\"'; t++);
			nfa* a = rpn_to_nfa(infix_to_rpn(substr(input, p + 1, t - 1)), nfas, names, *k);
			if (a) {
				(*nfas)[redef] = a;
			}
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
	else if (str_eq(command, "exit")) {
		return 1;
	}
	else {
		cout << "Inavalid input. See \"help\" command" << endl;
	}
	return 0;
}

// def leq "Ex0 ($sum(x0,x1,x2))"
// eval "$leq(2,3)"
// def test "Ax0 ($div3(x0) & $div2(x0))"
// eval "$test()"
// def test "Ax0 ($div2(x0) & ~$div2(x0))"
// eval "$test()"

void cli() {
	struct _finddata_t c_file;
	intptr_t hFile;

	nfa** nfas = (nfa**)malloc(0 * sizeof(nfa));
	char** names = (char**)malloc(0 * sizeof(char));
	char* path = (char*)automata_folder;
	//strcat(path, "*.txt");
	hFile = _findfirst("automata_lib\\*.txt", &c_file);
	int automata_path_len = strlen(automata_folder);
	int i = 0;
	do {
		names = (char**)realloc(names, (i + 1) * sizeof(char*));
		nfas = (nfa**)realloc(nfas, (i + 1) * sizeof(nfa*));
		int j = 0;
		for (; j < 0xFF && c_file.name[j] != '.'; j++);
		names[i] = substr(c_file.name, 0, j - 1);

		char* n = substr((char*)automata_folder, 0, automata_path_len);
		n = (char*)realloc(n, strlen(automata_folder) + (j + 4) * (sizeof(char)));
		for (int h = automata_path_len; h < automata_path_len + j; h++) {
			n[h] = c_file.name[h - automata_path_len];
		}
		n[automata_path_len + j] = '.', n[14 + j] = 't', n[15 + j] = '.x', n[16 + j] = 't', n[17 + j] = '\0';
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
		int x = parse_input(input, &nfas, &names, &k);
		if (x == 1) {
			break;
		}
	}
}

nfa* nfa_by_word(const char* word, int size)
{
	int symb;
	nfa* a = nfa_init(1, size + 1, node_get(0), node_get(size));

	for (int i = 0; i < size; ++i)
	{
		symb = (int)(word[size - i - 1] - '0');
		nfa_add(a, i, symb, i + 1);
	}
	return a;
}

void nfa_closure(nfa* a) //a+
{
	node* curr = NULL;

	for (node* end = a->end; end; end = end->next)
	{
		for (int symb = 0; symb < (1 << a->dim); ++symb)
		{
			for (node* start = a->start; start; start = start->next)
			{
				for (node* curr = a->g->adj_list[start->q].symbols[symb].head; curr; curr = curr->next)
					nfa_add(a, end->q, symb, curr->q);
			}
		}
	}
}


//* - 0 or more
//? - 0 or 1

nfa* add_empty_word(nfa* a) //a?
{
	nfa* b = nfa_init(0, 1, node_get(0), node_get(0)); //"empty" automata that accepts only empty word
	nfa* c = nfa_extend(b, 0);
	return nfa_union(a, c);
}

nfa* get_regex(const char* input, char** p)
{
	int len = strlen(input), size = 0, k = 0;
	char* str = NULL;
	nfa* curr = NULL;
	nfa** nfas = NULL;
	stack* s = stack_init();

	char* symb1 = (char*)malloc(1); //(
	char* symb2 = (char*)malloc(1); //|

	*symb1 = '(';
	*symb2 = '|';

	char* symb = (char *) input;

	for (; *symb != '\0'; ++symb)
	{

		switch (*symb)
		{
		case '(':
			stack_push(s, symb1);
			break;

		case ')':
			if (stack_is_empty(s))
			{
				*p = symb;

				for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
				free(nfas); free(symb1);
				free(symb2); stack_free(s);
				return NULL;
			}

			if (s->top->str[0] == '(')
			{
				if (!curr)
				{
					curr = nfa_by_word(str, size);
					free(str); str = NULL;
					size = 0;
				}
				else
				{
					*p = symb - size - 1;

					for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
					free(nfas); free(symb1);
					free(symb2); stack_free(s);
					return NULL;
				}
			}
			else
			{
				if (!curr) {
					*p = symb;

					for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
					free(nfas); free(symb1);
					free(symb2); stack_free(s);
					return NULL;
				}
				curr = nfa_union(curr, nfas[k - 1]);
				nfa_to_dot(curr, "test.dot");
				nfas = (nfa**)realloc(nfas, (k - 1) * sizeof(nfa*));
				k--;

				while (s->top->str[0] != '(' && !stack_is_empty(s))
					stack_pop(s);
				if (stack_is_empty(s))
				{
					*p = symb;

					for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
					free(nfas); free(symb1);
					free(symb2); stack_free(s);
					return NULL;
				}
			}
			stack_pop(s);
			break;

		case '|':
			if (s->top->str[0] == '|' && !curr)
			{
				*p = symb;
				for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
				free(nfas); free(symb1);
				free(symb2); stack_free(s);
				return NULL;
			}
			nfas = (nfa**)realloc(nfas, (k + 1) * sizeof(nfa*));
			nfas[k] = curr;
			curr = NULL;
			k++;
			stack_push(s, symb2);
			break;

		case '+':
			if (*(symb - 1) != ')')
			{
				*p = symb;
				for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
				free(nfas); free(symb1);
				free(symb2); stack_free(s);
				return NULL;
			}
			nfa_closure(curr);
			break;

		case '?':
			if (*(symb - 1) != ')')
			{
				*p = symb;
				for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
				free(nfas); free(symb1);
				free(symb2); stack_free(s);
				return NULL;
			}
			curr = add_empty_word(curr);
			break;

		case '*':
			if (*(symb - 1) != ')')
			{
				*p = symb;
				for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
				free(nfas); free(symb1);
				free(symb2); stack_free(s);
				return NULL;
			}
			nfa_closure(curr);
			curr = add_empty_word(curr);
			break;

		case ' ':
			break;

		default:
			if (!str)
				str = (char*)malloc(1);
			else
				str = (char*)realloc(str, size + 1);
			str[size] = *symb;
			size++;
			break;
		}
	}

	if (!stack_is_empty(s))
	{
		*p = symb;
		for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
		free(nfas); free(symb1);
		free(symb2); stack_free(s);
		return NULL;
	}

	for (int i = 0; i < size; ++i) nfa_free(nfas[k]);
	free(nfas); free(symb1);
	free(symb2); stack_free(s);

	return curr;
}

nfa* nfa_from_regex(const char* s)
{
	char* p = NULL;
	nfa* a = get_regex(s, &p);
	if (p)
	{
		cout << "Invalid input" << endl;
		cout << s << endl;
		for (char* symb = (char*)s; symb != p; ++symb) cout << '~';
		cout << '^';
		return NULL;
	}
	return a;
}