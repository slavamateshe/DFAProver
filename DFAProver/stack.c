#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "stack.h"


stack* stack_init()
{
	stack* s = (stack*)malloc(sizeof(stack));
	s->st = (unsigned char*)calloc(MAX_STACK_SIZE, 1);
	s->length = 0;
	return s;
}

void stack_push(stack* s, int a) {
	if (!(s->length < MAX_STACK_SIZE))
		s->st = (unsigned char*)realloc(s->st, s->length + 1);
	s->st[s->length] = a;
	s->length++;
}

void stack_pop(stack* s) {
	if (s) {
		unsigned char* new_st = (unsigned char*)calloc(MAX_STACK_SIZE, sizeof(int));
		for (int i = 0; i < s->length - 1; i++) new_st[i] = s->st[i];
		free(s->st);
		s->st = new_st;
		s->length--;
	}
}

void stack_free(stack* s) {
	if (s) {
		free(s->st);
		free(s);
	}
}

int stack_is_empty(stack* s) {
	return (s->st != NULL);
}
