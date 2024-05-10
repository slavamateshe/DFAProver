#pragma once

#define MAX_STACK_SIZE 100

typedef struct {
	unsigned char* st;
	int length;
}stack;

stack* stack_init();

void stack_push(stack* s, int a);

void stack_pop(stack* s);

void stack_free(stack* s);

int stack_is_empty(stack* s);
