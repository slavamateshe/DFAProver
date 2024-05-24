#pragma once

typedef struct node {
    int q;
    struct node* next;
} node;

typedef struct node_str {
    char* str;
    struct node_str* next;
} node_str;

typedef struct list_of_nodes {
    node* head;
} list_of_nodes;

typedef struct list_of_lists {
    list_of_nodes* symbols;
};

typedef struct graph {
    int count;
    list_of_lists* adj_list;
} graph;

typedef struct stack {
    node_str* top;
    int size;
} stack;

node* node_get(int q);
node* list_add(node* start, node* a);
void list_del(node* start, node* a);
bool node_in_list(node* a, node* list);
void list_free(node* n);

graph* graph_init(int n, int dim);
void graph_add_arc(graph* g, int q, int symb, int q_new);
void graph_del_arc(graph* g, int q, int symb, int q_del);
void graph_free(graph* g, int dim);

stack* stack_init();
void stack_push(stack* s, char* x);
void stack_pop(stack* s);
char* stack_top(stack* s);
int stack_is_empty(stack* s);