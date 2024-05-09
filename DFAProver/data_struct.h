#pragma once

typedef struct node {
    int q;
    struct node* next;
} node;

typedef struct node_2 {
    int* q;
    int n;
    struct node_2* next;
};

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

node* node_get(int q);
node_2* node_get2(int* q, int n);
node* list_add(node* start, node* a);
node_2* list_add2(node_2* start, node_2* a);
void list_del(node* start, node* a);
void list_free(node* n);

graph* graph_init(int n, int dim);
void graph_add_arc(graph* g, int q, int symb, int q_new);
void graph_del_arc(graph* g, int q, int symb, int q_del);
void graph_free(graph* g, int dim);