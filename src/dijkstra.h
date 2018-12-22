#ifndef DIJKSTRA_H
#define DIJKSTRA_H

typedef struct {
        int vertex;
        double weight;
} edge_t;

typedef struct {
        edge_t **edges;
        int edges_len;
        int edges_size;
        double dist;		/* needed for dij */
        int prev;		/* needed for dij */
        int visited;		/* needed for dij */
} vertex_t;

typedef struct {
        vertex_t **vertices;
        int vertices_len;
        int vertices_size;
} graph_t;

typedef struct {
        int *data;
        int *prio;
        int *index;
        int len;
        int size;
} heap_t;


typedef struct {
        int* path;
        int len;
        double score;
} dpath_t;

extern graph_t* alloc_graph(void);
extern void free_graph(graph_t* g);

extern int add_vertex (graph_t *g, int i);

extern int add_edge (graph_t *g, int a, int b, double w);

extern int remove_edge(graph_t *g, int a, int b);

extern void dijkstra (graph_t *g, int a, int b);

extern dpath_t* get_dijkstra_path (graph_t *g, int i);
extern void free_dijkstra_path(dpath_t* p);

#endif
