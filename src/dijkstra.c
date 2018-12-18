#include "ssrdt_net.h"

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
 
typedef struct {
        int vertex;
        float weight;
} edge_t;
 
typedef struct {
        edge_t **edges;
        int edges_len;
        int edges_size;
        float dist;		/* needed for dij */
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

void add_vertex (graph_t *g, int i);

void add_edge (graph_t *g, int a, int b, float w);

void remove_edge(graph_t *g, int a, int b);

void dijkstra (graph_t *g, int a, int b);
void dijkstra_local (graph_t *g, int a, int b,float max_distance);


graph_t * make_graph(struct rbtree_edge_struct** list,int num_edges);

void free_graph(graph_t* g);


int calculate_connectivity(graph_t* g, float* connectivity);


int progress_bar(float x);

/*
  Algorithm 1.
  NA algorithm
  Input: Aweighted graph G=(V,E), q and γ
  Output: Subgraph H⊂G
  1: Sort edges E by weights in an ascending order.
  2: F ← E
  3: n ← γ(|E|-(|V|−1))
  4: {Iteratively prune the weakest edge which does not cut the graph}
  5: i ← 1,j← 1 {j is an index to the sorted list of edges}
  6: while i ≤ n do
  7:    if C (u, v;F\{ej})is not −∞ then
  8:    F ← F \{ej}
  9:    i←i+1
  10:    j←j+1
  11: Return H=(V,F)
*/

int prune_edges_naive(struct rbtree_edge_struct** list,int num_edges,float gamma)
{
        struct rbtree_edge_struct* e = NULL;
        graph_t *g = NULL;
        int i;
        int a,b;
        int num_remove;
        int removed;
	
        RUNP(g = make_graph(list,num_edges));
        num_remove =  (int)(gamma * (float)(num_edges - (g->vertices_len-1)));
	
        LOG_MSG("Pruning.");
        LOG_MSG("%d edges.", num_edges);
        LOG_MSG("will remove %d edges.",num_remove);
	
        ///DECLARE_TIMER(t1);
        removed = 0;
        for(i = num_edges-1 ; i >= 0;i--){
                e = (struct rbtree_edge_struct*) list[i];
	         
                a = e->key >> 32ULL;
                b = e->key & 0xFFFFFFFF;
                //fprintf(stdout,"%d\t%d\t%f\n",a,b,e->dist);

                if(e->dist != -1){

                        // remove
                        remove_edge(g,a,b);
                        dijkstra(g, a,b);
                        if(g->vertices[b]->dist == FLT_MAX){
                                add_edge(g,a,b, e->dist );
                        }else{
                                e->dist = -1;
                                //		fprintf(stdout,"%d\t%d\t%f removed...\n",a,b,e->dist);
                                removed++;
				
					
                        }
		
                }
                if(removed == num_remove){
                        break;
                }
        }
	
        free_graph(g);

        return OK;
ERROR:
        return FAIL;
}




/*
  Algorithm 3.
  PS algorithm
  Input: A weighted graph G =(V,E),q and γ
  Output: Subgraph H⊂G
  1:  F←E
  2:  n←γ(|E|−(|V|−1))
  3:  { Iteratively prune the edge with the largest κ value. }
  4:  M←∅
  5:  for r=1 to n do
  6:  κlargest←−∞
  7:  elargest←null
  8:  for e={u, v} in F and e ∈ M do
  9:    Find path S such that q(S)=C(u,v;F\{e})
  10:   if q(S)≥q({e})then
  11:      κ←1
  12:      F←F\{e}
  13:      break
  14:   else if 0<q(S)<q({e})then
  15:      κ←q(S)/q({e})
  16:   else
  17:      κ←−∞
  18:      M←M+e
  19:   if κ>κlargest then
  20:      κlargest←κ
  21:      elargest←e
  22:   F←F\{elargest}
  23: Return H=(V,F)
*/

int prune_edges_path_simplification(struct rbtree_edge_struct** list,int num_edges,float gamma)
{
        struct rbtree_edge_struct* e = NULL;
        graph_t *g = NULL;
        int i,j;
        int a,b;
        int num_remove;
        int remove_target;
        float q_s;		/* path quality without edge */
        float q_e;		/* path quality with edge...  */
        float min_k;
        float k;
        
        RUNP(g = make_graph(list,num_edges));
        num_remove =  (int)(gamma * (float)(num_edges - (g->vertices_len-1)));
	
        LOG_MSG("Pruning.");
        LOG_MSG("%d edges.", num_edges);
        LOG_MSG("will remove %d edges.",num_remove);
	
        DECLARE_TIMER(t1);
        //num_remove = 2;
        for(i = 0; i < num_remove;i++){
                START_TIMER(t1);
                //RUN(calculate_connectivity(g, &connectivity));
                min_k = FLT_MAX;
                remove_target = -1;
                for(j = 0; j < num_edges;j++){
                        e = (struct rbtree_edge_struct*) list[j];
                        a = e->key >> 32ULL;
                        b = e->key & 0xFFFFFFFF;
                        if(e->dist != -1){
                                dijkstra_local(g, a,b,1.0);
                                q_e = g->vertices[b]->dist;
				
                                // remove
                                remove_edge(g,a,b);
				
                                // calculate new connectivity
                                dijkstra_local(g, a,b,1.0);
                                q_s = g->vertices[b]->dist;

                                //fprintf(stdout,"testing %d %d; org %f new %f. ",a,b,q_e,q_s);
					
                                if(q_s <= q_e){
                                        //	fprintf(stdout,"remove immediately\n");
                                        // remove immediately
                                        e->dist = -1;
			
                                        remove_target = -1;

                                        break;
                                }

                                if(q_s != -FLT_MAX){
                                        k = q_s / q_e;
                                        if(k < min_k){
						
                                                //fprintf(stdout,"new candidate\n");
                                                remove_target = j;
                                                min_k = k;
                                        }else{
						
                                                //		fprintf(stdout,"\n");
                                        }
                                }
                                add_edge(g,a,b, e->dist );
                        }
			
                        RUN(progress_bar( (float)(j+1)/ (float) num_edges));
                }
                if(remove_target != -1){
                        e = (struct rbtree_edge_struct*) list[remove_target];
                        a = e->key >> 32ULL;
                        b = e->key & 0xFFFFFFFF;
		
                        e->dist = -1;
                        remove_edge(g,a,b);
                }
                STOP_TIMER(t1);
	
                LOG_MSG("remove edge in %f seconds.", GET_TIMING(t1));
        }
	
        free_graph(g);
        return OK;
ERROR:
        return FAIL;
}


/*
  Algorithm 2.
  BF algorithm
  Input:
  Aweightedgraph G=(V,E),q and γ
  Output:
  Subgraph H ⊂ G
  1: F ← E
  2: n ← γ ( |E| − ( |V| − 1))
  3: { Iteratively prune the edge with the highest rk value. }
  4: M ←∅{ edges whose removal is known to cut the graph. }
  5: for r =1 to n do
  6: rk largest ←−∞
  7: e largest ← null
  8: for e={u, v} in F and  e ∈ M do
  9: if graph (V,F\{e}) is connected then
  10: compute rk(V,F,e)=C(V,F\{e})/ C(V,F)
  11: if rk(V,F,e) >rk largest then
  12: rklargest←rk(V,F,e)
  13: elargest←e
  14: else
  15: M←M+e
  16: F←F\{elargest}
  17: Return H=(V,F)
*/


int prune_edges_brute_force(struct rbtree_edge_struct** list,int num_edges,float gamma)
{
        struct rbtree_edge_struct* e = NULL;
        graph_t *g = NULL;
        int i,j;
        int a,b;
        int num_remove;
        int remove_target;
        float connectivity;
        float connectivity_minus_one;
        float min_rk;

        RUNP(g = make_graph(list,num_edges));
        num_remove =  (int)(gamma * (float)(num_edges - (g->vertices_len-1)));
	
		
        // 2: n ← γ ( |E| − ( |V| − 1))
        LOG_MSG("Pruning.");
        LOG_MSG("%d edges.", num_edges);
        LOG_MSG("will remove %d edges.",num_remove);

        DECLARE_TIMER(t1);
	
	
	
        //num_remove = 2;
        for(i = 0; i < num_remove;i++){
                START_TIMER(t1);
                RUN(calculate_connectivity(g, &connectivity));
                connectivity_minus_one = -FLT_MAX;
                min_rk = FLT_MAX;
                remove_target = -1;
                for(j = 0; j < num_edges;j++){
                        e = (struct rbtree_edge_struct*) list[j];
                        a = e->key >> 32ULL;
                        b = e->key & 0xFFFFFFFF;
                        if(e->dist != -1){

                                // remove
                                remove_edge(g,a,b);

                                // calculate new connectivity
                                connectivity_minus_one = 0;
                                RUN(calculate_connectivity(g, &connectivity_minus_one ));
                                //fprintf(stdout,"edge:%d %f %f (%f)\n",j,connectivity,connectivity_minus_one, connectivity_minus_one / connectivity);
                                // sort num diff
                                if( connectivity_minus_one / connectivity < min_rk){
                                        min_rk =  connectivity_minus_one / connectivity;
                                        remove_target = j;
                                }
                                add_edge(g,a,b, e->dist );
                        }
			
                        RUN(progress_bar( (float)(j+1)/ (float) num_edges));
                }
                e = (struct rbtree_edge_struct*) list[remove_target];
                a = e->key >> 32ULL;
                b = e->key & 0xFFFFFFFF;
		
                e->dist = -1;
                remove_edge(g,a,b);
		
                STOP_TIMER(t1);
	
                LOG_MSG("remove edge in %f seconds.", GET_TIMING(t1));
        }




	
        free_graph(g);
        return OK;
ERROR:
        return FAIL;
}

int progress_bar(float x)
{
        int barWidth = 70;
        fprintf(stdout,"[");
        int pos = barWidth * x;
        for (int i = 0; i < barWidth; ++i) {
                if (i < pos){
                        fprintf(stdout,"=");

                } else if (i == pos){
                        fprintf(stdout,">");
	
                }else{
                        fprintf(stdout," ");
                }
        }

	
        fprintf(stdout,"] %d \r", (int)(x * 100.0f));
        fflush(stdout);
        return OK;
}
 
int calculate_connectivity(graph_t* g, float* connectivity)
{
        int i,j;
        int len = g->vertices_len;
        float c = 0;
        for(i = 0; i < len-1;i++){
                for(j = i+1;j < len;j++){
                        dijkstra(g, i,j);
                        if(g->vertices[j]->dist == FLT_MAX){
                                *connectivity = FLT_MAX;
                                return OK;
                        }
                        c += g->vertices[j]->dist;
                }
        }
        c =  2.0 / ((float)(len *(len-1))) * c;

        *connectivity = c;
        return OK;
}

graph_t * make_graph(struct rbtree_edge_struct** list,int num_edges)
{
        struct rbtree_edge_struct* e = NULL;

        uint32_t* count_v = NULL;
        graph_t *g = NULL;

        int num_vertices = -1;
        int max = -1;

        int i;
        int a,b;
        RUNP(count_v = make_bitvector(num_edges*2));
        for(i = 0;i < num_edges;i++){
                e = (struct rbtree_edge_struct*) list[i];
                a = (uint32_t)(e->key >> 32ull);
                b = (uint32_t)(e->key & 0xFFFFFFFFull);

                bit_set(count_v,a);
                bit_set(count_v,b);
		
                if(a > max){
                        max = a;
                }
                if(b > max){
                        max = b;
                }
        }
        num_vertices = 0;
        for(i = 0; i <= max;i++){
                if(bit_test(count_v,i)){
                        num_vertices++;
                }
        }
	
        // build while graph..
        g = calloc(1, sizeof (graph_t));

        for(i = 0; i <= max;i++){
                if(bit_test(count_v,i)){
                        add_vertex(g,i);
                }else{
                        LOG_MSG("%d not found!!!\n",i);
			       
                }
        }

        for(i = 0;i < num_edges;i++){
                e = (struct rbtree_edge_struct*) list[i];
                a = e->key >> 32ULL;
                b = e->key & 0xFFFFFFFF;
                add_edge(g, a, b, e->dist );
        }
        MFREE(count_v);

        return g;
ERROR:
        return NULL;
		
	       
}

void free_graph(graph_t* g)
{
        int i,j;
        /* free stuff */
        for(i = 0; i < g->vertices_len;i++){
                vertex_t *v = g->vertices[i];
                for(j = 0; j < v->edges_len;j++){
                        free(v->edges[j]);
                }
                free(v->edges);
                free(v);
		
        }
        free(g->vertices);
        free(g);

}


void add_vertex (graph_t *g, int i)
{
        if (g->vertices_size < i + 1) {
                int size = g->vertices_size * 2 > i ? g->vertices_size * 2 : i + 4;
                g->vertices = realloc(g->vertices, size * sizeof (vertex_t *));
                for (int j = g->vertices_size; j < size; j++){
                        g->vertices[j] = NULL;
                }
                g->vertices_size = size;
        }
        if (!g->vertices[i]) {
                g->vertices[i] = calloc(1, sizeof (vertex_t));
                g->vertices_len++;
        }
}
 
void add_edge (graph_t *g, int a, int b, float w)
{

        vertex_t *v = NULL; //a = a - 'a';
        edge_t* e = NULL;
        //b = b - 'a';
        //add_vertex(g, a);
        //add_vertex(g, b);
	
        v = g->vertices[a];
        /*c = -1;
          for(i = 0; i < v->edges_len;i++){
          if(v->edges[i]->vertex == b){
          v->edges[i]->weight += w;
          c = i;
          break;
          }
          }*/
        //if(c == -1){
        if (v->edges_len >= v->edges_size) {
                v->edges_size = v->edges_size ? v->edges_size * 2 : 4;
                v->edges = realloc(v->edges, v->edges_size * sizeof (edge_t *));
        }
        e = calloc(1, sizeof (edge_t));
        e->vertex = b;
        e->weight = w;
        v->edges[v->edges_len++] = e;
        //}
        //v = g->vertices[b];
        /*c = -1;
          for(i = 0; i < v->edges_len;i++){
          if(v->edges[i]->vertex == a){
          v->edges[i]->weight += w;
          c = i;
          break;
          }
          }
          if(c == -1){
        */
        //if (v->edges_len >= v->edges_size) {
        //        v->edges_size = v->edges_size ? v->edges_size * 2 : 4;
        //        v->edges = realloc(v->edges, v->edges_size * sizeof (edge_t *));
        //}
        //e = calloc(1, sizeof (edge_t));
        //e->vertex = a;
        //e->weight = w;
        //v->edges[v->edges_len++] = e;
        //}
	
}

void remove_edge(graph_t *g, int a, int b)
{
        int i;
        int c;
        //a = a - 'a';
        //b = b - 'a';
        //add_vertex(g, a);
        //add_vertex(g, b);
        c = -1;
        vertex_t *v = g->vertices[a];
        for(i = 0; i < v->edges_len;i++){
                if(v->edges[i]->vertex == b){
                        //fprintf(stdout,"found v\n");
                        c = i;
                        break;
                }
		
        }
        if(c == -1){
                fprintf(stdout,"V not found!!!!\n");
        }
        free(v->edges[c]);
        v->edges[c] = v->edges[v->edges_len-1];
        v->edges_len--;

        /*c = -1;
        v = g->vertices[b];
        for(i = 0; i < v->edges_len;i++){
                if(v->edges[i]->vertex == a){
                        //fprintf(stdout,"found v\n");
                        c = i;
                        break;
                }
		
        }
        if(c == -1){
                fprintf(stdout,"V not found!!!!\n");
        }
        free(v->edges[c]);
        v->edges[c] = v->edges[v->edges_len-1];
        v->edges_len--;*/
}
 
heap_t *create_heap (int n) {
        heap_t *h = calloc(1, sizeof (heap_t));
        h->data = calloc(n + 1, sizeof (int));
        h->prio = calloc(n + 1, sizeof (int));
        h->index = calloc(n, sizeof (int));
        return h;
}

void free_heap(heap_t* h)
{
        if(h){
                free(h->data);
                free(h->prio);
                free(h->index);
                free(h);
        }
}

 
void push_heap (heap_t *h, int v, int p) {
        int i = h->index[v] == 0 ? ++h->len : h->index[v];
        int j = i / 2;
        while (i > 1) {
                if (h->prio[j] < p){
                        break;
                }
                h->data[i] = h->data[j];
                h->prio[i] = h->prio[j];
                h->index[h->data[i]] = i;
                i = j;
                j = j / 2;
        }
        h->data[i] = v;
        h->prio[i] = p;
        h->index[v] = i;
}
 
int min (heap_t *h, int i, int j, int k) {
        int m = i;
        if (j <= h->len && h->prio[j] < h->prio[m])
                m = j;
        if (k <= h->len && h->prio[k] < h->prio[m])
                m = k;
        return m;
}
 
int pop_heap (heap_t *h) {
        int v = h->data[1];
        int i = 1;
        while (1) {
                int j = min(h, h->len, 2 * i, 2 * i + 1);
                if (j == h->len){
                        break;
                }
                h->data[i] = h->data[j];
                h->prio[i] = h->prio[j];
                h->index[h->data[i]] = i;
                i = j;
        }
        h->data[i] = h->data[h->len];
        h->prio[i] = h->prio[h->len];
        h->index[h->data[i]] = i;
        h->len--;
        return v;
}
 
void dijkstra (graph_t *g, int a, int b) {
        int i, j;
        //a = a - 'a';
        //b = b - 'a';
        for (i = 0; i < g->vertices_len; i++) {
                vertex_t *v = g->vertices[i];
                v->dist = FLT_MAX;
                v->prev = 0;
                v->visited = 0;
        }
        vertex_t *v = g->vertices[a];
        v->dist = 0;
        heap_t *h = create_heap(g->vertices_len);
        push_heap(h, a, v->dist);
	
        while (h->len) {
                i = pop_heap(h);
                if (i == b){
                        break;
                }
                v = g->vertices[i];
                v->visited = 1;
                for (j = 0; j < v->edges_len; j++) {
                        edge_t *e = v->edges[j];
                        vertex_t *u = g->vertices[e->vertex];
			
                        if (!u->visited && v->dist + e->weight <= u->dist) {
                                u->prev = i;
                                u->dist = v->dist + e->weight;

                                push_heap(h, e->vertex, u->dist);
                        }
                }
        }
        free_heap(h);
}

void dijkstra_local (graph_t *g, int a, int b,float max_distance)
{
        int i, j;
        //a = a - 'a';
        //b = b - 'a';
        for (i = 0; i < g->vertices_len; i++) {
                vertex_t *v = g->vertices[i];
                v->dist = FLT_MAX;
                //v->prev = 0;
                v->visited = 0;
        }
        vertex_t *v = g->vertices[a];
        v->dist = 0;
        heap_t *h = create_heap(g->vertices_len);
        push_heap(h, a, v->dist);
	
        while (h->len) {
                i = pop_heap(h);
                if (i == b){
                        break;
                }
                v = g->vertices[i];
                v->visited = 1;
                for (j = 0; j < v->edges_len; j++) {
                        edge_t *e = v->edges[j];
                        vertex_t *u = g->vertices[e->vertex];
                        if(v->dist + e->weight <= max_distance){
                                if (!u->visited && v->dist + e->weight <= u->dist) {
                                        u->prev = i;
                                        u->dist = v->dist + e->weight;

                                        push_heap(h, e->vertex, u->dist);
                                }
                        }
                }
        }
        free_heap(h);
}


int print_path (graph_t *g, int i) {
        int n, j;
        vertex_t *v, *u;
        //i = i - 'a';
        v = g->vertices[i];
        /*if (v->dist == INT_MAX) {
                printf("no path\n");
                return 0;
        }
        int c = 0;
        for (n = 1, u = v; u->dist; u = g->vertices[u->prev], n++){
                fprintf(stdout,"%d ", u->prev);
                c++;
                if(c > 30){
                        fprintf(stdout,"\n");
                        exit(0);
                }
        }
        fprintf(stdout,"\n");*/
        /*char *path = malloc(n);
        path[n - 1] = 'a' + i;
        for (j = 0, u = v; u->dist; u = g->vertices[u->prev], j++){
                path[n - j - 2] = 'a' + u->prev;
        }

        printf("%d %.*s\n", v->dist, n, path);
        free(path);*/
        return v->dist;
}


#ifdef ITEST
int main () {
        graph_t *g = calloc(1, sizeof (graph_t));
        int i,j;
        char nodes[] = "abcdef";
        int num_nodes = 6;
        for(i = 0; i < num_nodes;i++){
                add_vertex(g,i);
        }
        //for(i = 0; i < num_nodes-1;i++){
        //	for(j = i+1;j < num_nodes;j++){
        //		add_edge(g,i,j,random_int_zero_to_x(100));
        //	}
        //}
        // random_sample_start
	
        add_edge(g, 0, 1, 30);
        add_edge(g, 0, 1, 30);
        add_edge(g, 0, 2, 9);
        add_edge(g, 0, 5, 14);
        add_edge(g, 1, 2, 10);
        add_edge(g, 1, 3, 15);
        add_edge(g, 2, 3, 11);
        add_edge(g, 2, 5, 2);
        add_edge(g, 3, 4, 6);
        add_edge(g, 4, 5, 9);

        DECLARE_TIMER(t1);
        for(i = 0; i < num_nodes-1;i++){
                START_TIMER(t1);
                for(j = i+1;j < num_nodes;j++){
                        dijkstra(g, i,j);
			
                }
                STOP_TIMER(t1);
                LOG_MSG("done in %f seconds.",GET_TIMING(t1));
        }
	
        for(i = 0; i < num_nodes;i++){
                for(j = 0; j < num_nodes;j++){
                        if(i != j){
                                dijkstra(g, i,j);
                                fprintf(stdout,"%d->%d %d ",i,j, print_path(g, j));
                        }else{
                                fprintf(stdout,"0 ");
                        }
                        fprintf(stdout,"\n");
			
                }

        }
        fprintf(stdout,"\n");
        remove_edge(g,3,4);
        remove_edge(g,4,5);
	
        for(i = 0; i < num_nodes;i++){
                for(j = 0; j < num_nodes;j++){
                        if(i != j){
                                dijkstra(g, i,j);
                                fprintf(stdout,"%d->%d %d ",i,j, print_path(g, j));
                                //fprintf(stdout,"%d ", print_path(g, j));
                        }else{
                                fprintf(stdout,"0 ");
                        }
                         fprintf(stdout,"\n");
                }

        }
        fprintf(stdout,"\n");

//dijkstra(g, 'a', 'f');
//print_path(g, 'f');
//dijkstra(g, 'a', 'f');
//print_path(g, 'f');

        for(i = 0; i < g->vertices_len;i++){
                vertex_t *v = g->vertices[i];
                for(j = 0; j < v->edges_len;j++){
                        free(v->edges[j]);
                }
                free(v->edges);
                free(v);
		
        }
        free(g->vertices);
        free(g);
        return 0;
}




#endif // ITEST

