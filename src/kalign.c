#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "math.h"
#include <string.h>
//#include "io.h"
#include <ctype.h>

#include <zlib.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"

#include "kalign.h"
#ifdef ITEST
#endif

struct sequence_info{
        struct sequence_info* next;
        char* name;
        char* value;
};

struct aln_tree_node{
        struct aln_tree_node** links;
        int* internal_lables;
        int* path;
        int* profile;
        int* seq;
        int len;
        int done;
        int num;
};

struct bignode{
        struct bignode *next;
        unsigned int pos[NODESIZE];
        unsigned int num;
};

struct bignode* big_insert_hash(struct bignode *n,const unsigned int pos);
void big_remove_nodes(struct bignode *n);
void big_print_nodes(struct bignode *n);

struct alignment{
        struct sequence_info** si;
        unsigned int** sip;
        unsigned int* nsip;
        unsigned int* sl;
        unsigned int* lsn;
        int** s;
        char**seq;
        char** sn;
        int numseq;
        int numprofiles;
};

struct hirsch_mem{
        struct states* f;
        struct states* b;
        int starta;
        int startb;
        int enda;
        int endb;
        int size;
        int len_a;
        int len_b;
};

struct dp_matrix{
        struct states* s;
        void* tb_mem;
        char** tb;
        int x;
        int y;
};

struct states{
        float a;
        float ga;
        float gb;
        float x;
};

struct aln_param{
        float** subm;
        float gpo;
        float gpe;
        float tgpe;
};


static int** dna_alignment(struct alignment* aln,int* tree,float**submatrix, int** map);


static int* mirror_hirsch_path(int* hirsch_path,int len_a,int len_b);
static int* add_gap_info_to_hirsch_path(int* hirsch_path,int len_a,int len_b);
static float* dna_update(const float* profa, const float* profb, float* newp,int* path,int sipa,int sipb);

static int* hirsch_dna_ss_dyn(float**subm, const int* seq1,const int* seq2,struct hirsch_mem* hm, int* hirsch_path);
static int* hirsch_align_two_dna_ss_vector(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[]);
static struct states* foward_hirsch_dna_ss_dyn(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm);
static struct states* backward_hirsch_dna_ss_dyn(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm);

static int* hirsch_dna_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm, int* hirsch_path,int sip);
static int* hirsch_align_two_dna_ps_vector(const float* prof1,const int* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],int sip);
static struct states* foward_hirsch_dna_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm,int sip);
static struct states* backward_hirsch_dna_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm,int sip);

static int* hirsch_dna_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm, int* hirsch_path);
static int* hirsch_align_two_dna_pp_vector(const float* prof1,const float* prof2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[]);
static struct states* foward_hirsch_dna_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm);
static struct states* backward_hirsch_dna_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm);


static struct hirsch_mem* hirsch_mem_alloc(struct hirsch_mem* hm,int x);
static struct hirsch_mem* hirsch_mem_realloc(struct hirsch_mem* hm,int x);
static void hirsch_mem_free(struct hirsch_mem* hm);

static float* dna_make_profile(float* prof, int* seq,int len,float** subm);
static void dna_set_gap_penalties(float* prof,int len,int nsip);
static float* dna_update(const float* profa, const float* profb, float* newp,int* path,int sipa,int sipb);
static struct alignment* make_seq(struct alignment* aln,int a,int b,int* path);
static void update_gaps(int old_len,int*gis,int new_len,int *newgaps);

static char* get_input_into_string(char* string,char* infile);
void fasta_output(struct alignment* aln,char* outfile);
static int count_sequences_fasta(char* string);
static struct alignment* make_dna(struct alignment* aln);
static struct alignment* aln_alloc(int numseq);
static void free_aln(struct alignment* aln);
static struct alignment* read_sequences(struct alignment* aln,char* string);
struct alignment* detect_and_read_sequences(char* infile);

static float dna_distance_calculation(struct bignode* hash[],int* p,int seqlen,int diagonals,float mode);
static float** dna_distance(struct alignment* si,float** dm);
static struct aln_tree_node* real_upgma(float **dm,int numseq);
static int* readtree(struct aln_tree_node* p,int* tree);


char** kalign_align(char** sequences, int numseq)
{
        struct alignment* aln = NULL;
        float** subm = NULL;
        int* tree = NULL;
        int i,j;
        int a, b, c, f, tmp;

        float** dm = NULL;
        struct aln_tree_node* tree2 = NULL;
        int** map = NULL;
        int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
        subm = malloc(sizeof (float*) * 32);
        for (i = 32;i--;){
                subm[i] = malloc(sizeof(float) * 32);
                for (j = 32;j--;){
                        subm[i][j] = 0.0; //0;//gpe << 1;//-5;// better on Balibase
                }
        }

        for(i = 0; i < 4;i++){
                for(j = 0; j < 4;j++){
                        subm[i][j] = 1;
                        if(i ==j){
                                subm[i][j] += 2;
                        }

                }

        }


        /*
          104	-95	-95	-113
          -95	139	-78	-95
          -95	-78	139	-95
          -113	-95	-95	104

        */
        //setup;
        if(numseq > 1){
                aln =  aln_alloc(numseq);
                for(i = 0; i < numseq;i++){
                        a = (int) strlen(sequences[i]);
                        aln->seq[i] = malloc(sizeof(char) * (a+1));
                        aln->s[i] = malloc(sizeof(int) * (a+1));
                        aln->sl[i] = a;
                        aln->sn[i] =malloc(sizeof(char) * (1)); // blank so free doesn't fall over later...
                        //fprintf(stderr,"Sequence %d len: %d\n", i, a);
                        //fprintf(stderr,"%s\n",sequences[i]);
                        for(j = 0; j < a;j++){
                                //	fprintf(stderr,"%c	%d\n",sequences[i][j],j);
                                aln->seq[i][j] = sequences[i][j];

                                aln->s[i][j] =      aacode[toupper(sequences[i][j])-65];
                        }
                        aln->seq[i][a] = 0;
                        aln->s[i][a] = 0;
                }

                aln = make_dna(aln);

                tree = malloc(sizeof(int)*(aln->numseq*3+1));
                for ( i = 1; i < (aln->numseq*3)+1;i++){
                        tree[i] = 0;
                }
                tree[0] = 1;

                if(aln->numseq < 10000){
                        //	fprintf(stdout,"Started estimating distances.\n");
                        dm =  dna_distance(aln,dm);
                        //	fprintf(stdout,"Started tree building\n");
                        tree2 = real_upgma(dm,aln->numseq);

                        tree = readtree(tree2,tree);
                        for (i = 0; i < (aln->numseq*3);i++){

                                tree[i] = tree[i+1];
                        }
                        free(tree2->links);
                        free(tree2->internal_lables);
                        free(tree2);
                }else{

                        tree[0] = 0;
                        tree[1] = 1;

                        c = aln->numseq;
                        tree[2] = c;
                        a = 2;
                        for ( i = 3; i < (aln->numseq-1)*3;i+=3){
                                tree[i] = c;
                                tree[i+1] = a;
                                c++;
                                tree[i+2] = c;
                                a++;
                        }
                }

                //fprintf(stdout,"Started alignment\n");

                map =  dna_alignment(aln,tree,subm, map);

                for(i = 0; i < aln->numseq;i++){
                        sequences[i] = realloc(sequences[i], sizeof(char) * (map[aln->numprofiles-1 ][0]+ 1) );
                }

                for (i = 0; i < aln->numseq;i++){
                        for (a = 0; a < aln->sl[i];a++){
                                aln->s[i][a] = 0;
                        }
                }
                //clear up
                for (i = 0; i < (aln->numseq-1)*3;i +=3){
                        a = tree[i];
                        b = tree[i+1];
                        aln = make_seq(aln,a,b,map[tree[i+2]]);
                }
                for (i = 0; i < aln->numseq;i++){
                        aln->nsip[i] = i;
                }

                // added -  make alignment in input sequence pointer
                for (i = 0; i < aln->numseq;i++){
                        f = aln->nsip[i];
                        //fprintf(fout,">%s\n",aln->sn[f]);
                        c = 0;
                        for (j = 0; j < aln->sl[f];j++){
                                tmp = aln->s[f][j];
                                while (tmp){
                                        sequences[i][c] = '-';
                                        c++;
                                        tmp--;
                                }
                                sequences[i][c] =aln->seq[f][j];
                                c++;
                        }
                        tmp = aln->s[f][aln->sl[f]];
                        while (tmp){
                                sequences[i][c] = '-';
                                c++;
                                tmp--;
                        }
                        sequences[i][c] = 0;
                }
                free(tree);
                free(map);
                free_aln(aln);
        }

        for (i = 32;i--;){
                free(subm[i]);// = malloc(sizeof(float) * 32);
        }
        free(subm);// = malloc(sizeof (float*) * 32);

        return sequences;
        //free(alignment_rotate_90);
}




float** dna_distance(struct alignment* si,float** dm)
{
        struct bignode* hash[1024];

        int *p = 0;
        int i,j,a;//,b;
        unsigned int hv;


        //fprintf(stderr,"Distance Calculation:\n");


        for (i = 0;i < 1024;i++){
                hash[i] = 0;
        }


        dm = malloc (sizeof(float*)*si->numseq);
        for (i = si->numseq;i--;){
                dm[i] = malloc (sizeof (float)*(si->numseq));
                for (j = si->numseq;j--;){
                        dm[i][j] = 0.0f;
                }
        }


        //b = (si->numseq*(si->numseq-1))/2;
        a = 1;

        for (i = 0; i < si->numseq-1;i++){
                p = si->s[i];
                for (j = si->sl[i]-5;j--;){
                        hv = ((p[j]&3)<<8) + ((p[j+1]&3)<<6) + ((p[j+2]&3)<<4)  + ((p[j+3]&3)<<2) + (p[j+4]&3);//ABCDE
                        hash[hv] = big_insert_hash(hash[hv],j);
                        hv = ((p[j]&3)<<8) + ((p[j+1]&3)<<6) + ((p[j+2]&3)<<4)  + ((p[j+3]&3)<<2) + (p[j+5]&3);//ABCDF
                        hash[hv] = big_insert_hash(hash[hv],j);
                        hv = ((p[j]&3)<<8) + ((p[j+1]&3)<<6) + ((p[j+2]&3)<<4)  + ((p[j+4]&3)<<2) + (p[j+5]&3);//ABCEF
                        hash[hv] = big_insert_hash(hash[hv],j);
                        hv = ((p[j]&3)<<8) + ((p[j+1]&3)<<6) + ((p[j+3]&3)<<4)  + ((p[j+4]&3)<<2) + (p[j+5]&3);//ABDEF
                        hash[hv] = big_insert_hash(hash[hv],j);
                        hv = ((p[j]&3)<<8) + ((p[j+2]&3)<<6) + ((p[j+3]&3)<<4) + ((p[j+4]&3)<<2) + (p[j+5]&3);//ACDEF
                        hash[hv] = big_insert_hash(hash[hv],j);
                }
                for (j = i+1; j < si->numseq;j++){


                        //min =  (si->sl[i] > si->sl[j]) ?si->sl[j] :si->sl[i];
                        dm[i][j] = dna_distance_calculation(hash,si->s[j],si->sl[j],si->sl[j]+si->sl[i],60);
                        dm[i][j] /= (si->sl[i] > si->sl[j]) ?si->sl[j] :si->sl[i];
                        dm[j][i] = dm[i][j];
                        //fprintf(stderr,"\r%8.0f percent done",(float)a /(float)b * 100);
                        a++;
                }

                for (j = 1024;j--;){
                        if (hash[j]){
                                big_remove_nodes(hash[j]);
                                hash[j] = 0;
                        }
                }
        }
        return dm;
}

float dna_distance_calculation(struct bignode* hash[],int* p,int seqlen,int diagonals,float mode)
{

        struct bignode* node_p;
        float out = 0.0;
        unsigned int* tmp = 0;
        unsigned int* d = 0;
        int i,j;
        unsigned int hv;

        d = malloc(sizeof(int)*diagonals);
        for (i = 0;i < diagonals;i++){
                d[i] = 0;
        }
        for (i = seqlen-5;i--;){

                hv = ((p[i]&3)<<8) + ((p[i+1]&3)<<6) + ((p[i+2]&3)<<4)  + ((p[i+3]&3)<<2) + (p[i+4]&3);//ABCDE
                if (hash[hv]){
                        node_p = hash[hv];
                        while(node_p){
                                tmp = node_p->pos;
                                for(j = 0;j < node_p->num;j++){
                                        d[tmp[j]]++;
                                }
                                node_p = node_p->next;
                        }
                }


                hv = ((p[i]&3)<<8) + ((p[i+1]&3)<<6) + ((p[i+2]&3)<<4)  + ((p[i+3]&3)<<2) + (p[i+5]&3);//ABCDF
                if (hash[hv]){
                        node_p = hash[hv];
                        while(node_p){
                                tmp = node_p->pos;
                                for(j = 0;j < node_p->num;j++){
                                        d[tmp[j]]++;
                                }
                                node_p = node_p->next;
                        }
                }
                hv = ((p[i]&3)<<8) + ((p[i+1]&3)<<6) + ((p[i+2]&3)<<4)  + ((p[i+4]&3)<<2) + (p[i+5]&3);//ABCEF
                if (hash[hv]){
                        node_p = hash[hv];
                        while(node_p){
                                tmp = node_p->pos;
                                for(j = 0;j < node_p->num;j++){
                                        d[tmp[j]]++;
                                }
                                node_p = node_p->next;
                        }
                }
                hv = ((p[i]&3)<<8) + ((p[i+1]&3)<<6) + ((p[i+3]&3)<<4)  + ((p[i+4]&3)<<2) + (p[i+5]&3);//ABDEF
                if (hash[hv]){
                        node_p = hash[hv];
                        while(node_p){
                                tmp = node_p->pos;
                                for(j = 0;j < node_p->num;j++){
                                        d[tmp[j]]++;
                                }
                                node_p = node_p->next;
                        }
                }
                hv = ((p[i]&3)<<8) + ((p[i+2]&3)<<6) + ((p[i+3]&3)<<4) + ((p[i+4]&3)<<2) + (p[i+5]&3);//ACDEF
                if (hash[hv]){
                        node_p = hash[hv];
                        while(node_p){
                                tmp = node_p->pos;
                                for(j = 0;j < node_p->num;j++){
                                        d[tmp[j]]++;
                                }
                                node_p = node_p->next;
                        }
                }

                d++;
        }
        //exit(0);
        d -= (seqlen-5);

        for (i = diagonals;i--;){
                //d[i] /= minlen;

                //printf("%d ",d[i]);

                if(d[i] > mode){
                        //fprintf(stderr,"%f	%d\n",d[i]/ minlen,d[i]);
                        out += d[i];
                }
        }
        free(d);
        return out;
}








struct alignment* make_seq(struct alignment* aln,int a,int b,int* path)
{
        int c;
        int i;
        int posa = 0;
        int posb = 0;

        int* gap_a = 0;
        int* gap_b = 0;

        gap_a = malloc ((path[0]+1)*sizeof(int));
        gap_b = malloc ((path[0]+1)*sizeof(int));

        for (i = path[0]+1;i--;){
                gap_a[i] = 0;
                gap_b[i] = 0;
        }
        c = 1;
        while(path[c] != 3){
                if (!path[c]){
                        posa++;
                        posb++;
                }
                if (path[c] & 1){
                        gap_a[posa] += 1;
                        posb++;
                }
                if (path[c] & 2){
                        gap_b[posb] += 1;
                        posa++;
                }
                c++;
        }


        for (i = aln->nsip[a];i--;){
                update_gaps(aln->sl[aln->sip[a][i]],aln->s[aln->sip[a][i]],path[0],gap_a);
        }

        for (i = aln->nsip[b];i--;){
                update_gaps(aln->sl[aln->sip[b][i]],aln->s[aln->sip[b][i]],path[0],gap_b);
        }
        free(gap_a);
        free(gap_b);
        free(path);
        return aln;
}


void update_gaps(int old_len,int*gis,int new_len,int *newgaps)
{
        unsigned int i,j;
        int add = 0;
        int rel_pos = 0;
        for (i = 0; i <= old_len;i++){
                add = 0;
                for (j = rel_pos;j <= rel_pos + gis[i];j++){
                        if (newgaps[j] != 0){
                                add += newgaps[j];
                        }
                }
                rel_pos += gis[i]+1;
                gis[i] += add;
        }
}


int** dna_alignment(struct alignment* aln,int* tree,float**submatrix, int** map)
{
        struct hirsch_mem* hm = 0;
        int i,j,g,a,b,c;
        int len_a;
        int len_b;
        float** profile = 0;


        profile = malloc(sizeof(float*)*aln->numprofiles);
        for ( i = 0;i< aln->numprofiles;i++){
                profile[i] = 0;
        }

        map = malloc(sizeof(int*)*aln->numprofiles);
        for ( i = 0;i < aln->numprofiles;i++){
                map[i] = 0;
        }

        hm = hirsch_mem_alloc(hm,1024);
        //fprintf(stderr,"\nAlignment:\n");
        for (i = 0; i < (aln->numseq-1);i++){
                a = tree[i*3];
                b = tree[i*3+1];
                c = tree[i*3+2];
                //fprintf(stderr,"\r%8.0f percent done",(float)(i) /(float)aln->numseq * 100);
                //fprintf(stderr,"Aligning:%d %d->%d	done:%0.2f\n",a,b,c,((float)(i+1)/(float)aln->numseq)*100);
                len_a = aln->sl[a];
                len_b = aln->sl[b];


                g = (len_a > len_b)? len_a:len_b;
                map[c] = malloc(sizeof(int) * (g+2));
                if(g > hm->size){
                        hm = hirsch_mem_realloc(hm,g);
                }

                for (j = 0; j < (g+2);j++){
                        map[c][j] = -1;
                }

                if (a < aln->numseq){
                        profile[a] = dna_make_profile(profile[a],aln->s[a],len_a,submatrix);
                }
                if (b < aln->numseq){
                        profile[b] = dna_make_profile(profile[b],aln->s[b],len_b,submatrix);
                }


                dna_set_gap_penalties(profile[a],len_a,aln->nsip[b]);
                dna_set_gap_penalties(profile[b],len_b,aln->nsip[a]);

                hm->starta = 0;
                hm->startb = 0;
                hm->enda = len_a;
                hm->endb = len_b;
                hm->len_a = len_a;
                hm->len_b = len_b;

                hm->f[0].a = 0.0;
                hm->f[0].ga =  -INFINITY;
                hm->f[0].gb = -INFINITY;
                hm->b[0].a = 0.0;
                hm->b[0].ga =  -INFINITY;
                hm->b[0].gb =  -INFINITY;
                //fprintf(stderr,"LENA:%d	LENB:%d	numseq:%d\n",len_a,len_b,aln->numseq);
                if(a < aln->numseq){
                        if(b < aln->numseq){
                                map[c] = hirsch_dna_ss_dyn(submatrix,aln->s[a],aln->s[b],hm,map[c]);
                        }else{
                                hm->enda = len_b;
                                hm->endb = len_a;
                                hm->len_a = len_b;
                                hm->len_b = len_a;
                                map[c] = hirsch_dna_ps_dyn(profile[b],aln->s[a],hm,map[c],aln->nsip[b]);
                                map[c] = mirror_hirsch_path(map[c],len_a,len_b);
                        }
                }else{
                        if(b < aln->numseq){
                                map[c] = hirsch_dna_ps_dyn(profile[a],aln->s[b],hm,map[c],aln->nsip[a]);
                        }else{
                                if(len_a < len_b){
                                        map[c] = hirsch_dna_pp_dyn(profile[a],profile[b],hm,map[c]);
                                }else{
                                        hm->enda = len_b;
                                        hm->endb = len_a;
                                        hm->len_a = len_b;
                                        hm->len_b = len_a;
                                        map[c] = hirsch_dna_pp_dyn(profile[b],profile[a],hm,map[c]);
                                        map[c] = mirror_hirsch_path(map[c],len_a,len_b);
                                }
                        }
                }


                map[c] = add_gap_info_to_hirsch_path(map[c],len_a,len_b);
                //fprintf(stderr,"%p\n",map[c]);
                if(i != aln->numseq-2){
                        profile[c] = malloc(sizeof(float)*22*(map[c][0]+2));
                        profile[c] = dna_update(profile[a],profile[b],profile[c],map[c],aln->nsip[a],aln->nsip[b]);
                }

                aln->sl[c] = map[c][0];

                aln->nsip[c] = aln->nsip[a] + aln->nsip[b];
                aln->sip[c] = malloc(sizeof(int)*(aln->nsip[a] + aln->nsip[b]));
                g =0;
                for (j = aln->nsip[a];j--;){
                        aln->sip[c][g] = aln->sip[a][j];
                        g++;
                }
                for (j = aln->nsip[b];j--;){
                        aln->sip[c][g] = aln->sip[b][j];
                        g++;
                }

                free(profile[a]);
                free(profile[b]);
                //	fprintf(stderr,"%p\n",map[c]);
        }

        ///fprintf(stderr,"\r%8.0f percent done\n",100.0);
        //free(profile[numprofiles-1]);
        free(profile);
        hirsch_mem_free(hm);
        return map;
}


float* dna_update(const float* profa, const float* profb, float* newp,int* path,int sipa,int sipb)
{
        int i,j,c;

        for (i = 22; i--;){
                newp[i] = profa[i] + profb[i];
        }

        profa += 22;
        profb += 22;
        newp += 22;


        c = 1;
        while(path[c] != 3){
                //Idea: limit the 'virtual' number of residues of one type to x.
                // i.e. only allow a maximum of 10 alanines to be registered in each column
                // the penalty for aligning a 'G' to this column will stay stable even when many (>10) alanines are present.
                // the difference in score between the 'correct' (all alanine) and incorrect (alanines + glycine) will not increase
                // with the number of sequences. -> see Durbin pp 140

                if (!path[c]){
                        //fprintf(stderr,"Align	%d\n",c);
                        for (i = 22; i--;){
                                newp[i] = profa[i] + profb[i];
                        }


                        profa += 22;
                        profb += 22;
                }
                if (path[c] & 1){
                        //fprintf(stderr,"Gap_A:%d\n",c);
                        //printf("open:%d	ext:%d	%d	%d\n",si->nsip[a] * gpo,si->nsip[a] * gpe,si->nsip[a] * profb[41],si->nsip[a] * profb[46]);
                        for (i = 22; i--;){
                                newp[i] = profb[i];
                        }
                        profb += 22;
                        if(!(path[c]&20)){
                                if(path[c]&32){
                                        newp[7] += sipa;//1;
                                        i = TGPE*sipa;
                                }else{
                                        newp[6] += sipa;//1;
                                        i = GPE*sipa;
                                }

                                for (j = 11; j < 16;j++){
                                        newp[j] -=i;
                                }
                        }else{
                                if (path[c] & 16){
                                        //			fprintf(stderr,"close_open");
                                        if(path[c]&32){
                                                newp[7] += sipa;//1;
                                                i = TGPE*sipa;
                                                newp[5] += sipa;//1;
                                                i += GPO*sipa;
                                        }else{
                                                newp[5] += sipa;//1;
                                                i = GPO*sipa;
                                        }

                                        for (j = 11; j < 16;j++){
                                                newp[j] -=i;
                                        }
                                }
                                if (path[c] & 4){
                                        //			fprintf(stderr,"Gap_open");
                                        if(path[c]&32){
                                                newp[7] += sipa;//1;
                                                i = TGPE*sipa;
                                                newp[5] += sipa;//1;
                                                i += GPO*sipa;
                                        }else{
                                                newp[5] += sipa;//1;
                                                i = GPO*sipa;
                                        }
                                        for (j = 11; j < 16; j++){
                                                newp[j] -=i;
                                        }
                                }
                        }


                }
                if (path[c] & 2){
                        //fprintf(stderr,"Gap_B:%d\n",c);
                        //printf("open:%d	ext:%d	%d	%d\n",si->nsip[b] * gpo,si->nsip[b] * gpe,profa[26],profa[27]);
                        for (i = 22; i--;){
                                newp[i] = profa[i];
                        }
                        profa+=22;
                        if(!(path[c]&20)){
                                if(path[c]&32){
                                        newp[7] += sipb;//1;
                                        i = TGPE*sipb;
                                }else{
                                        newp[6] += sipb;//1;
                                        i = GPE*sipb;
                                }
                                for (j = 11; j < 16;j++){
                                        newp[j] -=i;
                                }
                        }else{
                                if (path[c] & 16){
                                        //			fprintf(stderr,"close_open");
                                        if(path[c]&32){
                                                newp[7] += sipb;//1;
                                                i =  TGPE*sipb;
                                                newp[5] += sipb;//1;
                                                i +=  GPO*sipb;
                                        }else{
                                                newp[5] += sipb;//1;
                                                i =  GPO*sipb;
                                        }
                                        for (j = 11; j < 16;j++){
                                                newp[j] -=i;
                                        }
                                }
                                if (path[c] & 4){
                                        //			fprintf(stderr,"Gap_open");
                                        if(path[c]&32){
                                                newp[7] += sipb;//1;
                                                i = TGPE*sipb;
                                                newp[5] += sipb;//1;
                                                i += GPO*sipb;
                                        }else{
                                                newp[5] += sipb;//1;
                                                i = GPO*sipb;
                                        }

                                        for (j = 11; j < 16;j++){
                                                newp[j] -=i;
                                        }
                                }
                        }

                }
                newp += 22;
                c++;
        }
        for (i = 22; i--;){
                newp[i] =  profa[i] + profb[i];
        }
        newp -= (path[0]+1) *22;
        return newp;
}

int* mirror_hirsch_path(int* hirsch_path,int len_a,int len_b)
{
        int* np = 0;

        int i;
        np =malloc(sizeof(int)*(len_a+2));
        for(i =0; i < len_a+2;i++){
                np[i] = -1;
        }

        for(i = 1; i <= len_b;i++){
                if(hirsch_path[i] != -1){
                        np[hirsch_path[i]] = i;
                }
        }

        free(hirsch_path);
        return np;
}

int* add_gap_info_to_hirsch_path(int* hirsch_path,int len_a,int len_b)
{
        int i,j;
        int a = 0;
        int b = 0;

        int* np = 0;
        np =malloc(sizeof(int)*(len_a+len_b+2));
        for(i =0; i < len_a+len_b+2;i++){
                np[i] = 0;
        }

        j = 1;
        b = -1;
        if(hirsch_path[1] == -1){
                np[j] = 2;
                j++;
        }else{
                if(hirsch_path[1] != 1){
                        for ( a = 0;a < hirsch_path[1] -1;a++){
                                np[j] = 1;
                                j++;
                        }
                        np[j] = 0;
                        j++;
                }else{
                        np[j] = 0;
                        j++;
                }
        }
        b = hirsch_path[1];

        /*for ( i= 0;i <= len_a;i++){
          fprintf(stderr,"%d,",hirsch_path[i]);
          }
          fprintf(stderr,"\n");*/

        for(i = 2; i <= len_a;i++){

                if(hirsch_path[i] == -1){
                        np[j] = 2;
                        j++;
                }else{
                        if(hirsch_path[i]-1 != b && b != -1){
                                for ( a = 0;a < hirsch_path[i] - b-1;a++){
                                        np[j] = 1;
                                        j++;
                                }
                                np[j] = 0;
                                j++;
                        }else{
                                np[j] = 0;
                                j++;
                        }
                }
                b = hirsch_path[i];
        }





        if(hirsch_path[len_a] < len_b && hirsch_path[len_a] != -1){
                //	fprintf(stderr,"WARNING:%d	%d\n",hirsch_path[len_a],len_b);
                for ( a = 0;a < len_b - hirsch_path[len_a];a++){
                        np[j] = 1;
                        j++;
                }
        }
        np[0] = j-1;
        np[j] = 3;
        np = realloc(np,sizeof(int)* (np[0]+2));
        //for ( i= 0;i <= np[0];i++){
        //	fprintf(stderr,"%d,",np[i]);
        //}
        //fprintf(stderr,"\n");

        free(hirsch_path);

        //add gap info..
        i = 2;
        while(np[i] != 3){
                if ((np[i-1] &3) && !(np[i] & 3)){
                        if(np[i-1] & 8){
                                np[i-1] += 8;
                        }else{
                                np[i-1] |= 16;
                        }
                }else if (!(np[i-1] & 3) &&(np[i] &3)){
                        np[i] |= 4;
                }else if ((np[i-1] & 1) && (np[i] & 1)){
                        np[i] |= 8;
                }else if ((np[i-1] & 2) && (np[i] & 2)){
                        np[i] |= 8;
                }
                i++;
        }
        //add terminal gap...
        i = 1;
        while(np[i] != 0){
                np[i] |= 32;
                i++;
        }
        j = i;
        i = np[0];
        while(np[i] != 0){
                np[i] |= 32;
                i--;
        }
        //for ( i= 0;i <= np[0];i++){
        //	fprintf(stderr,"%d,",np[i]);
        //}
        //fprintf(stderr,"\n");
        return np;
}


int* hirsch_dna_ss_dyn(float**subm, const int* seq1,const int* seq2,struct hirsch_mem* hm, int* hirsch_path)
{
        int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
        float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
        int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};

        if(hm->starta  >= hm->enda){
                return hirsch_path;
        }
        if(hm->startb  >= hm->endb){
                return hirsch_path;
        }


        hm->enda = mid;

        //fprintf(stderr,"Forward:%d-%d	%d-%d\n",hm->starta,hm->enda,hm->startb,hm->endb);
        hm->f = foward_hirsch_dna_ss_dyn(subm,seq1,seq2,hm);

        hm->starta = mid;
        hm->enda = old_cor[1];
        //fprintf(stderr,"Backward:%d-%d	%d-%d\n",hm->starta,hm->enda,hm->startb,hm->endb);
        hm->b = backward_hirsch_dna_ss_dyn(subm,seq1,seq2,hm);


        hirsch_path = hirsch_align_two_dna_ss_vector(subm,seq1,seq2,hm,hirsch_path,input_states,old_cor);
        return hirsch_path;
}

int* hirsch_align_two_dna_ss_vector(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[])
{
        struct states* f = hm->f;
        struct states* b = hm->b;
        int i,c;
        int transition = -1;


        //code:
        // a -> a = 1
        // a -> ga = 2
        // a -> gb = 3
        // ga ->ga = 4
        // ga -> a = 5
        //gb->gb = 6;
        //gb->a = 7;

        //int max = -INFTY;
        float max = -INFINITY;
        float middle =  (hm->endb - hm->startb)/2 + hm->startb;
        float sub = 0.0;

        i = hm->startb;
        c = -1;
        for(i = hm->startb; i < hm->endb;i++){
                sub = fabsf(middle -i);
                sub /= 1000;
                //	fprintf(stderr,"%d-%d	%f\n",hm->startb,hm->endb,sub);
                if(f[i].a+b[i].a-sub > max){
                        max = f[i].a+b[i].a-sub;
                        //		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
                        transition = 1;
                        c = i;
                }
                if(f[i].a+b[i].ga-GPO-sub > max){
                        max = f[i].a+b[i].ga-GPO-sub;
                        //		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
                        transition = 2;
                        c = i;
                }
                if(f[i].a+b[i].gb -GPO-sub > max){
                        max = f[i].a+b[i].gb - GPO-sub;
                        //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
                        transition = 3;
                        c = i;
                }
                if(f[i].ga+b[i].a - GPO-sub > max){
                        max = f[i].ga+b[i].a - GPO-sub;
                        //		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
                        transition = 5;
                        c = i;
                }


                if(hm->startb == 0){
                        if(f[i].gb+b[i].gb - TGPE-sub > max){
                                max = f[i].gb+b[i].gb -TGPE-sub;
                                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                                transition = 6;
                                c = i;
                        }
                }else{
                        if(f[i].gb+b[i].gb - GPE -sub> max){
                                max = f[i].gb+b[i].gb - GPE-sub;
                                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                                transition = 6;
                                c = i;
                        }
                }
                if(f[i].gb+b[i].a - GPO-sub > max){
                        max = f[i].gb+b[i].a - GPO-sub;
                        //		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
                        transition = 7;
                        c = i;
                }
        }
        i = hm->endb;
        sub = fabsf(middle -i);
        sub /= 1000;

        if(f[i].a+b[i].gb-GPO-sub > max){
                max = f[i].a+b[i].gb - GPO-sub;
                //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
                transition = 3;
                c = i;
        }
        if(hm->endb == hm->len_b){
                if(f[i].gb+b[i].gb -TGPE-sub > max){
                        max = f[i].gb+b[i].gb - TGPE-sub;
                        //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                        transition = 6;
                        c = i;
                }
        }else{
                if(f[i].gb+b[i].gb - GPE-sub > max){
                        max = f[i].gb+b[i].gb - GPE-sub;
                        //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                        transition = 6;
                        c = i;
                }
        }


        //fprintf(stderr,"Transition:%d	at:%d\n",transition,c);

        //j = hirsch_path[0];
        switch(transition){
        case 1: //a -> a = 1

                hirsch_path[old_cor[4]] = c;
                hirsch_path[old_cor[4]+1] = c+1;

                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -INFINITY;
                hm->b[0].gb = -INFINITY;
                //		fprintf(stderr,"Using this for start:%d	%d	%d\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -INFINITY;
                hm->f[0].gb = -INFINITY;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
                break;
        case 2:// a -> ga = 2

                hirsch_path[old_cor[4]] = c;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -INFINITY;
                hm->b[0].gb = -INFINITY;


                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);


                //backward:
                hm->starta = old_cor[4];
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = -INFINITY;
                hm->f[0].ga = 0.0;
                hm->f[0].gb = -INFINITY;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
                break;
        case 3:// a -> gb = 3

                hirsch_path[old_cor[4]] = c;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -INFINITY;
                hm->b[0].gb = -INFINITY;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c;
                hm->endb = old_cor[3];
                hm->f[0].a = -INFINITY;
                hm->f[0].ga = -INFINITY;
                hm->f[0].gb = 0.0;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
                break;
        case 5://ga -> a = 5
                hirsch_path[old_cor[4]+1] = c+1;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -INFINITY;
                hm->b[0].ga = 0.0;
                hm->b[0].gb = -INFINITY;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4];

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -INFINITY;
                hm->f[0].gb = -INFINITY;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
                break;
        case 6://gb->gb = 6;

                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -INFINITY;
                hm->b[0].ga = -INFINITY;
                hm->b[0].gb = 0.0;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;
                hm->startb = old_cor[2];
                hm->endb = c;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c;
                hm->endb = old_cor[3];
                hm->f[0].a = -INFINITY;
                hm->f[0].ga = -INFINITY;
                hm->f[0].gb = 0.0;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
                break;
        case 7://gb->a = 7;

                hirsch_path[old_cor[4]+1] = c+1;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -INFINITY;
                hm->b[0].ga = -INFINITY;
                hm->b[0].gb = 0.0;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;
                hm->startb = old_cor[2];
                hm->endb = c;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -INFINITY;
                hm->f[0].gb = -INFINITY;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
                break;
        }

        return hirsch_path;
}


struct states* foward_hirsch_dna_ss_dyn(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm)
{
        struct states* s = hm->f;
        float *subp = 0;
        const int starta = hm->starta;
        const int enda = hm->enda;
        const int startb = hm->startb;
        const int endb = hm->endb;

        register float pa = 0;
        register float pga = 0;
        register float pgb = 0;
        register float ca = 0;
        register int i = 0;
        register int j = 0;


        s[startb].a = s[0].a;
        s[startb].ga = s[0].ga;
        s[startb].gb = s[0].gb;
        if(startb == 0){
                for (j = startb+1; j < endb;j++){

                        s[j].a = -INFINITY;
                        s[j].ga = KALIGN_MAX(s[j-1].ga,s[j-1].a)-TGPE;

                        s[j].gb = -INFINITY;
                }
        }else{

                for (j = startb+1; j < endb;j++){
                        s[j].a = -INFINITY;
                        s[j].ga = KALIGN_MAX(s[j-1].ga - GPE,s[j-1].a-GPO);
                        s[j].gb = -INFINITY;
                }
        }
        s[endb].a = -INFINITY;
        s[endb].ga = -INFINITY;
        s[endb].gb = -INFINITY;
        seq2--;

        for (i = starta;i < enda;i++){
                subp = subm[seq1[i]];

                pa = s[startb].a;
                pga = s[startb].ga;
                pgb = s[startb].gb;
                s[startb].a = -INFINITY;
                s[startb].ga = -INFINITY;
                if(startb == 0){
                        s[startb].gb = KALIGN_MAX(pgb,pa) - TGPE;
                }else{
                        s[startb].gb = KALIGN_MAX(pgb - GPE,pa - GPO);
                }
                for (j = startb+1; j < endb;j++){
                        ca = s[j].a;
                        pa = KALIGN_MAX3(pa,pga-GPO,pgb-GPO);
                        pa += subp[seq2[j]];

                        s[j].a = pa;

                        pga = s[j].ga;

                        s[j].ga = KALIGN_MAX(s[j-1].ga-GPE,s[j-1].a-GPO);

                        pgb = s[j].gb;

                        s[j].gb = KALIGN_MAX(pgb-GPE ,ca-GPO);

                        pa = ca;
                }
                ca = s[j].a;
                pa = KALIGN_MAX3(pa,pga-GPO,pgb-GPO);
                pa += subp[seq2[j]];

                s[j].a = pa;

                s[j].ga = -INFINITY;//MAX(s[j-1].ga-gpe,s[j-1].a-gpo);
                if (endb != hm->len_b){
                        s[j].gb = KALIGN_MAX(s[j].gb-GPE ,ca-GPO);
                }else{
                        s[j].gb = KALIGN_MAX(s[j].gb,ca)-TGPE;
                }
        }
        return s;
}

struct states* backward_hirsch_dna_ss_dyn(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm)
{

        struct states* s = hm->b;
        float *subp = 0;
        const int starta = hm->starta;
        const int enda = hm->enda;
        const int startb = hm->startb;
        const int endb = hm->endb;
        register float pa = 0;
        register float pga = 0;
        register float pgb = 0;
        register float ca = 0;
        register int i = 0;
        register int j = 0;

        s[endb].a = s[0].a ;
        s[endb].ga = s[0].ga;
        s[endb].gb = s[0].gb;


        //init of first row;

        //j = endb-startb;
        if(endb == hm->len_b){
                for(j = endb-1;j > startb;j--){
                        s[j].a = -INFINITY;
                        s[j].ga = KALIGN_MAX(s[j+1].ga,s[j+1].a)-TGPE;
                        s[j].gb = -INFINITY;
                }
        }else{
                for(j = endb-1;j > startb;j--){
                        s[j].a = -INFINITY;
                        s[j].ga = KALIGN_MAX(s[j+1].ga-GPE,s[j+1].a-GPO);
                        s[j].gb = -INFINITY;
                }
        }


        s[startb].a = -INFINITY;
        s[startb].ga = -INFINITY;
        s[startb].gb = -INFINITY;

        i = enda-starta;
        seq1+= starta;
        while(i--){
                subp = subm[seq1[i]];
                pa = s[endb].a;
                pga = s[endb].ga;
                pgb = s[endb].gb;
                s[endb].a = -INFINITY;
                s[endb].ga = -INFINITY;

                if(endb == hm->len_b){
                        s[endb].gb = KALIGN_MAX(pgb,pa)-TGPE;
                }else{
                        s[endb].gb = KALIGN_MAX(pgb-GPE,pa-GPO);
                }

                for(j = endb-1;j > startb;j--){

                        ca = s[j].a;
                        pa = KALIGN_MAX3(pa,pga - GPO,pgb-GPO);

                        pa += subp[seq2[j]];

                        s[j].a = pa;

                        pga = s[j].ga;

                        s[j].ga = KALIGN_MAX(s[j+1].ga-GPE,s[j+1].a-GPO);

                        pgb = s[j].gb;

                        s[j].gb = KALIGN_MAX(pgb-GPE,ca-GPO);

                        pa = ca;
                }
                ca = s[j].a;

                pa = KALIGN_MAX3(pa,pga - GPO,pgb-GPO);

                pa += subp[seq2[j]];

                s[j].a = pa;

                s[j].ga = -INFINITY;//MAX(s[j+1].ga-gpe,s[j+1].a-gpo);

                if(startb){
                        s[j].gb = KALIGN_MAX(s[j].gb-GPE,ca-GPO);
                }else{
                        s[j].gb = KALIGN_MAX(s[j].gb,ca)-TGPE;
                }
        }
        return s;
}




int* hirsch_dna_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm, int* hirsch_path,int sip)
{
        int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
        float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
        int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};


        if(hm->starta  >= hm->enda){
                return hirsch_path;
        }
        if(hm->startb  >= hm->endb){
                return hirsch_path;
        }

        hm->enda = mid;
        hm->f = foward_hirsch_dna_ps_dyn(prof1,seq2,hm,sip);

        /*int i;
          fprintf(stderr,"FOWARD\n");
          for (i = hm->startb; i <= hm->endb;i++){
          fprintf(stderr,"%d	%d	%d\n",hm->f[i].a,hm->f[i].ga,hm->f[i].gb);
          }*/

        hm->starta = mid;
        hm->enda = old_cor[1];
        hm->b = backward_hirsch_dna_ps_dyn(prof1,seq2,hm,sip);

        /*fprintf(stderr,"BaCKWARD\n");
          for (i = hm->startb; i <= hm->endb;i++){
          fprintf(stderr,"%d	%d	%d\n",hm->b[i].a,hm->b[i].ga,hm->b[i].gb);
          }*/

        hirsch_path = hirsch_align_two_dna_ps_vector(prof1,seq2,hm,hirsch_path,input_states,old_cor,sip);
        return hirsch_path;
}



int* hirsch_align_two_dna_ps_vector(const float* prof1,const int* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],int sip)
{
        struct states* f = hm->f;
        struct states* b = hm->b;
        int i,c;
        int transition = -1;
        const int open = GPO * sip;
        //code:
        // a -> a = 1
        // a -> ga = 2
        // a -> gb = 3
        // ga ->ga = 4
        // ga -> a = 5
        //gb->gb = 6;
        //gb->a = 7;
        //int max = -INFTY;
        float max = -INFINITY;
        float middle =  (hm->endb - hm->startb)/2 + hm->startb;
        float sub = 0.0;

        prof1+= (22 * (old_cor[4]+1));
        i = hm->startb;
        c = -1;
        for(i = hm->startb; i < hm->endb;i++){
                sub = fabsf(middle -i);
                sub /= 1000;
                if(f[i].a+b[i].a-sub> max){
                        max = f[i].a+b[i].a-sub;
                        //		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
                        transition = 1;
                        c = i;
                }
                if(f[i].a+b[i].ga-open-sub > max){
                        max = f[i].a+b[i].ga-open-sub;
                        //		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
                        transition = 2;
                        c = i;
                }
                if(f[i].a+b[i].gb+prof1[8]-sub > max){
                        max = f[i].a+b[i].gb+prof1[8]-sub;
                        //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
                        transition = 3;
                        c = i;
                }
                if(f[i].ga+b[i].a-open-sub > max){
                        max = f[i].ga+b[i].a-open-sub;
                        //		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
                        transition = 5;
                        c = i;
                }


                if(hm->startb == 0){
                        if(f[i].gb+b[i].gb+prof1[10]-sub > max){
                                max = f[i].gb+b[i].gb+prof1[10]-sub;
                                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                                transition = 6;
                                c = i;
                        }
                }else{
                        if(f[i].gb+b[i].gb+prof1[9]-sub > max){
                                max = f[i].gb+b[i].gb+prof1[9]-sub;
                                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                                transition = 6;
                                c = i;
                        }
                }
                if(f[i].gb+b[i].a+prof1[8-22]-sub > max){
                        max = f[i].gb+b[i].a+prof1[8-22]-sub;
                        //		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
                        transition = 7;
                        c = i;
                }
        }
        i = hm->endb;
        sub = fabsf(middle -i);
        sub /= 1000;
        if(f[i].a+b[i].gb+prof1[8]-sub > max){
                max = f[i].a+b[i].gb+prof1[8]-sub;
                //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
                transition = 3;
                c = i;
        }
        if(hm->endb == hm->len_b){
                if(f[i].gb+b[i].gb+prof1[10]-sub > max){
                        max = f[i].gb+b[i].gb+prof1[10]-sub;
                        //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                        transition = 6;
                        c = i;
                }
        }else{
                if(f[i].gb+b[i].gb+prof1[9]-sub > max){
                        max = f[i].gb+b[i].gb+prof1[0]-sub;
                        //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                        transition = 6;
                        c = i;
                }
        }



        prof1-= (22 * (old_cor[4]+1));

        //fprintf(stderr,"Transition:%d	at:%d\n",transition,c);

        //j = hirsch_path[0];
        switch(transition){
        case 1: //a -> a = 1

                hirsch_path[old_cor[4]] = c;
                hirsch_path[old_cor[4]+1] = c+1;

                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -INFINITY;
                hm->b[0].gb = -INFINITY;
                //		fprintf(stderr,"Using this for start:%d	%d	%d\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -INFINITY;
                hm->f[0].gb = -INFINITY;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
                break;
        case 2:// a -> ga = 2

                hirsch_path[old_cor[4]] = c;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -INFINITY;
                hm->b[0].gb = -INFINITY;


                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

                //backward:
                hm->starta = old_cor[4];
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = -INFINITY;
                hm->f[0].ga = 0.0;
                hm->f[0].gb = -INFINITY;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
                break;
        case 3:// a -> gb = 3

                hirsch_path[old_cor[4]] = c;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -INFINITY;
                hm->b[0].gb = -INFINITY;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c;
                hm->endb = old_cor[3];
                hm->f[0].a = -INFINITY;
                hm->f[0].ga = -INFINITY;
                hm->f[0].gb = 0.0;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
                break;
        case 5://ga -> a = 5
                hirsch_path[old_cor[4]+1] = c+1;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -INFINITY;
                hm->b[0].ga = 0.0;
                hm->b[0].gb = -INFINITY;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4];

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -INFINITY;
                hm->f[0].gb = -INFINITY;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
                break;
        case 6://gb->gb = 6;

                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -INFINITY;
                hm->b[0].ga = -INFINITY;
                hm->b[0].gb = 0.0;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;
                hm->startb = old_cor[2];
                hm->endb = c;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);


                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c;
                hm->endb = old_cor[3];
                hm->f[0].a = -INFINITY;
                hm->f[0].ga = -INFINITY;
                hm->f[0].gb = 0.0;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
                break;
        case 7://gb->a = 7;

                hirsch_path[old_cor[4]+1] = c+1;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -INFINITY;
                hm->b[0].ga = -INFINITY;
                hm->b[0].gb = 0.0;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;
                hm->startb = old_cor[2];
                hm->endb = c;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -INFINITY;
                hm->f[0].gb = -INFINITY;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
                break;
        }

        return hirsch_path;
}

struct states* foward_hirsch_dna_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm,int sip)
{
        //unsigned int freq[26];
        struct states* s = hm->f;
        const int starta = hm->starta;
        const int enda = hm->enda;
        const int startb = hm->startb;
        const int endb = hm->endb;

        register float pa = 0;
        register float pga = 0;
        register float pgb = 0;
        register float ca = 0;
        register int i = 0;
        register int j = 0;

        const float open = GPO * sip;
        const float ext = GPE *sip;
        const float text = TGPE * sip;



        prof1 += (starta) * 22;
        s[startb].a = s[0].a;
        s[startb].ga = s[0].ga;
        s[startb].gb = s[0].gb;
        if(startb == 0){
                for (j = startb+1; j < endb;j++){
                        s[j].a = -INFINITY;
                        s[j].ga = KALIGN_MAX(s[j-1].ga,s[j-1].a) - text;
                        s[j].gb = -INFINITY;
                }
        }else{
                for (j = startb+1; j < endb;j++){
                        s[j].a = -INFINITY;
                        s[j].ga = KALIGN_MAX(s[j-1].ga-ext,s[j-1].a-open);
                        s[j].gb = -INFINITY;
                }
        }


        s[endb].a = -INFINITY;
        s[endb].ga = -INFINITY;
        s[endb].gb = -INFINITY;
        seq2--;

        for (i = starta;i < enda;i++){
                prof1 += 22;
                pa = s[startb].a;
                pga = s[startb].ga;
                pgb = s[startb].gb;
                s[startb].a = -INFINITY;
                s[startb].ga = -INFINITY;
                if(startb == 0){
                        s[startb].gb = KALIGN_MAX(pgb,pa)+prof1[10];
                }else{
                        s[startb].gb = KALIGN_MAX(pgb+prof1[9],pa+prof1[8]);
                }
                for (j = startb+1; j < endb;j++){
                        ca = s[j].a;
                        pa = KALIGN_MAX3(pa,pga -open,pgb + prof1[-14]);
                        pa += prof1[11 + seq2[j]];


                        s[j].a = pa;

                        pga = s[j].ga;

                        s[j].ga = KALIGN_MAX(s[j-1].ga-ext,s[j-1].a-open);

                        pgb = s[j].gb;

                        s[j].gb = KALIGN_MAX(pgb+prof1[9],ca+prof1[8]);

                        pa = ca;
                }
                ca = s[j].a;

                pa = KALIGN_MAX3(pa,pga -open,pgb + prof1[-14]);

                pa += prof1[11 + seq2[j]];


                s[j].a = pa;

                s[j].ga = -INFINITY;//MAX(s[j-1].ga-ext,s[j-1].a-open);

                if (hm->endb != hm->len_b){
                        s[j].gb = KALIGN_MAX(s[j].gb+prof1[9] ,ca+prof1[8]);
                }else{
                        s[j].gb = KALIGN_MAX(s[j].gb,ca)+ prof1[10];
                }
        }
        prof1 -= 22 * (enda);
        return s;
}

struct states* backward_hirsch_dna_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm,int sip)
{
        //unsigned int freq[26];
        struct states* s = hm->b;
        const int starta = hm->starta;
        const int enda = hm->enda;
        const int startb = hm->startb;
        const int endb = hm->endb;

        register float pa = 0;
        register float pga = 0;
        register float pgb = 0;
        register float ca = 0;
        register int i = 0;
        register int j = 0;

        const float open = GPO * sip;
        const float ext = GPE *sip;
        const float text = TGPE * sip;


        prof1 += (enda+1) * 22;

        s[endb].a = s[0].a;
        s[endb].ga = s[0].ga;
        s[endb].gb = s[0].gb;


        //init of first row;
        //j = endb-startb;
        if(endb == hm->len_b){
                for(j = endb-1;j > startb;j--){
                        s[j].a = -INFINITY;
                        s[j].ga = KALIGN_MAX(s[j+1].ga,s[j+1].a)-text;
                        s[j].gb = -INFINITY;
                }
        }else{
                for(j = endb-1;j > startb;j--){
                        s[j].a = -INFINITY;
                        s[j].ga = KALIGN_MAX(s[j+1].ga-ext,s[j+1].a-open);
                        s[j].gb = -INFINITY;
                }
        }

        s[startb].a = -INFINITY;
        s[startb].ga = -INFINITY;
        s[startb].gb = -INFINITY;

        i = enda-starta;
        while(i--){
                prof1 -= 22;

                pa = s[endb].a;
                pga = s[endb].ga;
                pgb = s[endb].gb;
                s[endb].a = -INFINITY;
                s[endb].ga = -INFINITY;

                if(endb == hm->len_b){
                        s[endb].gb = KALIGN_MAX(pgb,pa) +prof1[10];
                }else{
                        s[endb].gb = KALIGN_MAX(pgb+prof1[9],pa+prof1[8]);
                }

                for(j = endb-1;j > startb;j--){
                        ca = s[j].a;
                        pa = KALIGN_MAX3(pa,pga - open,pgb +prof1[30]);
                        pa += prof1[11 + seq2[j]];

                        s[j].a = pa;

                        pga = s[j].ga;

                        s[j].ga = KALIGN_MAX(s[j+1].ga-ext,s[j+1].a-open);

                        pgb = s[j].gb;

                        s[j].gb = KALIGN_MAX(pgb+prof1[9],ca+prof1[8]);

                        pa = ca;
                }
                ca = s[j].a;

                pa = KALIGN_MAX3(pa,pga - open,pgb +prof1[30]);
                pa += prof1[11 + seq2[j]];

                s[j].a = pa;


                s[j].ga = -INFINITY;//MAX(s[j+1].ga-ext,s[j+1].a-open);
                if(hm->startb){
                        s[j].gb = KALIGN_MAX(s[j].gb+prof1[9], ca+prof1[8]);
                }else{
                        s[j].gb = KALIGN_MAX(s[j].gb,ca)+prof1[10];
                }
        }
        return s;
}




int* hirsch_dna_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm, int* hirsch_path)
{
        int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
        float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
        int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};


        //fprintf(stderr,"starta:%d enda:%d startb:%d endb:%d mid:%d\n",hm->starta,hm->enda,hm->startb,hm->endb,mid);


        if(hm->starta  >= hm->enda){
                return hirsch_path;
        }
        if(hm->startb  >= hm->endb){
                return hirsch_path;
        }

        hm->enda = mid;
        hm->f = foward_hirsch_dna_pp_dyn(prof1,prof2,hm);
        /*int i;
          fprintf(stderr,"FOWARD\n");
          for (i = hm->startb; i <= hm->endb;i++){
          fprintf(stderr,"%d	%d	%d\n",hm->f[i].a,hm->f[i].ga,hm->f[i].gb);
          }*/

        hm->starta = mid;
        hm->enda = old_cor[1];
        hm->b = backward_hirsch_dna_pp_dyn(prof1,prof2,hm);
        /*fprintf(stderr,"BaCKWARD\n");

          for (i = hm->startb; i <= hm->endb;i++){
          fprintf(stderr,"%d	%d	%d\n",hm->b[i].a,hm->b[i].ga,hm->b[i].gb);
          }*/

        hirsch_path = hirsch_align_two_dna_pp_vector(prof1,prof2,hm,hirsch_path,input_states,old_cor);
        return hirsch_path;
}



int* hirsch_align_two_dna_pp_vector(const float* prof1,const float* prof2,struct hirsch_mem* hm,int* hirsch_path, float input_states[],int old_cor[])
{
        struct states* f = hm->f;
        struct states* b = hm->b;
        int i,c;
        int transition = -1;


        //code:
        // a -> a = 1
        // a -> ga = 2
        // a -> gb = 3
        // ga ->ga = 4
        // ga -> a = 5
        //gb->gb = 6;
        //gb->a = 7;

        //int max = -INFTY;
        float max = -INFINITY;
        float middle =  (hm->endb - hm->startb)/2 + hm->startb;
        float sub = 0.0;


        prof1+= (22 * (old_cor[4]+1));
        prof2 += (22 * (hm->startb));

        i = hm->startb;
        c = -1;
        for(i = hm->startb; i < hm->endb;i++){
                sub = fabsf(middle -i);
                sub /= 1000;
                prof2 += 22;
                if(f[i].a+b[i].a-sub > max){
                        max = f[i].a+b[i].a-sub;
                        //		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
                        transition = 1;
                        c = i;
                }
                if(f[i].a+b[i].ga+prof2[8]-sub > max){
                        max = f[i].a+b[i].ga+prof2[8]-sub;
                        //		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
                        transition = 2;
                        c = i;
                }
                if(f[i].a+b[i].gb+prof1[8]-sub > max){
                        max = f[i].a+b[i].gb+prof1[8]-sub;
                        //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
                        transition = 3;
                        c = i;
                }
                if(f[i].ga+b[i].a+prof2[-14]-sub > max){
                        max = f[i].ga+b[i].a+prof2[-14]-sub;
                        //		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
                        transition = 5;
                        c = i;
                }


                if(hm->startb == 0){
                        if(f[i].gb+b[i].gb+prof1[10]-sub > max){
                                max = f[i].gb+b[i].gb+prof1[10]-sub;
                                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                                transition = 6;
                                c = i;
                        }
                }else{
                        if(f[i].gb+b[i].gb+prof1[9]-sub > max){
                                max = f[i].gb+b[i].gb+prof1[9]-sub;
                                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                                transition = 6;
                                c = i;
                        }
                }
                if(f[i].gb+b[i].a+prof1[-14]-sub > max){
                        max = f[i].gb+b[i].a+prof1[-14]-sub;
                        //		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
                        transition = 7;
                        c = i;
                }
        }
        i = hm->endb;
        sub = fabsf(middle -i);
        sub /= 1000;
        if(f[i].a+b[i].gb+prof1[8]-sub > max){
                max = f[i].a+b[i].gb+prof1[8]-sub;
                //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
                transition = 3;
                c = i;
        }
        if(hm->endb == hm->len_b){
                if(f[i].gb+b[i].gb+prof1[10]-sub > max){
                        max = f[i].gb+b[i].gb+prof1[10]-sub;
                        //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                        transition = 6;
                        c = i;
                }
        }else{
                if(f[i].gb+b[i].gb+prof1[9]-sub > max){
                        max = f[i].gb+b[i].gb+prof1[9]-sub;
                        //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                        transition = 6;
                        c = i;
                }
        }



        prof1-= (22 * (old_cor[4]+1));
        prof2 -= (hm->endb *22);

        //fprintf(stderr,"Transition:%d	at:%d\n",transition,c);

        //j = hirsch_path[0];
        switch(transition){
        case 1: //a -> a = 1

                hirsch_path[old_cor[4]] = c;
                hirsch_path[old_cor[4]+1] = c+1;

                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -INFINITY;
                hm->b[0].gb = -INFINITY;
                //		fprintf(stderr,"Using this for start:%d	%d	%d\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -INFINITY;
                hm->f[0].gb = -INFINITY;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);
                break;
        case 2:// a -> ga = 2

                hirsch_path[old_cor[4]] = c;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -INFINITY;
                hm->b[0].gb = -INFINITY;


                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4];
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = -INFINITY;
                hm->f[0].ga = 0.0;
                hm->f[0].gb = -INFINITY;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);
                break;
        case 3:// a -> gb = 3

                hirsch_path[old_cor[4]] = c;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -INFINITY;
                hm->b[0].gb = -INFINITY;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c;
                hm->endb = old_cor[3];
                hm->f[0].a = -INFINITY;
                hm->f[0].ga = -INFINITY;
                hm->f[0].gb = 0.0;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);
                break;
        case 5://ga -> a = 5
                hirsch_path[old_cor[4]+1] = c+1;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -INFINITY;
                hm->b[0].ga = 0.0;
                hm->b[0].gb = -INFINITY;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4];

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -INFINITY;
                hm->f[0].gb = -INFINITY;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);
                break;
        case 6://gb->gb = 6;

                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -INFINITY;
                hm->b[0].ga = -INFINITY;
                hm->b[0].gb = 0.0;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;
                hm->startb = old_cor[2];
                hm->endb = c;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c;
                hm->endb = old_cor[3];
                hm->f[0].a = -INFINITY;
                hm->f[0].ga = -INFINITY;
                hm->f[0].gb = 0.0;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);
                break;
        case 7://gb->a = 7;

                hirsch_path[old_cor[4]+1] = c+1;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -INFINITY;
                hm->b[0].ga = -INFINITY;
                hm->b[0].gb = 0.0;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;
                hm->startb = old_cor[2];
                hm->endb = c;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -INFINITY;
                hm->f[0].gb = -INFINITY;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);
                break;
        }

        return hirsch_path;
}

struct states* foward_hirsch_dna_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm)
{
        struct states* s = hm->f;

        register float pa = 0;
        register float pga = 0;
        register float pgb = 0;
        register float ca = 0;
        register int i = 0;
        register int j = 0;



        prof1 += (hm->starta) * 22;
        prof2 +=  (hm->startb) * 22;
        s[hm->startb].a = s[0].a;
        s[hm->startb].ga = s[0].ga;
        s[hm->startb].gb = s[0].gb;
        if(hm->startb == 0){
                for (j = hm->startb+1; j < hm->endb;j++){
                        prof2+=22;
                        s[j].a = -INFINITY;
                        s[j].ga = KALIGN_MAX(s[j-1].ga,s[j-1].a)+prof2[10];
                        s[j].gb = -INFINITY;
                }
                prof2 += 22;
        }else{

                for (j = hm->startb+1; j < hm->endb;j++){
                        prof2+=22;
                        s[j].a = -INFINITY;
                        s[j].ga = KALIGN_MAX(s[j-1].ga+prof2[9],s[j-1].a+prof2[8]);
                        s[j].gb = -INFINITY;
                }
                prof2 += 22;
        }

        prof2 -= (hm->endb-hm->startb) * 22;

        s[hm->endb].a = -INFINITY;
        s[hm->endb].ga = -INFINITY;
        s[hm->endb].gb = -INFINITY;


        for (i = hm->starta;i < hm->enda;i++){
                prof1 += 22;

                pa = s[hm->startb].a;
                pga = s[hm->startb].ga;
                pgb = s[hm->startb].gb;
                s[hm->startb].a = -INFINITY;
                s[hm->startb].ga = -INFINITY;
                if(hm->startb == 0){
                        s[hm->startb].gb = KALIGN_MAX(pgb,pa)+ prof1[10];
                }else{
                        s[hm->startb].gb = KALIGN_MAX(pgb+prof1[9],pa+prof1[8]);
                }
                for (j = hm->startb+1; j < hm->endb;j++){
                        prof2 += 22;
                        ca = s[j].a;
                        pa = KALIGN_MAX3(pa,pga + prof2[-14],pgb + prof1[-14]);

                        prof2 += 11;

                        pa += prof1[0]*prof2[0];
                        pa += prof1[1]*prof2[1];
                        pa += prof1[2]*prof2[2];
                        pa += prof1[3]*prof2[3];
                        pa += prof1[4]*prof2[4];
                        pa += prof1[5]*prof2[5];
                        pa += prof1[6]*prof2[6];
                        pa += prof1[7]*prof2[7];


                        prof2 -= 11;

                        s[j].a = pa;

                        pga = s[j].ga;

                        s[j].ga = KALIGN_MAX(s[j-1].ga+prof2[9],s[j-1].a+prof2[8]);

                        pgb = s[j].gb;

                        s[j].gb = KALIGN_MAX(pgb+prof1[9] ,ca+prof1[8]);

                        pa = ca;
                }
                prof2 += 22;
                ca = s[j].a;

                pa = KALIGN_MAX3(pa,pga + prof2[-14],pgb + prof1[-14]);
                prof2 += 11;

                pa += prof1[0]*prof2[0];
                pa += prof1[1]*prof2[1];
                pa += prof1[2]*prof2[2];
                pa += prof1[3]*prof2[3];
                pa += prof1[4]*prof2[4];
                pa += prof1[5]*prof2[5];
                pa += prof1[6]*prof2[6];
                pa += prof1[7]*prof2[7];

                prof2 -= 11;

                s[j].a = pa;

                s[j].ga = -INFINITY;

                if (hm->endb != hm->len_b){
                        s[j].gb = KALIGN_MAX(s[j].gb+prof1[9] ,ca+prof1[8]);
                }else{
                        s[j].gb = KALIGN_MAX(s[j].gb,ca)+ prof1[10];
                }


                prof2 -= (hm->endb-hm->startb) * 22;

        }
        prof1 -= 22 * (hm->enda);
        return s;
}

struct states* backward_hirsch_dna_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm)
{
        struct states* s = hm->b;
        register float pa = 0;
        register float pga = 0;
        register float pgb = 0;
        register float ca = 0;
        register int i = 0;
        register int j = 0;

        prof1 += (hm->enda+1) * 22;
        prof2 += (hm->endb+1) * 22;
        s[hm->endb].a = s[0].a;
        s[hm->endb].ga = s[0].ga;
        s[hm->endb].gb = s[0].gb;


        //init of first row;
        //j = endb-startb;
        if(hm->endb == hm->len_b){

                for(j = hm->endb-1;j > hm->startb;j--){
                        prof2 -= 22;
                        s[j].a = -INFINITY;
                        s[j].ga = KALIGN_MAX(s[j+1].ga,s[j+1].a)+prof2[10];
                        s[j].gb = -INFINITY;
                }
                prof2 -= 22;
        }else{
                for(j = hm->endb-1;j > hm->startb;j--){
                        prof2 -= 22;
                        s[j].a = -INFINITY;
                        s[j].ga = KALIGN_MAX(s[j+1].ga+prof2[9],s[j+1].a+prof2[8]);
                        s[j].gb = -INFINITY;
                }
                prof2 -= 22;
        }

        s[hm->startb].a = -INFINITY;
        s[hm->startb].ga = -INFINITY;
        s[hm->startb].gb = -INFINITY;

        i = hm->enda-hm->starta;
        while(i--){
                prof1 -= 22;

                pa = s[hm->endb].a;
                pga = s[hm->endb].ga;
                pgb = s[hm->endb].gb;
                s[hm->endb].a = -INFINITY;
                s[hm->endb].ga = -INFINITY;

                if(hm->endb == hm->len_b){
                        s[hm->endb].gb = KALIGN_MAX(pgb,pa)+prof1[10];
                }else{
                        s[hm->endb].gb = KALIGN_MAX(pgb+prof1[9] ,pa+prof1[8]);
                }
                //j = endb-startb;
                prof2 += (hm->endb-hm->startb) *22;
                //while(j--){
                for(j = hm->endb-1;j > hm->startb;j--){
                        prof2 -= 22;
                        ca = s[j].a;

                        pa = KALIGN_MAX3(pa,pga + prof2[30],pgb + prof1[30]);

                        prof2 += 11;
                        pa += prof1[0]*prof2[0];
                        pa += prof1[1]*prof2[1];
                        pa += prof1[2]*prof2[2];
                        pa += prof1[3]*prof2[3];
                        pa += prof1[4]*prof2[4];
                        pa += prof1[5]*prof2[5];
                        pa += prof1[6]*prof2[6];
                        pa += prof1[7]*prof2[7];
                        prof2 -= 11;

                        s[j].a = pa;

                        pga = s[j].ga;

                        s[j].ga = KALIGN_MAX(s[j+1].ga+prof2[9], s[j+1].a+prof2[8]);

                        pgb = s[j].gb;

                        s[j].gb = KALIGN_MAX(pgb+prof1[9], ca+prof1[8]);

                        pa = ca;
                }
                prof2 -= 22;
                ca = s[j].a;

                pa = KALIGN_MAX3(pa,pga + prof2[30],pgb + prof1[30]);

                prof2 += 11;
                pa += prof1[0]*prof2[0];
                pa += prof1[1]*prof2[1];
                pa += prof1[2]*prof2[2];
                pa += prof1[3]*prof2[3];
                pa += prof1[4]*prof2[4];
                pa += prof1[5]*prof2[5];
                pa += prof1[6]*prof2[6];
                pa += prof1[7]*prof2[7];
                prof2 -= 11;

                s[j].a = pa;

                //pga = s[j].ga;
                s[j].ga = -INFINITY;//MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);

                //pgb = s[j].gb;
                if(hm->startb){
                        s[j].gb = KALIGN_MAX(s[j].gb+prof1[9], ca+prof1[8]);
                }else{
                        s[j].gb = KALIGN_MAX(s[j].gb,ca)+prof1[10];
                }
        }
        return s;
}





float* dna_make_profile(float* prof, int* seq,int len,float** subm)
//int* make_profile(int* prof, int* seq,int len)
{
        int i,j,c;
        prof = malloc(sizeof(float)*(len+2)*22);
        prof +=  (22 *(len+1));
        //fprintf(stderr,"Len:%d	%d\n",len,64*len);
        //for (i = 64;i--;){
        for (i = 0;i < 22;i++){
                prof[i] = 0;
        }
        //here I set the penalties...
        prof[5+11] = -1 * GPO;
        prof[6+11] = -1 * GPE;
        prof[7+11] = -1 * TGPE;


        i = len;
        while(i--){
                prof -= 22;
                //fprintf(stderr,"-64\n");
                //for (j = 64; j--;){
                for (j = 0;j < 22;j++){
                        prof[j] = 0;
                }
                c = seq[i];

                prof[c] += 1;

                //n = feature[i];
                //prof[n+23] = 1;

                prof += 11;
                for(j = 5;j--;){
                        prof[j] = subm[c][j];
                }
                prof[5] = -1 * GPO;
                prof[6] = -1 * GPE;
                prof[7] = -1 * TGPE;
                prof -= 11;
        }
        prof -= 22;
        for (i = 0;i < 22;i++){
                prof[i] = 0;
        }
        prof[5+11] = -1 * GPO;
        prof[6+11] = -1 * GPE;
        prof[7+11] = -1 * TGPE;

        return prof;
}

void dna_set_gap_penalties(float* prof,int len,int nsip)
{
        int i;

        prof +=  (22 *(len+1));
        prof[8] = prof[16]*nsip;//gap open or close
        prof[9] = prof[17]*nsip;//gap extention

        prof[10] = prof[18]*nsip;//gap open or close
        //prof[30] = prof[58]*nsip;//gap extention


        i = len+1;
        while(i--){
                prof -= 22;
                prof[8] = prof[16]*nsip;//gap open or close
                prof[9] = prof[17]*nsip;//gap extention

                prof[10] = prof[18]*nsip;//gap open or close
                //	prof[30] = prof[58]*nsip;//gap extention
        }
}



struct hirsch_mem* hirsch_mem_alloc(struct hirsch_mem* hm,int x)
{

        // a=((typeof(a))(((int)(((void *)malloc(c+15))+15))&-16)).
        hm = (struct hirsch_mem *) malloc(sizeof(struct hirsch_mem));
        hm->starta = 0;
        hm->startb = 0;
        hm->enda = 0;
        hm->endb = 0;
        hm->size = x;
        hm->len_a = 0;
        hm->len_b = 0;
        hm->f = malloc(sizeof(struct states)* (x+1));
        hm->b = malloc(sizeof(struct states)* (x+1));
        return hm;
}

struct hirsch_mem* hirsch_mem_realloc(struct hirsch_mem* hm,int x)
{
        hm->starta = 0;
        hm->startb = 0;
        hm->enda = 0;
        hm->endb = 0;
        hm->len_a = 0;
        hm->len_b = 0;
        hm->size = x;
        hm->f = realloc(hm->f,sizeof(struct states)* (x+1));
        hm->b = realloc(hm->b,sizeof(struct states)* (x+1));
        return hm;
}

void hirsch_mem_free(struct hirsch_mem* hm)
{
        free(hm->f);
        free(hm->b);
        free(hm);
}


struct alignment* detect_and_read_sequences(char* infile)
{

        char **input = 0;
        unsigned short int* input_type = 0;
        unsigned short int* input_numseq = 0;
        struct alignment* aln = 0;
        int num_input = 0;
        int i = 0;

        int free_read = 1;
        int numseq = 0;

        //aln->numseq = 0;
        num_input = 1;

        input = malloc(sizeof(char*) * num_input);
        input_type = malloc(sizeof(unsigned short int) * num_input);
        input_numseq = malloc(sizeof(unsigned short int) * num_input);

        for (i = 0; i < num_input;i++){
                input[i] = 0;
                input_type[i] = 0;
                input_numseq[i] = 0;
        }

        free_read = 0;


        input[0] = get_input_into_string(input[0],infile);
        if(input[0]){
                free_read++;

                input_numseq[0]  = count_sequences_fasta(input[0]);
                input_type[0] = 0;

                fprintf(stderr,"found %d sequences\n",input_numseq[0]);

                if(input_numseq[0] < 1){
                        free(input[0]);
                        input[0] = 0;
                }else{
                        numseq += input_numseq[0];
                }
        }



        if(numseq < 2){

                if(!aln->numseq){
                        fprintf(stderr,"\nWARNING: No sequences found.\n\n");
                }else{
                        fprintf(stderr,"\nWARNING: Only one sequence found.\n\n");
                }
                exit(EXIT_FAILURE);
        }




        //aln->numprofiles = (numseq << 1) - 1;
        aln = aln_alloc(numseq);
        aln = read_sequences(aln,input[0]);


        free(input);
        free(input_type);
        free(input_numseq);
        return aln;
}


int count_sequences_fasta(char* string)
{
        int nbytes;
        int i;
        int n = 0;
        int stop = 0;
        nbytes = (int)strlen(string);
        for (i =0;i < nbytes;i++){
                if (string[i] == '>'&& stop == 0){
                        stop = 1;
                        n++;
                }
                if (string[i] == '\n'){
                        stop = 0;
                }
        }
        if(!n){
                return 0;
        }
        return n;
}


char* get_input_into_string(char* string,char* infile)
{
        int i = 0;
        int string_length = 2;
        char c = 0;
        FILE *file = 0;
        size_t bytes_read = 0;
        if(infile){
                if (!(file = fopen( infile, "r" ))){

                        fprintf(stderr,"Cannot open file '%s'\n", infile);
                        exit(EXIT_FAILURE);
                }
                if (fseek(file,0,SEEK_END) != 0){
                        (void)fprintf(stderr, "ERROR: fseek failed\n");
                        exit(EXIT_FAILURE);
                }
                i= (int)ftell (file);
                if (fseek(file,0,SEEK_START) != 0){
                        (void)fprintf(stderr, "ERROR: fseek failed\n");
                        exit(EXIT_FAILURE);
                }
                string = malloc ((i+1)* sizeof(char));
                bytes_read= fread(string,sizeof(char), i, file);
                if(!bytes_read){
                        fprintf(stderr,"Reading from file:%s failed \n", infile );
                        exit(EXIT_FAILURE);
                }

                string[i] = 0;
                fclose(file);
        }else{
                if (!isatty(0)){
                        string = malloc(sizeof(char*)*string_length);
                        while (!feof (stdin)){
                                c = getc(stdin);
                                if (i == string_length){
                                        string_length <<=1;
                                        string = realloc(string,sizeof(char)*string_length);
                                }
                                string[i] = c;
                                i++;
                        }
                        string[i-1] = 0;
                }else{
                        return 0;
                }
        }
        return string;
}

struct alignment* read_sequences(struct alignment* aln,char* string)
{
        int c = 0;
        int n = 0;
        int i = 0;
        int j = 0;
        int stop = 0;
        int start = 0;
        int nbytes;
        int local_numseq = 0;				// O	12				//U17
        int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
        nbytes = (int)strlen(string);

        //aln = (struct alignment*) malloc(sizeof(struct alignment));
        for (i =0;i < nbytes;i++){
                if (string[i] == '>'&& stop == 0){
                        stop = 1;
                        local_numseq++;
                }
                if (string[i] == '\n'){
                        stop = 0;
                }
        }
        /*
          if(!numseq){
          fprintf(stderr,"No sequences found!\n");
          exit(1);
          }
          numprofiles = (numseq << 1) - 1;
          aln->s = malloc(sizeof(int*) * (numseq));
          aln->seq = malloc(sizeof(char*) * (numseq));
          aln->ft = 0;
          aln->si = 0;
          aln->sl = malloc(sizeof(int) * (numprofiles));
          aln->sip = malloc(sizeof(int*)* numprofiles);
          aln->nsip = malloc(sizeof(int)* numprofiles);
          aln->sn = malloc(sizeof(char*) * numseq);
          aln->lsn = malloc(sizeof(int) * numseq);

          for (i =0;i < numprofiles;i++){
          aln->sip[i] = 0;
          aln->nsip[i] = 0;
          }*/
        start = 0;
        while(aln->sl[start]){
                start++;
        }
        j = start;

        for (i =0;i < nbytes;i++){
                if (string[i] == '>' && stop == 0){
                        stop = 1;
                        aln->sl[j] =c;
                        j++;
                        c = 0;
                }
                if (string[i] == '\n'){
                        if(stop == 1){
                                aln->lsn[j-1] = n;
                                n = 0;
                        }
                        stop = 0;
                }
                if (stop == 1 && string[i] != '\n' && string[i] != 0 && string[i] != '>' ){
                        n++;
                }
                if (stop == 0 && string[i] != '\n' && string[i] != 0 ){
                        if (isalpha((int)string[i])){
                                c++;
                        }
                }
        }
        aln->sl[j] = c;

        for (i =1+start;i < local_numseq+1+start;i++){
                if(!aln->sl[i]){
                        fprintf(stderr,"Sequence %d has a length of 0!!\n",i-1);
                        exit(1);
                }
                aln->sl[i-1] = aln->sl[i];
        }
        aln->sl[start+local_numseq] = 0;

        //for (i = numseq;i--;){
        for (i = start; i < local_numseq+start;i++){
                aln->s[i] = malloc(sizeof(int)*(aln->sl[i]+1));
                aln->seq[i] = malloc(sizeof(char)*(aln->sl[i]+1));
                aln->sn[i] = malloc(sizeof(char)*(aln->lsn[i]+1));
                //aln->sip[i] = malloc(sizeof(int)*1);
                //aln->nsip[i] = 1;
                //aln->sip[i][0] = i;
        }

        stop = 0;
        j = start;
        for (i =0;i < nbytes;i++){
                if (string[i] == '>' && stop == 0 ){
                        stop = 1;
                        j++;
                        c = 0;
                }
                if (string[i] == '\n'){
                        if(stop == 1){
                                n = 0;
                        }
                        stop = 0;
                }
                if (stop == 1 &&string[i] != '\n' && string[i] != 0 && string[i] != '>' ){
                        aln->sn[j-1][n] = string[i];
                        n++;
                }
                if (stop == 0 && string[i] != '\n' && string[i] != 0 ){
                        if(isalpha((int)string[i])){
                                aln->s[j-1][c] = aacode[toupper(string[i])-65];
                                aln->seq[j-1][c] = string[i];
                                c++;
                        }
                }
        }

        for (i = start;i< local_numseq+start;i++){
                aln->s[i][aln->sl[i]] = 0;
                aln->seq[i][aln->sl[i]] = 0;
                aln->sn[i][aln->lsn[i]] = 0;
        }

        free(string);
        return aln;
}



void fasta_output(struct alignment* aln,char* outfile)
{
        int i,j,c,f;
        int tmp;
        FILE *fout = NULL;
        if(outfile){
                if ((fout = fopen(outfile, "w")) == NULL){
                        fprintf(stderr,"can't open output\n");
                        exit(0);
                }
        }else{
                fout = stderr;
        }
        for (i = 0; i < aln->numseq;i++){
                f = aln->nsip[i];
                //fprintf(fout,">%s\n",aln->sn[f]);
                c = 0;
                for (j = 0; j < aln->sl[f];j++){
                        tmp = aln->s[f][j];
                        while (tmp){
                                fprintf(fout,"-");
                                c++;
                                if(c == 160){
                                        fprintf(fout,"\n");
                                        c = 0;
                                }
                                tmp--;
                        }
                        fprintf(fout,"%c",aln->seq[f][j]);
                        c++;
                        if(c == 160){
                                fprintf(fout,"\n");
                                c = 0;
                        }
                }
                tmp = aln->s[f][aln->sl[f]];
                while (tmp){
                        fprintf(fout,"-");
                        c++;
                        if(c == 160){
                                fprintf(fout,"\n");
                                c = 0;
                        }
                        tmp--;
                }
                fprintf(fout,"\n");
        }
        if(outfile){
                fclose(fout);
        }
}



struct alignment* aln_alloc(int numseq)
{
        int i;

        struct alignment* aln = NULL;
        int numprofiles =(numseq << 1) - 1;
        aln = (struct alignment*) malloc(sizeof(struct alignment));

        aln->s = malloc(sizeof(int*) * (numseq ));
        aln->seq = malloc(sizeof(char*) * (numseq ));
        aln->si  =  malloc(sizeof(struct sequence_information* ) * (numseq));
        aln->sl = malloc(sizeof(unsigned int) * (numprofiles));
        aln->sip = malloc(sizeof(unsigned int*)* numprofiles);

        aln->nsip = malloc(sizeof(unsigned int)* numprofiles);
        aln->sn = malloc(sizeof(char*) * numseq);
        aln->lsn = malloc(sizeof(unsigned int) * numseq);
        for (i =0;i < numprofiles;i++){
                aln->sip[i] = 0;
                aln->nsip[i] = 0;
                aln->sl[i] = 0;
        }

        for(i =0;i < numseq;i++){
                aln->lsn[i] = 0;
                aln->si[i] = 0;
                aln->sip[i] = malloc(sizeof(int)*1);
                aln->nsip[i] = 1;
                aln->sip[i][0] = i;
        }
        aln->numseq = numseq;
        aln->numprofiles =  (aln->numseq << 1) - 1;
        return aln;
}



struct alignment* make_dna(struct alignment* aln)
{

        //int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
        int i,j;
        int* p;

        for(i = 0;i < aln->numseq;i++){
                p = aln->s[i];
                for (j = 0; j < aln->sl[i];j++){
                        switch(p[j]){
                        case 2: //C
                                p[j] = 1;
                                break;
                        case 6: //G
                                p[j] = 2;
                                break;
                        case 17: //T  or U
                                p[j] = 3;
                                break;
                        case 12: // N
                                p[j] = 4;
                                break;
                        case 20: // X
                                p[j] = 4;
                                break;
                        case 23://O whatever that is...
                                p[j] = 4;
                                break;
                        }
                        //	printf("%d\n",p[j]);
                }
        }
        return aln;
}


void free_aln(struct alignment* aln)
{
        int i;
        for (i = aln->numseq;i--;){
                free(aln->s[i]);
                free(aln->seq[i]);
                free(aln->sn[i]);
        }


        if(aln->si){
                free(aln->si);
        }

        for (i = aln->numprofiles;i--;){
                if(aln->sip[i]){
                        free(aln->sip[i]);
                }
        }
        free(aln->seq);
        free(aln->s);
        free(aln->sn);
        free(aln->sl);
        free(aln->lsn);
        free(aln->sip);
        free(aln->nsip);
        free(aln);
}


struct bignode* big_insert_hash(struct bignode *n,const unsigned int pos)
{
        struct bignode* p = 0;
        if(n){
                if(n->num < NODESIZE){
                        n->pos[n->num] = pos;
                        n->num++;
                        return n;
                }else{
                        p = (struct bignode*) malloc(sizeof(struct bignode));
                        p->pos[0] = pos;
                        p->num = 1;
                        p->next = n;
                }
        }else{
                p = (struct bignode*) malloc(sizeof(struct bignode));
                p->pos[0] = pos;
                p->num = 1;
                p->next = n;
        }
        return p;
}

void big_remove_nodes(struct bignode *n)
{
        struct bignode* p;
        while(n){
                p = n;
                n = n->next;
                free(p);
        }
}

void big_print_nodes(struct bignode *n)
{
        int i;
        while(n){
                for (i = 0; i < n->num;i++){
                        fprintf(stderr,"%d ",n->pos[i]);
                }
                n = n->next;
        }
}



struct aln_tree_node* real_upgma(float **dm,int numseq)
{
        int i,j;
        int *as = 0;
        float max;
        int node_a = 0;
        int node_b = 0;
        int cnode = numseq;
        int ntree = 2;
        struct aln_tree_node** tree = 0;
        struct aln_tree_node* tmp = 0;

        int numprofiles =  (numseq << 1) - 1;

        as = malloc(sizeof(int)*numseq);
        for (i = numseq; i--;){
                as[i] = i+1;
        }

        tree = malloc(sizeof(struct aln_tree_node*)*numseq);
        for (i=0;i < numseq;i++){
                tree[i] = malloc(sizeof(struct aln_tree_node));
                tree[i]->done = 1;
                tree[i]->num = i;
                tree[i]->path = 0;
                tree[i]->profile = 0;
                tree[i]->seq = 0;//seq[i];
                tree[i]->len = 0;//len[i];
                /*
                  Needs to be +2 because:
                  at n = 3 is is possible to get a perfectly balanced binary tree with 4 sequences at intermediate nodes
                */
                /*tree[i]->links = malloc(sizeof(struct aln_tree_node*)*2);

                  for ( j =0;j < 2;j++){
                  tree[i]->links[j] = 0;
                  }*/

                tree[i]->internal_lables = malloc(sizeof(int)*(ntree+(ntree-1)));
                tree[i]->links = malloc(sizeof(struct aln_tree_node*)*(ntree+(ntree-1)));

                for ( j =0;j < (ntree+(ntree-1));j++){
                        tree[i]->links[j] = 0;
                        tree[i]->internal_lables[j] = 0;
                }
        }

        while (cnode != numprofiles){
                max =   -FLT_MAX;
                for (i = 0;i < numseq-1; i++){
                        if (as[i]){
                                for ( j = i + 1;j < numseq;j++){
                                        if (as[j]){
                                                if (dm[i][j] > max){
                                                        max = dm[i][j];
                                                        node_a = i;
                                                        node_b = j;
                                                }
                                        }
                                }
                        }
                }
                //fprintf(stderr,"%d	%d	-> %d\n",node_a,node_b,cnode);
                tmp = malloc(sizeof(struct aln_tree_node));
                tmp->done = 0;
                tmp->path = 0;
                tmp->profile = 0;
                tmp->num = cnode;
                tmp->seq = 0;
                tmp->len = 0;

                tmp->links = malloc(sizeof(struct aln_tree_node*)*(ntree+(ntree-1)));
                tmp->internal_lables = malloc(sizeof(int)*(ntree+(ntree-1)));
                tmp->links[0] = tree[node_a];
                tmp->links[1] = tree[node_b];
                tmp->internal_lables[0] = cnode;
                tmp->internal_lables[1] = 0;

                for ( i =2;i < (ntree+(ntree-1));i++){
                        tmp->links[i] = 0;
                        tmp->internal_lables[i] = 0;

                }


                tree[node_a] = tmp;
                tree[node_b] = 0;

                /*deactivate  sequences to be joined*/
                as[node_a] = cnode+1;
                as[node_b] = 0;
                cnode++;

                /*calculate new distances*/
                for (j = numseq;j--;){
                        if (j != node_b){
                                dm[node_a][j] = (dm[node_a][j] + dm[node_b][j])*0.5;
                        }
                }
                dm[node_a][node_a] = 0.0f;
                for (j = numseq;j--;){
                        dm[j][node_a] = dm[node_a][j];
                        dm[j][node_b] = 0.0f;
                        dm[node_b][j] = 0.0f;
                }
        }
        tmp = tree[node_a];

        for (i = numseq;i--;){
                free(dm[i]);
        }
        free(dm);


        free(tree);
        free(as);
        return tmp;
}


int* readtree(struct aln_tree_node* p,int* tree)
{
        if(p->links[0]){
                tree = readtree(p->links[0],tree);
        }
        if(p->links[1]){
                tree = readtree(p->links[1],tree);
        }

        if(p->links[0]){
                if(p->links[1]){
                        tree[tree[0]] = p->links[0]->num;
                        tree[tree[0]+1] = p->links[1]->num;
                        tree[tree[0]+2] = p->num;
                        tree[0] +=3;
                        free(p->links[0]->internal_lables);
                        free(p->links[0]->links);
                        free(p->links[0]);
                        free(p->links[1]->internal_lables);
                        free(p->links[1]->links);
                        free(p->links[1]);
                }
        }
        return tree;
}




#ifdef ITEST
int main(int argc, char * argv[])
{

        char *tmp_seq[119] = {
                "ACAGGCTAAAGGAGGGGGCAGTCCCCA",
                "AGGCTAAAGGAGGGGGCAGTCCCCACC",
                "AGGCTAAAGGAGGGGGCAGTCCCCACC",
                "AGTCCCCACCATATTTGAGTCTTTCTC",
                "AGTGGATATCACAGGCTAAAGGAGGGG",
                "AGTGGATATCACAGGCTAAAGGAGGGG",
                "AGTGGATATCACAGGCTAAAGGAGGGG",
                "AGTGGATATCACAGGCTAAAGGAGGGG",
                "AGTGGATATCACAGGCTAAAGGAGGGT",
                "CTAAAGGAGGGGGCAGTCCCCACCATA",
                "GAGGCTAAAGGAGGGGGCAGTCCCCAT",
                "GAGGGGGCAGTCCCCACCATATTTGAG",
                "GAGGGGGCAGTCCCCACCATATTTGAG",
                "GAGGGGGCAGTCCCCACCATATTTGAG",
                "GAGGGGGCAGTCCCCACCATATTTGAT",
                "GAGGGGGCAGTCCCCACCATATTTGAT",
                "GAGGGGGCAGTCCCCACCATATTTGAT",
                "GAGTCTTTCTCCAAGTTGCGCCGGACA",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTC",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTCCCCACCATATTTGAGTCTTTT",
                "GCAGTGGATATCACAGGCTAAAGGAGT",
                "GCTAAAGGAGGGGGCAGTCCCCACCAT",
                "GGAGGGGGCAGTCCCCACCATATTTGA",
                "GGAGGGGGCAGTCCCCACCATATTTGA",
                "GGATATCACAGGCTAAAGGAGGGGGCA",
                "GGCAGTCCCCACCATATTTGAGTCTTC",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGCAGTCCCCACCATATTTGAGTCTTT",
                "GGGCAGTCACCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGGCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGCAGTCCCCACCATATTTGAGTCTT",
                "GGGGCAGTCCCCACCATATTTGAGTCT",
                "GGGGCAGTCCCCACCATATTTGAGTCT",
                "GGGGCAGTCCCCACCATATTTGAGTCT",
                "GGGGGCAGTCCCCACCATATTTGAGTC",
                "GGGGGCAGTCCCCACCATATTTGAGTC",
                "GTCCCCACCATATTTGAGTCTTTCTCT",
                "TCCCCACCATATTTGAGTCTTTCTCCA",
                "TCCCCACCATATTTGAGTCTTTCTCCA",
                "TGAGTCTTTCTCCAAGTTGCGCCGGAT",
                "TGGATATCACAGGCTAAAGGAGGGGGC"};

        int i,c;
        int num_seq = 119;
        int max_len = 0;
        char** sequences = NULL;

        for(i =0; i < num_seq;i++){
                c = strlen(tmp_seq[i]);
                if(c > max_len){
                        max_len = c;
                }
        }
        max_len++;              /* for terminating character., */

        MMALLOC(sequences, sizeof(char*)* num_seq);
        for(i = 0; i< num_seq;i++){
                sequences[i] = NULL;
                MMALLOC(sequences[i], sizeof(char) * ( max_len));
                c = strlen(tmp_seq[i]);
                strncpy(sequences[i], tmp_seq[i],max_len);
                sequences[i][c] = 0;
        }

        sequences = kalign_align(sequences, num_seq);

        for(i = 0; i< num_seq;i++){
                fprintf(stdout,"%s\n",sequences[i]);
        }
        for(i = 0; i< num_seq;i++){
                MFREE(sequences[i]);
        }
        MFREE(sequences);


return EXIT_SUCCESS;
ERROR:
return EXIT_FAILURE;
}
#endif
