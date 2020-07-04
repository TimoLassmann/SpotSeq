#ifndef PST_STRUCTS_H
#define PST_STRUCTS_H


/* flat pst structure  */
/* idea : 4*N ints and 4*N floats  */

struct fpst{
        float** prob;
        int** links;
        int l;                  /* length */
        int m;                  /* what is malloced */
};

struct pst {
        struct fpst* fpst_root;
        double* background;
        float* lbg;
        float p_min;
        float gamma_min;
        float mean;
        double var;
        double a,b;
        int L;
        int len;
};


#endif
