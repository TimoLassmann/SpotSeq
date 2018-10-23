#ifndef INFOCLUST_H
#define INFOCLUST_H

struct paraclu_cluster{
        int* state_sequence;
        float max_d;
        float kl_divergence;
        int seq_id;
        int start;
        int stop;
        int len;
};

struct motif_list{
        struct paraclu_cluster** plist;
        int num_items;
        int alloc_items;
};


#endif
