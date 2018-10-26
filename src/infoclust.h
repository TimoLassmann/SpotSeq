#ifndef INFOCLUST_H
#define INFOCLUST_H

struct paraclu_cluster{
        int* state_sequence;
        float min_density;
        float max_density;
        float total;
        int start_in_sa;
        int end_in_sa;
        int seq_id;
        int start;
        int stop;
        int len;
        int count;
};

struct motif_list{
        struct paraclu_cluster** plist;
        int num_items;
        int alloc_items;
};


#endif
