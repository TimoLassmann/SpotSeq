#ifndef INFOCLUST_H
#define INFOCLUST_H


struct hit{
        int seq;
        int pos;
};

struct paraclu_cluster{
        struct hit** hits;
        double** freq_matrix;
        double** count_matrix;
        int* state_sequence;
        int* present_in_seq;
        int num_present_in_seq;

        int num_seq;
        double min_density;
        double max_density;
        double total;
        double log_likelihood;
        int start_in_sa;
        int end_in_sa;
        int seq_id;
        int start;
        int stop;
        int len;
        int count;
        int num_hits;
};

struct motif_list{
        struct paraclu_cluster** plist;
        int num_items;
        int alloc_items;


};


#endif
