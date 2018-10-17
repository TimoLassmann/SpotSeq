#ifndef INFOCLUST_H
#define INFOCLUST_H

struct paraclu_cluster{
        int seq_id;
        int start;
        int stop;
        float max_d;
        float kl_divergence;

        float min_d_parameter;
        int min_len_parameter;
};

#endif
