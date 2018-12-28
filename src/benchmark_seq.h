#ifndef BENCHMARK_SEQ_H
#define BENCHMARK_SEQ_H



struct sequence{
        uint8_t* seq;
        char* name;
        int malloc_len;
        int seq_len;
};

struct seq_buffer{
        struct sequence** seqs;
        int malloc_num;
        int num_seq;
};


/* seqbuffer stuff  */
extern  struct seq_buffer* alloc_seq_buffer(int num_seq);
extern int write_sequences_to_file(struct seq_buffer* sb,char* filename);
extern int reset_sb(struct seq_buffer* sb);
extern void free_sb(struct seq_buffer* sb);
#endif
