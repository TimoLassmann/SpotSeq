#ifndef SEQUENCE_IO_H
#define SEQUENCE_IO_H



#ifdef SEQUENCE_IO_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

#define IHMM_SEQ_READ_ALL 0
#define IHMM_SEQ_READ_ONLY_SEQ 1

struct seq_buffer;

EXTERN int read_sequences_file(struct seq_buffer** seq_buf,char* filename );
EXTERN struct seq_buffer* get_sequences_from_hdf5_model(char* filename, int mode);
EXTERN int add_sequences_to_hdf5_model(char* filename,struct seq_buffer* sb, int num_models);


#undef SEQUENCE_IO_IMPORT
#undef EXTERN

#endif
