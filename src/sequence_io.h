#ifndef SEQUENCE_IO_H
#define SEQUENCE_IO_H



#ifdef SEQUENCE_IO_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int read_sequences_file(struct seq_buffer** seq_buf,char* filename );

#undef SEQUENCE_IO_IMPORT
#undef EXTERN

#endif
