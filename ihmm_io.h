#ifndef IHMM_IO_HEADER

#define IHMM_IO_HEADER

int write_ihmm_parameters_to_file(struct iHMM_model* model,char* filename);
struct iHMM_model* read_ihmm_parameters_from_file(char* filename);

#endif
