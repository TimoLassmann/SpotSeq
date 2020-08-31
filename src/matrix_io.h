
#ifndef MATRIX_IO_HEADER

#define MATRIX_IO_HEADER 


#include <stdlib.h>
#include <stdio.h>



#include "tldevel.h"


struct double_matrix{
        double** matrix;
        char** col_names;
        char** row_names;
        void* matrix_mem;
        uint8_t* label;
        int nrow;
        int ncol;
        int name_len;
        int real_sample;
};


struct double_matrix* read_double_matrix(char* filename, int has_col_names,int has_row_names);
struct double_matrix* alloc_double_matrix(int ncol,int nrow, int name_len);
struct double_matrix* transpose_double_matrix(struct double_matrix* m);
int add_rows_double_matrix(struct double_matrix* m, int n_extra);
int remove_rows_double_matrix(struct double_matrix* m, int n_extra);
int fill_random_matrix(struct double_matrix* m);
int shuffle_double_matrix(struct double_matrix* m);
int print_double_matrix(struct double_matrix* m,FILE* file, int has_col_names,int has_row_names);
void free_double_matrix(struct double_matrix* m);

#endif
