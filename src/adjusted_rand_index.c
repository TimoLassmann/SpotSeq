
#include "adjusted_rand_index.h"
/* Calculates adjusted rand index  */
int ari(long long** n,int x, int y, double* ret)
{
        long long index_sum = 0;

        long long expected_index_x = 0;
        long long expected_index_y = 0;

        double max_index = 0.0;
        double expected_index = 0.0;
        double ari = 0.0;
        long long sum_n = 0;

        long long* a = NULL;
        long long* b = NULL;


        int i,j;
        /* expected index is sumover all  (n_ij* (n_ij-1)  */

        MMALLOC(a, sizeof(long long)* x);
        MMALLOC(b, sizeof(long long)* y);
        for(j = 0; j < y;j++){
                b[j] = 0;
        }

        for(i = 0;i < x;i++){
                a[i] = 0;
                for(j = 0; j < y;j++){
                        a[i] += n[i][j];
                        b[j] += n[i][j];
                        sum_n += n[i][j];
                }
        }

        for(j = 0; j < y;j++){
                expected_index_y +=  b[j]*(b[j]-1);
        }

        for(i = 0;i<x ;i++){
                //if(a[i]){
                expected_index_x += a[i]*(a[i]-1);
                for(j = 0; j < y;j++){
                        index_sum +=  n[i][j]*(n[i][j] -1);
                }
                //}
        }
        max_index = (double)(expected_index_x + expected_index_y) / 2.0;

        expected_index = (double)(expected_index_x * expected_index_y) / (double) (sum_n *(sum_n-1));

        ari = ((double)index_sum - expected_index) /(max_index - expected_index);
        *ret = ari;

        MFREE(a);
        MFREE(b);
        return OK;
ERROR:
        if(a){
                MFREE(a);
        }
        if(b){
                MFREE(b);
        }

        return FAIL;
}


#ifdef ITEST


int main (int argc, char *argv[])
{
        long long** matrix = NULL;
         double adjusted_rand_index = 0.0;
        int i,j;
        MMALLOC(matrix , sizeof(long long*) *3);

        for(i =0 ; i < 3;i++){
                matrix[i] = NULL;
                MMALLOC(matrix[i] , sizeof(long long) *3);
        }
        matrix[0][0] = 1;
        matrix[0][1] = 1;
        matrix[0][2] = 0;

        matrix[1][0] = 1;
        matrix[1][1] = 2;
        matrix[1][2] = 1;

        matrix[2][0] = 0;
        matrix[2][1] = 0;
        matrix[2][2] = 4;






        RUN(ari(matrix,3,3,& adjusted_rand_index));

        LOG_MSG("adjusted Rand Index:%f",  adjusted_rand_index);
        for(i = 0 ; i  < 3; i++){

                MFREE(matrix[i]);
        }
        MFREE(matrix);
        matrix = NULL;

        MMALLOC(matrix , sizeof(long long*) *4);

        for(i =0 ; i < 4;i++){
                matrix[i] = NULL;
                MMALLOC(matrix[i] , sizeof(long long) *4);
        }
        matrix[0][0] = 55;
        matrix[0][1] = 1;
        matrix[0][2] = 1;
        matrix[0][3] = 1;

        matrix[1][0] = 10;
        matrix[1][1] = 76;
        matrix[1][2] = 1;
        matrix[1][3] = 1;

        matrix[2][0] = 3;
        matrix[2][1] = 2;
        matrix[2][2] = 26;
        matrix[2][3] = 1;

        matrix[3][0] = 6;
        matrix[3][1] = 2;
        matrix[3][2] = 4;
        matrix[3][3] = 45;

        RUN(ari(matrix,4,4,& adjusted_rand_index));

        LOG_MSG("adjusted Rand Index:%f",  adjusted_rand_index);

        for(i = 0 ; i  < 4; i++){

                MFREE(matrix[i]);
        }
        MFREE(matrix);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}


#endif
