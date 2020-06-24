


#include "tldevel.h"

#include <string.h>

#define RANDOMKIT_TL_ADD_IMPORT
#include "randomkit_TL_add.h"



#ifdef TEST_COPY_RK
int main(void)
{
        int i;
        rk_state a,b;
        rk_randomseed(&a);

        RUN(copy_rk_state(&a,&b));
        for (i = 0; i < 100;i++){
                fprintf(stdout,"%f %f\n",rk_double(&a),rk_double(&b));
        }
        return OK;
ERROR:
        return FAIL;
}

#endif

int copy_rk_state(rk_state* source, rk_state* target)
{
        //int i;

        memcpy((void *)target, (void *)source, sizeof (struct rk_state_));
        return OK;
}

int compare_rk_state(rk_state* a, rk_state* b)
{
        int i;
        for(i = 0; i< RK_STATE_LEN;i++){
                if(a->key[i] != b->key[i]){
                        WARNING_MSG("Random state differs: key %d: %lud\t%lud", i, a->key[i],b->key[i]);
                }
        }
        if(a->pos != b->pos){
                WARNING_MSG("Random state differs: key %d: %d\t%d", i, a->pos, b->pos);
        }
        if(a->gauss != b->gauss){
                WARNING_MSG("Random state differs: key %d: %d\t%d", i, a->gauss,b->gauss);
        }
        if(a->psave != b->psave){
                WARNING_MSG("Random state differs: key %d: %d\t%d", i, a->psave,b->psave);
        }
        return OK;
}
/*
   unsigned long key[RK_STATE_LEN];
    int pos;
    int has_gauss;
    double gauss;

    int has_binomial;
    double psave;
    long nsave;
    double r;
    double q;
    double fm;
    long m;
    double p1;
    double xm;
    double xl;
    double xr;
    double c;
    double laml;
    double lamr;
    double p2;
    double p3;
    double p4;

*/
