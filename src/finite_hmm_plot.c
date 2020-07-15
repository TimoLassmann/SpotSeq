#include "tldevel.h"
#include "tllogsum.h"

#include "global.h"             /* need for START_STATE etc... macros */

#define FINITE_HMM_PLOT_IMPORT
#include "finite_hmm_plot.h"


/* this is intended for debugging */
int plot_finite_hmm_dot(struct fhmm* fhmm,char* filename)
{

        int i,j,c;

        FILE* f_ptr = NULL;

        RUNP(f_ptr = fopen(filename, "w"));

        fprintf(f_ptr,"digraph G\n{\nranksep = 1.0; size = \"10,10\";\n{\nnode [shape = plaintext, fontsize = 20];\nA -> B -> C ;}\n");

        fprintf(f_ptr,"{\n");

        for(i = 0; i < fhmm->K;i++){
                c = fhmm->tindex[i][0];
                for(j = 1; j < c;j++){
                        //if(i < 7 && fhmm->tindex[i][j] < 7){[ label = "SS(B)" ];
                        fprintf(f_ptr,"%d -> %d [label = \"%0.2f\"];\n",i,fhmm->tindex[i][j],scaledprob2prob( fhmm->t[i][fhmm->tindex[i][j]]));
                                //}
                }
        }

        fprintf(f_ptr,"\n}\n");



        fprintf(f_ptr,"}\n");
        fclose(f_ptr);


        return OK;
ERROR:
        return FAIL;
}
