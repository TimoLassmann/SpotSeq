#include "tldevel.h"
#include "tllogsum.h"

#include "global.h"             /* need for START_STATE etc... macros */

#define FINITE_HMM_PLOT_IMPORT
#include "finite_hmm_plot.h"


/* this is intended for debugging */
int plot_finite_hmm_dot(struct fhmm* fhmm,char* filename,float thres)
{

        int i,j,c;

        FILE* f_ptr = NULL;

        RUNP(f_ptr = fopen(filename, "w"));

        fprintf(f_ptr,"digraph finite_hmm {\nrankdir=LR;\nsize=\"8,5\"\n");

        /* define nodes  */
        fprintf(f_ptr,"node [ penwidth=3,shape = circle,color=deeppink,fontcolor=deeppink]; S B E T;\n");

        fprintf(f_ptr,"node [ penwidth=3,shape = diamond,color=deeppink,fontcolor=deeppink]; N J C;\n");
        for(i = 0; i < fhmm->K;i++){
                fprintf(f_ptr,"node [ penwidth=3,shape = box,color=black,fontcolor=black]; m%d;\n",i);
        }
        /* draw edges */
        fprintf(f_ptr,"edge [color=deeppink,fontcolor=deeppink];\n");
        fprintf(f_ptr,"S -> N [label =\"%0.2f\"];\n", scaledprob2prob(fhmm->tSN));
        fprintf(f_ptr,"N -> N [label =\"%0.2f\"];\n", scaledprob2prob(fhmm->tNN));
        fprintf(f_ptr,"N -> B [label =\"%0.2f\"];\n", scaledprob2prob(fhmm->tNB));
        fprintf(f_ptr,"E -> C [label =\"%0.2f\"];\n",scaledprob2prob(fhmm->tEC));
        fprintf(f_ptr,"E -> J [label =\"%0.2f\"];\n", scaledprob2prob(fhmm->tEJ));

        fprintf(f_ptr,"C -> C [label =\"%0.2f\"];\n", scaledprob2prob(fhmm->tCC));
        fprintf(f_ptr,"C -> T [label =\"%0.2f\"];\n", scaledprob2prob(fhmm->tCT));
        fprintf(f_ptr,"J -> J [label =\"%0.2f\"];\n", scaledprob2prob(fhmm->tJJ));
        fprintf(f_ptr,"J -> B [label =\"%0.2f\"];\n", scaledprob2prob(fhmm->tJB));

        for(i = 0; i < fhmm->K;i++){
                fprintf(f_ptr,"B -> m%d [label =\"%0.2f\"];\n",i, scaledprob2prob(fhmm->tBX));
                fprintf(f_ptr,"m%d -> E [label =\"%0.2f\"];\n",i, scaledprob2prob(fhmm->tXE));
        }

        fprintf(f_ptr,"edge [color=black,fontcolor=black];\n");


        for(i = 0; i < fhmm->K;i++){
                c = fhmm->tindex[i][0];
                for(j = 1; j < c;j++){
                        if(scaledprob2prob( fhmm->t[i][fhmm->tindex[i][j]]) >= thres){
                        fprintf(f_ptr,"m%d -> m%d [label = \"%0.2f\"];\n",i,fhmm->tindex[i][j],scaledprob2prob( fhmm->t[i][fhmm->tindex[i][j]]));
                        }
                }
        }

        fprintf(f_ptr,"\n}\n");
        fclose(f_ptr);


        return OK;
ERROR:
        return FAIL;
}
