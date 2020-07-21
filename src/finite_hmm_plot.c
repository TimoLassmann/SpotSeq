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

                float tSN;
        float tNN;
        float tNB;

        float tBX;
        float tXE;

        float tEC;
        float tCC;
        float tCT;

        float tEJ;
        float tJJ;
        float tJB;
        float p,q;
        int len = 10;           /* placeholder... */
        if(1){
                q = 0.5f;
                p = (float) len / ((float)len + 3.0f);
        }else{
                q = 0.0f;
                p = (float) len / ((float)len + 2.0f);

        }

        tSN = prob2scaledprob(1.0f);
        tNN = prob2scaledprob(p);
        tNB = prob2scaledprob(1.0f - p);
        tBX = prob2scaledprob(2.0f / (float) (fhmm->K * ( fhmm->K + 1.0f)));
        tXE = prob2scaledprob(1.0f);
        //LOG_MSG("%f", scaledprob2prob(fhmm->tBX));
        tEC = prob2scaledprob(1.0f - q);
        tCC = prob2scaledprob(p);
        tCT = prob2scaledprob(1.0f - p);

        tEJ = prob2scaledprob(q);
        tJJ = prob2scaledprob(p);
        tJB = prob2scaledprob(1.0f - p);


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
        fprintf(f_ptr,"S -> N [label =\"%0.2f\"];\n", scaledprob2prob(tSN));
        fprintf(f_ptr,"N -> N [label =\"%0.2f\"];\n", scaledprob2prob(tNN));
        fprintf(f_ptr,"N -> B [label =\"%0.2f\"];\n", scaledprob2prob(tNB));
        fprintf(f_ptr,"E -> C [label =\"%0.2f\"];\n",scaledprob2prob(tEC));
        fprintf(f_ptr,"E -> J [label =\"%0.2f\"];\n", scaledprob2prob(tEJ));

        fprintf(f_ptr,"C -> C [label =\"%0.2f\"];\n", scaledprob2prob(tCC));
        fprintf(f_ptr,"C -> T [label =\"%0.2f\"];\n", scaledprob2prob(tCT));
        fprintf(f_ptr,"J -> J [label =\"%0.2f\"];\n", scaledprob2prob(tJJ));
        fprintf(f_ptr,"J -> B [label =\"%0.2f\"];\n", scaledprob2prob(tJB));

        for(i = 0; i < fhmm->K;i++){
                fprintf(f_ptr,"B -> m%d [label =\"%0.2f\"];\n",i, scaledprob2prob(tBX));
                fprintf(f_ptr,"m%d -> E [label =\"%0.2f\"];\n",i, scaledprob2prob(tXE));
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
