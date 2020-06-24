
#include "tldevel.h"

#define NULL_MODEL_EMISSION_IMPORT
#include "null_model_emission.h"

int get_null_model_emissions(double** b, int L)
{
        double sum;
        int i;
        double* background = NULL;

        RUN(galloc(&background,L));


        if(L == 20){
                /* taken from hmmer  */
                /* Function:  p7_AminoFrequencies()
                 * Incept:    SRE, Fri Jan 12 13:46:41 2007 [Janelia]
                 *
                 * Purpose:   Fills a vector <f> with amino acid background frequencies,
                 *            in [A..Y] alphabetic order, same order that Easel digital
                 *            alphabet uses. Caller must provide <f> allocated for at
                 *            least 20 floats.
                 *
                 *            These were updated 4 Sept 2007, from Swiss-Prot 50.8,
                 *            (Oct 2006), counting over 85956127 (86.0M) residues.
                 *
                 * Returns:   <eslOK> on success.
                 */
                background[0] = 0.0787945;		/* A */
                background[1] = 0.0151600;		/* C */
                background[2] = 0.0535222;		/* D */
                background[3] = 0.0668298;		/* E */
                background[4] = 0.0397062;		/* F */
                background[5] = 0.0695071;		/* G */
                background[6] = 0.0229198;		/* H */
                background[7] = 0.0590092;		/* I */
                background[8] = 0.0594422;		/* K */
                background[9] = 0.0963728;		/* L */
                background[10]= 0.0237718;		/* M */
                background[11]= 0.0414386;		/* N */
                background[12]= 0.0482904;		/* P */
                background[13]= 0.0395639;		/* Q */
                background[14]= 0.0540978;		/* R */
                background[15]= 0.0683364;		/* S */
                background[16]= 0.0540687;		/* T */
                background[17]= 0.0673417;		/* V */
                background[18]= 0.0114135;		/* W */
                background[19]= 0.0304133;		/* Y */

        }else{
                for(i = 0; i < L;i++){
                        background[i] = 1.0 / (double) L;
                }
        }
        sum = 0.0;
        for(i = 0; i < L;i++){
                sum += background[i];
        }

        ASSERT(sum == 1.0,"background sum is != 1.0");

        *b = background;
        return OK;
ERROR:
        return FAIL;
}
