#include "ihmm_seq.h"

#include "finite_hmm.h"
#include "finite_hmm_io.h"
#include "run_score.h"
#include "model.h"

#include <getopt.h>
#include <libgen.h>

#include "tllogsum.h"
#include "tlhdf5wrap.h"

#define BUFFER_LEN 128

struct parameters{
        char* input;
        char* out;
        unsigned long seed;
        rk_state rndstate;
};

static int print_help(char **argv);

static int extract_motifs(struct parameters* param);

static int calc_per_state_rel_entrophy(struct fhmm* fhmm, double* rel_entropy);

static int traverse_motifs(struct fhmm* fhmm, double* rel_e, int* motif,int state, int depth);

int main(int argc, char *argv[])
{
        struct parameters* param = NULL;
        int c;

        //print_program_header(argv, "Extracts motifs based on sequence labelling.");

        MMALLOC(param, sizeof(struct parameters));
        param->out = NULL;
        param->input = NULL;
        param->seed = 0;
        while (1){
                static struct option long_options[] ={
                        {"model",required_argument,0,'m'},
                        {"out",required_argument,0,'o'},
                        {"seed",required_argument,0,'s'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"m:o:s:",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
                case 's':
                        param->seed = atoi(optarg);
                        break;

                case 'm':
                        param->input = optarg;
                        break;
                case 'o':
                        param->out = optarg;
                        break;
                case 'h':
                        RUN(print_help(argv));
                        MFREE(param);
                        exit(EXIT_SUCCESS);
                        break;
                default:
                        ERROR_MSG("not recognized");
                        break;
                }
        }

        LOG_MSG("Starting run");

        if(!param->input){
                RUN(print_help(argv));
                ERROR_MSG("No input file! use -i <blah.h5>");

        }
        if(!param->out){
                RUN(print_help(argv));
                ERROR_MSG("No output file! use -o <blah.h5>");
        }
        if(param->seed){
                rk_seed(param->seed, &param->rndstate);
        }else{
                rk_randomseed(&param->rndstate);
        }


        RUN(extract_motifs(param));
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}


int extract_motifs(struct parameters* param)
{
        struct fhmm* fhmm = NULL;
        double* rel_entropy = NULL;
        int i,j,c;

        LOG_MSG("Read best model");
        RUNP(fhmm =  read_best_fmodel(param->input, &c));
        RUN(convert_fhmm_scaled_to_prob(fhmm));
        MMALLOC(rel_entropy, sizeof(double) * fhmm->K);

        /* Step one: calculate relative entropy for each state */

        RUN(calc_per_state_rel_entrophy(fhmm, rel_entropy));

        for(i = 0; i < fhmm->K;i++){
                fprintf(stdout,"%d : %f\n",i, rel_entropy[i]);

        }

        /* remove transitions with a probability of less than 1%  */
        for(i = 0; i < fhmm->K;i++){
                c = 0;
                for(j = 0; j < fhmm->K;j++){
                        if(fhmm->t[i][j] >= 0.01){
                                fhmm->tindex[i][c+1] = j;
                                c++;
                        }
                }
                fhmm->tindex[i][0] = c+1;
        }

        int f;
        double* trans;
        for(j = 0; j < fhmm->K;j++){

                trans = fhmm->t[j];
                for(c = 1; c < fhmm->tindex[j][0];c++){
                        f = fhmm->tindex[j][c];
                        fprintf(stdout,"t:%d -> %d : %f\n",j,f,trans[f]);

                }
        }
        //exit(0);
        int motif[18];
        for(j = 2; j < fhmm->K;j++){
                traverse_motifs(fhmm, rel_entropy,  motif, j, 0);
        }



        free_fhmm(fhmm);
        return OK;
ERROR:
        free_fhmm(fhmm);
        return FAIL;
}


int traverse_motifs(struct fhmm* fhmm, double* rel_e, int* motif,int state, int depth)
{
        int i,j;
        double score;
        double prob;
        /* arrived at state motif length depth   */
        motif[depth] = state;

        if(depth > 6){
                score = 0.0;
                for(i = 0; i <= depth;i++){
                        score += rel_e[  motif[i]];
                }
                score /= (double)(depth+1);
                prob = 1.0;///= (double)(depth+1);
                for(i = 1; i <= depth;i++){
                        prob *= fhmm->t[  motif[i-1]][motif[i]];
                }

                prob = log(prob) - log( pow(1.0/ (double) fhmm->K , depth));




                if(prob < 15){

                        return OK;
                }else{
                        if(score < 1.1){
                                return OK;

                        }else{
                                for(i = 0; i <= depth;i++){
                                        fprintf(stdout,"%d ", motif[i]);
                                }
                                fprintf(stdout,"\t%f prob:%f\n",score,prob);
                        }

                }


        }
        if(depth ==11){

                return OK;
        }


        for(i = 1; i < fhmm->tindex[state][0];i++){
                j = fhmm->tindex[state][i];

                traverse_motifs(fhmm, rel_e, motif, j, depth+1);

        }

        return OK;
}



int print_help(char **argv)
{
        const char usage[] = " -m <h5 model> -out <h5 out>";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");
        /*fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--len","Minimum pattern length." ,"[10]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--frac","Minimum fraction of seq containing pattern." ,"[0.5]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-d","Resolution." ,"[100]"  );*/
        return OK;
}


int calc_per_state_rel_entrophy(struct fhmm* fhmm, double* rel_entropy)
{
        int i,j;
        double* background = NULL;
        double x;
        ASSERT(rel_entropy != NULL, "No entropy array malloced");
        ASSERT(fhmm != NULL, "No fhmm");
        background = fhmm->background;
        for(i = 0; i < fhmm->K;i++){
                x = 0.0;
                for(j = 0; j < fhmm->L;j++){
                        LOG_MSG("Emission:%d %d: %f",i,j, fhmm->e[i][j]);
                        if(fhmm->e[i][j] > 0.0){
                                x += fhmm->e[i][j] * log2f(fhmm->e[i][j] / background[j]);
                        }
                }
                rel_entropy[i] = x;
                LOG_MSG("Rel Entropy %f",rel_entropy[i]);
        }
        return OK;
ERROR:
        return FAIL;
}
