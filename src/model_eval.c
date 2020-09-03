#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <libgen.h>

#include "tldevel.h"
#include "tlmisc.h"
#include "tllogsum.h"
#include "tlrng.h"

#include "tlalphabet.h"

//#include "model.h"
#include "sequence_struct.h"
#include "sequence_alloc.h"
#include "sequence_io.h"
#include "sequence_prep.h"
#include "sequences_sim.h"

#include "thread_data.h"

#include "null_model_emission.h"


#include "finite_hmm.h"
#include "finite_hmm_io.h"
#include "finite_hmm_alloc.h"
#include "finite_hmm_score.h"
#include "finite_hmm_plot.h"

#include "run_score.h"

struct parameters{
        char* in_model;
        char* output;
        char* plot_output;
        int num_threads;
        float plot_thres;
        struct rng_state* rng;
};

#define OPT_PLOT_THRES 1

static int double_cmp(const void *a, const void *b);
static int print_help(char **argv);
static int free_parameters(struct parameters* param);

int main (int argc, char *argv[])
{
        FILE* fptr = NULL;
        double** out_table = NULL;
        double* score_arr = NULL;
        struct parameters* param = NULL;
        struct fhmm* fhmm = NULL;
        struct seq_buffer* sb = NULL;


        struct seqer_thread_data** td = NULL;

        double p;
        int i,j,c;
        int num_test_seq = 10000;

        int test_lengths[3] = {100,400,1600};
        //print_program_header(argv, "Scores sequences.");

        MMALLOC(param, sizeof(struct parameters));
        param->in_model = NULL;
        param->output = NULL;
        param->plot_output = NULL;
        param->num_threads = 8;
        param->rng = NULL;
        param->plot_thres = 0.01;

        while (1){
                static struct option long_options[] ={
                        {"model",required_argument,0,'m'},
                        {"out",required_argument,0,'o'},
                        {"plot",required_argument,0,'p'},
                        {"ethres",required_argument,0,OPT_PLOT_THRES},
                        {"nthreads",required_argument,0,'t'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"o:t:m:p:",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
                case OPT_PLOT_THRES:
                        param->plot_thres = atof(optarg);
                        break;
                case 'p':
                        param->plot_output = optarg;
                        break;

                case 'o':
                        param->output = optarg;
                        break;
                case 't':
                        param->num_threads = atoi(optarg);
                        break;
                case 'm':
                        param->in_model = optarg;
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

        if(1){
                //rk_seed(42, &param->rndstate);
                RUNP(param->rng = init_rng(0));
        }else{
                //rk_randomseed(&param->rndstate);
        }
        LOG_MSG("Starting run");


        if(!param->in_model){
                RUN(print_help(argv));
                ERROR_MSG("No model file! use -m  <blah.h5>");
        }else{
                if(!my_file_exists(param->in_model)){
                        RUN(print_help(argv));
                        ERROR_MSG("The file <%s> does not exist.",param->in_model);
                }
        }

        if(!param->output){
                RUN(print_help(argv));
                ERROR_MSG("No output file! use -o   <blah.csv>");
        }else{
                if(my_file_exists(param->output)){
                       WARNING_MSG("The file %s will be over-written.",param->output);
                }
        }


        LOG_MSG("Loading model.");


        init_logsum();

        /* struct tl_seq_buffer* hits = NULL; */


        /* RUN(convert_tl_seq_buf_into_ihmm_seq_buf(hits, &sb)); */
        /* free_tl_seq_buffer(hits); */

        /* LOG_MSG("Found %d putative hits", sb->num_seq); */

        LOG_MSG("Read search fhmm");
        //sb->sequences[0]->score_arr
        RUN(read_searchfhmm(param->in_model, &fhmm));

        if(param->plot_output){
                plot_finite_hmm_dot(fhmm, param->plot_output, param->plot_thres);
        }


        RUN(galloc(&score_arr, num_test_seq));
        RUN(galloc(&out_table, num_test_seq,9));

        for(i = 0; i < 3;i++){
                LOG_MSG("Testing sequences of length: %d", test_lengths[i]);
                RUN(sim_sequences(num_test_seq, fhmm->L,test_lengths[i] ,&sb, param->rng));
                if(td){
                        RUN(resize_seqer_thread_data(td  ,(sb->max_len+2)  , fhmm->K+1));
                }else{
                        RUN(create_seqer_thread_data(&td,param->num_threads,(sb->max_len+2)  , fhmm->K+1, NULL));
                }
                RUN(run_score_sequences(fhmm,sb, td));
                for(j = 0; j < num_test_seq;j++){
                        score_arr[j] = sb->sequences[j]->score;
                }
                qsort(score_arr, num_test_seq, sizeof(double),double_cmp);
                for(j = 0; j < num_test_seq;j++){
                        p = esl_exp_surv(score_arr[j] , fhmm->tau, fhmm->lambda);
                        //fprintf(stdout,"%f %f %f\n", score_arr[j],p,p* (double) num_test_seq);
                        out_table[j][i*3+0] = score_arr[j];
                        out_table[j][i*3+1] = p;
                        out_table[j][i*3+2] = p* (double) num_test_seq;
                }
        }

        RUNP(fptr = fopen(param->output,"w"));

        for(i = 0; i < 3;i++){
                if(i){
                        fprintf(fptr,",");
                }
                fprintf(fptr,"Score_L%d,",test_lengths[i]);
                fprintf(fptr,"Pvalue_L%d,",test_lengths[i]);
                fprintf(fptr,"Evalue_L%d",test_lengths[i]);
        }
        fprintf(fptr,"\n");

        for(j = 0; j < num_test_seq;j++){
                for(i = 0; i < 9;i++){
                        if(i){
                                fprintf(fptr,",");
                        }
                        fprintf(fptr,"%f",out_table[j][i]);
                }
                fprintf(fptr,"\n");
        }
        fclose(fptr);

        /*
mat = read.table("testmodelteststats",sep = ",",header = T)
f = mat[,c(1,2)]
f$group = "L100"
colnames(f) = c("Score","Pvalue","Group")
x = f
f = mat[,c(4,5)]
f$group = "L400"
colnames(f) = c("Score","Pvalue","Group")
x = rbind(x,f)
f = mat[,c(7,8)]
f$group = "L1600"
colnames(f) = c("Score","Pvalue","Group")
x = rbind(x,f)
ggplot(x,aes(x=Score, y=Pvalue, group=Group)) + geom_line(aes(color=Group))+ geom_point(aes(color=Group)) + coord_trans(y = "log10") + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +  theme_bw()

max= min(dim(mat)[1], 1000)
f = as.data.frame(mat[1:max,3])
colnames(f) = "Evalue"
f$Group = "L100"
f$Rank = 1:max;
x = f

f = as.data.frame(mat[1:max,6])
colnames(f) = "Evalue"
f$Group = "L400"
f$Rank = 1:max;
x = rbind(x,f)

f = as.data.frame(mat[1:max,9]);
colnames(f) = "Evalue"
f$Group = "L1600"
f$Rank = 1:max;
x = rbind(x,f)
ggplot(x,aes(x=Rank, y=Evalue, group=Group)) + geom_point(aes(color=Group)) + coord_trans(y = "log10",x="log10") + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +  theme_bw()

 */
        free_ihmm_sequences(sb);
        free_seqer_thread_data(td);
        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        if(fptr){
                fclose(fptr);
        }

        free_seqer_thread_data(td);
        //
        //free_ihmm_sequences(sb);
        free_fhmm(fhmm);
        free_parameters(param);
        return EXIT_FAILURE;
}


int free_parameters(struct parameters* param)
{
        ASSERT(param != NULL, " No param found - free'd already???");
        if(param->rng){
                free_rng(param->rng);
        }
        MFREE(param);
        return OK;
ERROR:
        return FAIL;
}

int print_help(char **argv)
{
        const char usage[] = " -m <model.h5> -o <outfile> ";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nthreads","Number of threads." ,"[8]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-p <>","Create <>dot file of hmm ." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-ethres","Threshold to include edges in plotting." ,"[0.01]"  );
        return OK;
}

int double_cmp(const void *a, const void *b)
{
        const double *ia = (const double *)a;
        const double *ib = (const double *)b;
        if(*ia < *ib ){
                return 1;
        }else{
                return -1;
        }
}
