
#include <omp.h>
#include <getopt.h>


#include "tldevel.h"
#include "tlmisc.h"
#include "tlrng.h"
#include "tlseqbuffer.h"

#include "randomkit.h"

#include "sequence_struct.h"
#include "sequence_alloc.h"
#include "sequence_io.h"

#include "pst_build.h"

#include "model_struct.h"
#include "model_io.h"
#include "model_alloc.h"

#include "bias_model.h"

#include "thread_data.h"

#include "hmm_conversion.h"

#include "finite_hmm_stats.h"
#include "finite_hmm_alloc.h"
#include "finite_hmm_io.h"
#include "finite_hmm_score.h"


#include "run_score.h"

struct parameters{
        char* in_model;
        char* out_model;
        char* seq_db;
        char* cmd_line;
        unsigned long seed;
        rk_state rndstate;
        struct rng_state* rng;
        int num_threads;
};

#define OPT_SEQDB 1
#define OPT_SEED 2

static int run_bsm(struct parameters* param);

static int calibrate_all(struct model_bag* mb,struct seqer_thread_data** td);

static void* do_calibrate_per_model(void* threadarg);
//static int find_best_model(struct model_bag*mb, struct seq_buffer* sb, int* best);
static int find_best_model(struct model_bag*mb, struct tl_seq_buffer* sb, int* best);
static int print_help(char **argv);
static void free_param(struct parameters* param);


int main(int argc, char *argv[])
{
        struct parameters* param = NULL;
        int c;

        //print_program_header(argv, "Build HDPHMM model(s).");

        MMALLOC(param, sizeof(struct parameters));

        param->in_model = NULL;
        param->out_model = NULL;
        param->seq_db = NULL;
        param->cmd_line = NULL;
        param->seed = 0;
        param->num_threads = 8;

        param->rng = NULL;
        while (1){
                static struct option long_options[] ={
                        {"in",required_argument,0,'i'},
                        {"out",required_argument,0,'o'},
                        {"seqdb",required_argument,0,OPT_SEQDB},
                        {"seed",required_argument,0,OPT_SEED},
                        {"nthreads",required_argument,0,'t'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"i:o:t:h",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
                case OPT_SEQDB:
                        param->seq_db = optarg;
                        break;
                case OPT_SEED:
                        param->seed = atoi(optarg);
                        break;
                case 'i':
                        param->in_model = optarg;
                        break;
                case 'o':
                        param->out_model = optarg;
                        break;
                case 't':
                        param->num_threads = atoi(optarg);
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


        if(!param->in_model){
                RUN(print_help(argv));
                ERROR_MSG("No input file! use --in <model.h5>");
        }else{
                if(!my_file_exists(param->in_model)){
                        ERROR_MSG("File %s does not exist.", param->in_model);
                }
        }
        if(!param->out_model){
                RUN(print_help(argv));
                ERROR_MSG("No output file! use --out <searchmodel.h5>");
        }else{
                if(my_file_exists(param->out_model)){
                        ERROR_MSG("File %s already exists.", param->out_model);
                }
        }

        if(!param->seq_db){
                RUN(print_help(argv));
                ERROR_MSG("No seqDB use --seqdb <blah.fa>");

        }else{
                if(!my_file_exists(param->seq_db)){
                        RUN(print_help(argv));
                        ERROR_MSG("The file <%s> does not exist.",param->seq_db);
                }
        }

        if(param->seed){
                RUNP(param->rng = init_rng(param->seed));
                rk_seed(param->seed, &param->rndstate);
        }else{
                RUNP(param->rng = init_rng(0));
                rk_randomseed(&param->rndstate);
        }

        RUN(make_cmd_line(&param->cmd_line,argc,argv));
        RUN(run_bsm(param));
        free_param(param);
        return EXIT_SUCCESS;

ERROR:
        free_param(param);
        return EXIT_FAILURE;
}

int run_bsm(struct parameters* param)
{
        struct model_bag* model_bag = NULL;

        struct tl_seq_buffer* sb = NULL;
        struct seqer_thread_data** td = NULL;

        struct fhmm* bias = NULL;
        double* s = NULL;
        int i;
        int best;
        /* read sequences from in model */

        RUNP(sb = get_sequences_from_hdf5_model(param->in_model, IHMM_SEQ_READ_ONLY_SEQ));

        //RUN(convert_ihmm_seq_buf_into_tl_seq_buf(s, &sb));

        /*LOG_MSG("%d",sb->L);
        for(i = 0; i < sb->num_seq;i++){
                LOG_MSG("%s",sb->sequences[i]->name);

                }*/
        /* train PST */
        RUN(create_pst_model(param->rng,sb, NULL, param->seq_db, param->out_model,0.00001, 0.01, 20.0));

        //sb = NULL;
        /* read all models */
        RUNP(model_bag = read_model_bag_hdf5(param->in_model ));
        RUN(create_seqer_thread_data(&td,param->num_threads, 1024 , 128, &param->rndstate));
        /* convert to fhmmm  */
        RUN(convert_ihmm_to_fhmm_models(model_bag));

        /* calibrate */
        RUN(calibrate_all(model_bag, td));

        /* WARNING NEED TO ADD STORAGE FOR SCORES !!!! */
        //RUN(add_multi_model_label_and_u(s, model_bag->num_models));
        for(i = 0; i < sb->num_seq;i++){
                s = NULL;
                MMALLOC(s, sizeof(double) * model_bag->num_models);
                sb->sequences[i]->data = s;
        }
        /* score all training sequences */
        RUN(run_score_sequences( model_bag->finite_models,sb, td, model_bag->num_models, FHMM_SCORE_P_LODD));
        /* assign best */
        RUN(find_best_model(model_bag, sb, &best));
        LOG_MSG("Best model: %d",best);
        //write
        RUN(build_bias_model(model_bag->finite_models[best], &bias));
        RUN(write_biashmm(param->out_model, bias));
        RUN(write_searchfhmm(param->out_model, model_bag->finite_models[best]));

        for(i = 0; i < sb->num_seq;i++){
                MFREE(sb->sequences[i]->data);
                sb->sequences[i]->data = NULL;
        }
        free_tl_seq_buffer(sb);
        free_model_bag(model_bag);
        free_seqer_thread_data(td);
        free_fhmm(bias);
        return OK;
ERROR:
        return FAIL;
}

int find_best_model(struct model_bag*mb, struct tl_seq_buffer* sb, int* best)
{
        double* total_e = NULL;
        double* s;
        int i,j;
        double min;

        RUN(galloc(&total_e, mb->num_models));

        for(j = 0; j < mb->num_models;j++){
                total_e[j] = 0.0;
        }
        for(i= 0 ;i < sb->num_seq;i++){
                s = sb->sequences[i]->data;
                for(j = 0; j < mb->num_models;j++){
                        total_e[j] += s[j];
                        //fprintf(stdout,"%f %f ", s->sequences[i]->score_arr[j], esl_exp_logsurv(s->sequences[i]->score_arr[j], mb->finite_models[j]->tau,mb->finite_models[j]->lambda));
                }
                //fprintf(stdout,"\n");
        }
        j = -1;
        min = 1.0;
        for(i = 0; i < mb->num_models;i++){
                LOG_MSG(" Model %d: %d states: %f", i, mb->finite_models[i]->K, total_e[i]);
                if(total_e[i] < min){
                        min = total_e[i];
                        j =i;
                }
        }
        gfree(total_e);
        *best = j;
        return OK;
ERROR:
        return FAIL;
}


/* calibrate all models  */
int calibrate_all(struct model_bag* mb,struct seqer_thread_data** td)
{
        int i,j,c;
        int num_threads = td[0]->num_threads;
        int run;


        ASSERT(mb != NULL, "No models");
        c = 0;

        for(run = 0; run < mb->num_models;run+= num_threads){
                j = 0;
                for(i = 0; i < num_threads;i++){
                        td[i]->thread_ID = i;
                        td[i]->model_ID = c;
                        td[i]->fhmm = mb->finite_models;
                        //td[i]->sb = sb;
                        LOG_MSG("Cal %d",c);
                        j++;
                        c++;
                        if(c == mb->num_models){
                                break;
                        }
                }
#ifdef HAVE_OPENMP
                omp_set_num_threads( MACRO_MIN(num_threads,j));
#pragma omp parallel shared(td) private(i)
                {
#pragma omp for schedule(dynamic) nowait
#endif
                        for(i = 0; i < j;i++){
                                do_calibrate_per_model(td[i]);
                        }
#ifdef HAVE_OPENMP
                }
#endif
        }
        return OK;
ERROR:
        return FAIL;
}

void* do_calibrate_per_model(void* threadarg)
{
        struct seqer_thread_data *data;
        data = (struct seqer_thread_data *) threadarg;
        int r;
        r = rk_random(&data->rndstate);
        fhmm_calibrate(data->fhmm[data->model_ID], data->fmat, r);
        //LOG_MSG("Model %d: %f %f",   data->model_ID, data->fhmm[data->model_ID]->lambda,data->fhmm[data->model_ID]->tau);

        return NULL;
}


int print_help(char **argv)
{
        const char usage[] = " -i  <ihmm model> -out <search model>";
        char* tmp = NULL;

        RUN(tlfilename(argv[0], &tmp));
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",tmp,usage);
        fprintf(stdout,"Options:\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--seqdb","Reference database." ,"[8]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nthreads","Number of threads." ,"[8]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--seed","Seed" ,"[NA]"  );
        MFREE(tmp);
        return OK;
ERROR:
        MFREE(tmp);
        return FAIL;
}

void free_param(struct parameters* param)
{

        if(param){
                if(param->cmd_line){
                        gfree(param->cmd_line);
                }
                if(param->rng){
                        free_rng(param->rng);
                }
                MFREE(param);
        }
}
