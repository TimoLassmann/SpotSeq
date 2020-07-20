
#include "ihmm_seq.h"
#include "finite_hmm.h"
#include "finite_hmm_io.h"
#include "finite_hmm_alloc.h"
#include "run_score.h"
//#include "model.h"
#include "model_struct.h"
#include "model_alloc.h"
#include "model_io.h"
#include <getopt.h>
#include <libgen.h>

#include "tllogsum.h"

#define BUFFER_LEN 128

struct parameters{
        char* input;
        char* out;
        int num_threads;
};

int run_create_motif_for_each_sequence(struct parameters* param);

int free_parameters(struct parameters* param);
int print_help(char **argv);

int make_logo(int** matrix,int len, int L, char* outname);

int main (int argc, char *argv[])
{
        struct parameters* param = NULL;
        int c;

        //print_program_header(argv, "Extracts motifs based on sequence labelling.");

        MMALLOC(param, sizeof(struct parameters));
        param->out = NULL;
        param->input = NULL;
        param->num_threads = 4;
        while (1){
                static struct option long_options[] ={
                        {"model",required_argument,0,'m'},
                        {"out",required_argument,0,'o'},
                        {"nthreads",required_argument,0,'t'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"m:o:t:",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
                case 't':
                        param->num_threads = atoi(optarg);
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


        RUN(run_create_motif_for_each_sequence(param));
        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        free_parameters(param);
        return EXIT_FAILURE;
}


int run_create_motif_for_each_sequence(struct parameters* param)
{
        char buffer[BUFFER_LEN];
        struct seq_buffer* sb = NULL;
        struct ihmm_model* model = NULL;
        struct ihmm_sequence* s = NULL;
        struct fhmm* fhmm = NULL;

        int i,j,c,a;
        int** count_mat = NULL;




        //int best = 0;
        ASSERT(param!= NULL, "No parameters found.");


        init_logsum();

        ASSERT(param != NULL, "No parameters.");

        /* read in sequences */
        LOG_MSG("Load sequences");
        RUNP(sb = get_sequences_from_hdf5_model(param->input, IHMM_SEQ_READ_ONLY_SEQ));
        ASSERT(sb != NULL, "No sequence Buffer");


        //RUNP(fhmm_log = alloc_fhmm());
        /* get HMM parameters  */
        LOG_MSG("Read best model");
        RUNP(fhmm =  read_best_fmodel(param->input, &c));
        LOG_MSG("Labelling sequences");
        RUN(run_label_sequences(fhmm,sb, param->num_threads ));
        LOG_MSG("Done");

        RUN(convert_fhmm_scaled_to_prob(fhmm));

        RUNP(model = read_best_imodel(param->input, &c));

        for(i = 0; i < sb->num_seq;i++){
                s = sb->sequences[i];

                RUN(galloc(&count_mat,s->seq_len,fhmm->L));
                for(j = 0;j < s->seq_len;j++){
                        for(a = 0; a < fhmm->L;a++){
                                count_mat[j][a] = model->emission_counts[a][s->label[j]]  ;
                        }
                }
                LOG_MSG("Creating motif for sequence %d", i+1);
                snprintf(buffer, BUFFER_LEN, "seqmotif_%s.png", s->name);
                RUN(make_logo(count_mat, s->seq_len, model->L,buffer));
                gfree(count_mat);
                count_mat = NULL;
        }
        free_fhmm(fhmm);
        free_ihmm_sequences(sb);
        free_ihmm_model(model);
        return OK;

ERROR:
        free_fhmm(fhmm);
        free_ihmm_sequences(sb);
        free_ihmm_model(model);
        return FAIL;

}


int free_parameters(struct parameters* param)
{
        ASSERT(param != NULL, " No param found - free'd already???");

        MFREE(param);
        return OK;
ERROR:
        return FAIL;

}

int print_help(char **argv)
{
        const char usage[] = " -m <h5 model> -out <h5 out>";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nthreads","Number of threads." ,"[8]"  );

        return OK;
}
