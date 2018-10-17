#include "infoclust.h"
#include "ihmm_seq.h"
#include "finite_hmm.h"
#include <getopt.h>

struct parameters{
        char* input;
        char* out;
};

static int run_infoclust(struct parameters* param);

int calc_per_state_rel_entrophy(struct fhmm* fhmm, float* rel_entropy);

struct paraclu_cluster* init_paraclu_cluster(float min_d,int num_sequences);
void free_paraclu_cluster(struct paraclu_cluster* p);

int max_score_segment(float* x , int start ,int end, float min_density,struct paraclu_cluster* p);
float weakestPrefix(float* x , int start ,int end, int* min_prefix_pos, float* min_prefix);
void weakestSuffix(float* x , int start ,int end, int* min_suffix_pos, float* min_suffix);
int free_parameters(struct parameters* param);
int print_help(char **argv);


/* Here I plan to:
1) project information content of states onto labelled sequences
2) visualise
3) use a paraclu like algorithm to look for clusters of high information content
4) islands of high IC should somehow be compared across sequences to find shared patterns
 */

int main (int argc, char *argv[])
{
        struct parameters* param = NULL;
        int c;

        tlog.echo_build_config();

        MMALLOC(param, sizeof(struct parameters));
        param->out = NULL;
        param->input = NULL;


        while (1){
                static struct option long_options[] ={
                        {"model",required_argument,0,'m'},
                        {"out",required_argument,0,'o'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"m:o:",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
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

        RUN(run_infoclust(param));
        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        free_parameters(param);
        return EXIT_FAILURE;
}


int run_infoclust(struct parameters* param)
{
        struct seq_buffer* sb = NULL;
        struct ihmm_sequence* s = NULL;
        struct fhmm* fhmm = NULL;
        struct paraclu_cluster* p =  NULL;
        float* rel_entropy = NULL;
        int i,j,c;
        ASSERT(param != NULL, "No parameters.");

        /* read in sequences */
        RUNP(sb = get_sequences_from_hdf5_model(param->input));
        ASSERT(sb != NULL, "No sequence Buffer");

        /* read in finite hmm */
        RUNP(fhmm = alloc_fhmm());
        /* get HMM parameters  */
        RUN(read_hmm_parameters(fhmm,param->input));

        MMALLOC(rel_entropy, sizeof(float) * fhmm->K);

        /* Step one: calculate relative entropy for each state */
        RUN(calc_per_state_rel_entrophy(fhmm, rel_entropy));

        /* Don't need fhmm anymore */
        free_fhmm(fhmm);
        /* Step two: label sequences with rel entropy perstate */
        for(i = 0; i < sb->num_seq;i++){
                s = sb->sequences[i];
                for(j = 0; j < s->seq_len;j++){
                        s->u[j] = rel_entropy[s->label[j]];
                        //fprintf(stdout,"%d %d %f\n",j, s->seq[j],s->u[j]);
                }
                //fprintf(stdout,"\n");
        }

        p = init_paraclu_cluster(1.0, sb->num_seq);

         for(i = 0; i < sb->num_seq;i++){
                 s = sb->sequences[i];
                 for(j = 0; j < 5;j++){
                         RUN(max_score_segment(sb->sequences[i]->u, 0, sb->sequences[i]->seq_len,1,p));
                         if(p->start == -1){
                                 break;
                         }
                         /* print cluster  */
                         fprintf(stderr,"CLUSTER:%d:%d	%d	%d	%d	%f %f %f\n",i,j,p->start, p->stop,p->stop-p->start, p->kl_divergence,p->max_d ,p->kl_divergence / (float)(p->stop-p->start) );
                         /* delete motif */
                         for(c = p->start;c < p->stop;c++){
                                 sb->sequences[i]->u[c] = 0.0f;
                         }

                         p->max_d = 0;
                         p->start = -1;
                         p->stop = -1;
                 }
         }
        //run_paraclu_clustering(sb, "GAGA");
        /* clean up */
        free_ihmm_sequences(sb);

        MFREE(rel_entropy);
        return OK;
ERROR:
        if(fhmm){
                free_fhmm(fhmm);
        }
        if(sb){
                free_ihmm_sequences(sb);
        }
        MFREE(rel_entropy);
        return FAIL;
}

int calc_per_state_rel_entrophy(struct fhmm* fhmm, float* rel_entropy)
{
        int i,j;
        float* background = NULL;
        float x;
        ASSERT(rel_entropy != NULL, "No entropy array malloced");
        ASSERT(fhmm != NULL, "No fhmm");
        background = fhmm->background;
        for(i = 0; i < fhmm->K;i++){
                x = 0.0f;
                for(j = 0; j < fhmm->L;j++){
                        if(fhmm->e[i][j] > 0.0f){
                                x += fhmm->e[i][j] * log2f(fhmm->e[i][j] / background[j]);
                        }
                }
                rel_entropy[i] = x;
                /*fprintf(stdout,"State %d: %f\t",i, x);
                for(j = 0; j < fhmm->L;j++){
                        fprintf(stdout,"%f ", fhmm->e[i][j]);
                }
                fprintf(stdout,"\n");*/
        }
        return OK;
ERROR:
        return FAIL;
}

struct paraclu_cluster* init_paraclu_cluster(float min_d,int num_sequences)
{
        struct paraclu_cluster* p = NULL;

        MMALLOC(p, sizeof(struct paraclu_cluster));
        p->max_d = 0;
        p->seq_id = -1;
        p->start = -1;
        p->stop = -1;
        p->kl_divergence = -1;
        p->min_d_parameter = 1.0;
        p->min_len_parameter = 6;
        return p;
ERROR:
        return NULL;
}

void free_paraclu_cluster(struct paraclu_cluster* p)
{
        if(p){
                MFREE(p);
        }
}

int max_score_segment(float* x , int start ,int end, float min_density,struct paraclu_cluster* p)
{
        float max_density;
        float new_min_density;
        float min_prefix;
        float min_suffix;
        float total;
        int min_prefix_pos;
        int min_suffix_pos;
        int mid;
        //fprintf(stderr,"Looking at %d	%d\n",start,end);
        if (start == end){
                return OK;
        }
        total = weakestPrefix(x,start,end,&min_prefix_pos, &min_prefix);

        if (total < 0.0){
                return OK;
        }
        weakestSuffix(x,start,end,&min_suffix_pos,  &min_suffix);

        max_density =  (min_prefix < min_suffix) ? min_prefix : min_suffix;


        if (max_density > min_density &&  end-start >= p->min_len_parameter ){
                if(max_density /  min_density >= p->max_d){
                        p->max_d =  max_density /  min_density;
                        p->start = start;
                        p->stop = end;
                        p->kl_divergence = total;
                        //fprintf(stderr,"CLUSTER:	%d	%d	%d	%f	%e	%e	%f\n",start,end,end-start, total/(float)(end-start), min_density,max_density,max_density /  min_density);
                }
                //fprintf(stderr,"CLUSTER:	%d	%d	%d	%f	%e	%e	%f\n",start,end,end-start, total/(float)(end-start), min_density,max_density,max_density /  min_density);
        }

        if (max_density < 1e100) {
                mid = (min_prefix < min_suffix) ?  min_prefix_pos: min_suffix_pos;
                new_min_density =  (max_density > min_density) ? max_density : min_density;
                //Ci mid = (minPrefixDensity < minSuffixDensity) ? minPrefix : minSuffix;
                //fprintf(stderr,"test:%f	%f\n",new_min_density,max_density);
                RUN(max_score_segment(x ,start , mid , new_min_density, p));
                //fprintf(stderr,"test:%f	%f\n",new_min_density,max_density);
                //writeClusters(beg, mid, newMinDensity);
                RUN(max_score_segment(x ,mid , end , new_min_density, p));
                //writeClusters(mid, end, newMinDensity);
        }
        return OK;
ERROR:
        return FAIL;
}

float weakestPrefix(float* x , int start ,int end, int* min_prefix_pos, float* min_prefix)
{
        int origin = start;
        float density  = 0.0;
        //minPrefix = beg;
        *min_prefix = 1e100;
        float totalValue = x[start];
        density = totalValue / (float)(start  - origin);
        if (density < *min_prefix) {
                *min_prefix_pos = start;
                *min_prefix = density;
        }
        ++start;

        while (start < end) {
                //fprintf(stderr,"%d	%d	%f\n",origin, start,  totalValue / (float)(start  - origin));
                density = totalValue / (float)(start  - origin);
                if (density < *min_prefix) {
                        *min_prefix_pos = start;
                        *min_prefix = density;
                }
                totalValue += x[start];
                start++;
        }
        //fprintf(stderr,"%f	return\n",totalValue);
        return totalValue;
}

void weakestSuffix(float* x , int start ,int end, int* min_suffix_pos, float* min_suffix)
{

        --end;
        int origin = end;
        float density  = 0.0;
        //minSuffix = end + 1;
        *min_suffix = 1e100;
        float totalValue = x[end];

        while (end > start) {
                --end;
                //fprintf(stderr,"%d	%d	%f\n",origin, end,  totalValue / (float)( origin-end));
                density = totalValue / (origin - end);
                if (density < *min_suffix) {
                        *min_suffix_pos = end + 1;
                        *min_suffix = density;
                }
                totalValue += x[end];//end->value;
        }
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
        /*fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--len","Minimum pattern length." ,"[10]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--frac","Minimum fraction of seq containing pattern." ,"[0.5]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"-d","Resolution." ,"[100]"  );*/
        return OK;
}


float* shuffle_arr_r(float* arr,int n,unsigned int* seed)
{
        int i,j;
        float tmp;
        for (i = 0; i < n - 1; i++) {
                j = i +  (int) (rand_r(seed) % (int) (n-i));
                tmp = arr[j];
                arr[j] = arr[i];
                arr[i] = tmp;
        }
        return arr;
}
