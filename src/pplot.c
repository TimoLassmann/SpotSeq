
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <cairo.h>
#include "tldevel.h"

#include "ihmm_seq.h"

//#include "beam_sample.h"
#include "model.h"
#include "fast_hmm_param.h"
#include "hmm_conversion.h"


/*char* protein_colors[ALPHABET_PROTEIN] = {
        "#FF9966",
        "#009999",
        "#FF0000",
        "#CC0033",
        "#00FF00",
        "#f2f20c",
        "#660033",
        "#CC9933",
        "#663300",
        "#FF9933",
        "#CC99CC",
        "#336666",
        "#0099FF",
        "#6666CC",
        "#990000",
        "#0000FF",
        "#00FFFF",
        "#FFCC33",
        "#66CC66",
        "#006600"};

 {"#cbf751", "#5ec0cc", "#ffdf59", "#b51f16"};
*/

int nuc_colors[4] =  {
        0xcbf751,
        0x5ec0cc,
        0xffdf59,
        0xb51f16
};

int protein_colors[20] = {
        0xFF9966,
        0x009999,
        0xFF0000,
        0xCC0033,
        0x00FF00,
        0xf2f20c,
        0x660033,
        0xCC9933,
        0x663300,
        0xFF9933,
        0xCC99CC,
        0x336666,
        0x0099FF,
        0x6666CC,
        0x990000,
        0x0000FF,
        0x00FFFF,
        0xFFCC33,
        0x66CC66,
        0x006600};



#define WIDTH 48
#define HEIGHT 142


#define STRIDE (WIDTH * 4)

#define BASE_FONT_SIZE 60


static int sort_logo_letter_by_height(const void *a, const void *b);

struct logo_letter{
        double scale;
        double freq;
        char c;
        int letter;
};


struct logo_data{
        struct logo_letter*** letters;
        int len;
        int L;
};

struct parameters{
        char* input;
        char* output;
        float edge_threshold;
        int node_count_cutoff;
};

int free_parameters(struct parameters* param);

int print_help(char **argv);


int run_pplot_ihmm(struct parameters* param);
int make_pretty_plot_file(struct fhmm* fhmm, struct ihmm_model* model, struct parameters* param);
unsigned char image[STRIDE*HEIGHT];
int write_letter(cairo_t *cr, double x, double y,int c, double scale, int nuc);
int  run_draw_logo(struct logo_data* logo, char* outname);
int get_rgb_color(int color, double*r, double *g, double *b);



int main (int argc, char *argv[])
{
        struct parameters* param = NULL;
        int c;

        print_program_header(argv, "Generates iHMM image using a <.dot> file and motif png's.");

        MMALLOC(param, sizeof(struct parameters));
        param->input = NULL;
        param->output = NULL;
        param->edge_threshold = 0.5f;
        param->node_count_cutoff = 0;
        while (1){
                static struct option long_options[] ={
                        {"model",required_argument,0,'m'},
                        {"out",required_argument,0,'o'},
                        {"ethres",required_argument,0,'e'},
                        {"nthres",required_argument,0,'n'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
                int option_index = 0;
                c = getopt_long_only (argc, argv,"hm:o:e:n:",long_options, &option_index);

                if (c == -1){
                        break;
                }
                switch(c) {
                case 'n':
                        param->node_count_cutoff = atoi(optarg);
                        break;
                case 'e':
                        param->edge_threshold = atof(optarg);
                        break;
                case 'm':
                        param->input = optarg;
                        break;
                case 'o':
                        param->output = optarg;
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
                ERROR_MSG("No input file! use --in <blah.fa>");

        }else{
                if(!my_file_exists(param->input)){
                        RUN(print_help(argv));
                        ERROR_MSG("The file <%s> does not exist.",param->input);
                }
        }

        if(!param->output){
                RUN(print_help(argv));
                ERROR_MSG("No output file! use --out <blah.fa>");
        }else{
                if(my_file_exists(param->output)){
                        WARNING_MSG("Will overwrite: %s.",param->output);
                }
        }
        init_logsum();

        RUN(run_pplot_ihmm(param));

        //RUN(plot_model_entropy(param));

        //RUN(run_plot_positional_state_distribution(param));

        RUN(free_parameters(param));
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        free_parameters(param);
        return EXIT_FAILURE;
}

int run_pplot_ihmm(struct parameters* param)
{
        struct fast_hmm_param* ft = NULL;
        struct ihmm_model* model = NULL;
        struct fhmm* fhmm = NULL;

        int best = 0;
        ASSERT(param!= NULL, "No parameters found.");

        RUNP(fhmm = read_best_fmodel(param->input, &best));
        RUNP(model = read_best_imodel(param->input, &best));

        RUN(convert_fhmm_scaled_to_prob(fhmm));
        //RUNP(model = read_model_hdf5(param->input));


        //RUNP(ft = alloc_fast_hmm_param(initial_states,model->L));
        //RUN(print_fast_hmm_params(ft));
        //RUN(fill_background_emission_from_model(ft,model));
        //RUN(fill_fast_transitions_only_matrices(model,ft));
        RUN( make_pretty_plot_file( fhmm, model, param));

        free_fhmm(fhmm);
        free_ihmm_model(model);
        free_fast_hmm_param(ft);

        return OK;
ERROR:
        free_fast_hmm_param(ft);
        free_ihmm_model(model);
        return FAIL;
}

int make_pretty_plot_file(struct fhmm* fhmm, struct ihmm_model* model, struct parameters* param)
{
        int i,j,a,b;

        double sum;
        double max;
        double threshold;

        /* check if incoming edges are skewed towards one state */


        for(j =0; j < fhmm->K;j++){
                sum = 0.0;
                for(i =0; i < fhmm->K;i++){
                        sum +=  fhmm->t[i][j];
                }
                max= 0.0;

                for(i =0; i < fhmm->K;i++){
                        if(fhmm->t[i][j] / sum > max){
                                max = fhmm->t[i][j] / sum;
                                a = i;
                                b = j;
                        }
                }
                if(max > 0.7){
                        fprintf(stdout,"State %d -> %d \n",a,b);
                }
        }
        return OK;
ERROR:
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
        const char usage[] = " -m <input model> -o <output dot file>";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nthres","Distance between seeds." ,"[0]"  );

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--ethres","Edge threshold." ,"[0.5]"  );
        return OK;
}





int get_rgb_color(int color, double*r, double *g, double *b)
{


        *r = (double)((color >> 16) & 0xFF) / (double) 0xFF;
        *g = (double)((color >> 8) & 0xFF) / (double) 0xFF;
        *b = (double)(color & 0xFF) / (double) 0xFF;
        return OK;
}

struct logo_data* alloc_logo(int len, int L);
void free_logo(struct logo_data* logo);

struct logo_data* alloc_logo(int len, int L)
{
        struct logo_data* logo = NULL;
        int i,j;


        ASSERT(L > 3, "weird alphabet");
        ASSERT(len > 0 , "length 0 motif????");

        MMALLOC(logo, sizeof(struct logo_data));

        logo->L = L;
        logo->len =  len;
        logo->letters = NULL;

        MMALLOC(logo->letters, sizeof(struct logo_letter**)* logo->len);

        for(i = 0; i < logo->len;i++ ){
                logo->letters[i] = NULL;
                MMALLOC(logo->letters[i], sizeof(struct logo_letter*) * logo->L);
                for(j = 0; j < logo->L;j++){
                        logo->letters[i][j] = NULL;
                        MMALLOC(logo->letters[i][j], sizeof(struct logo_letter) );
                        logo->letters[i][j]->c = 'A';
                        logo->letters[i][j]->scale = 1.0;
                        logo->letters[i][j]->freq = 0.0;
                        logo->letters[i][j]->letter = j;

                }

        }


        return logo;
ERROR:
        return NULL;

}

void free_logo(struct logo_data* logo)
{
        int i,j;
        if(logo){
                for(i = 0; i < logo->len;i++ ){
                        for(j = 0; j < logo->L;j++){
                                MFREE(logo->letters[i][j]);
                        }
                        MFREE(logo->letters[i]);
                }
                MFREE(logo->letters);
                MFREE(logo);

        }

}

/*int main (void)
{
        struct logo_data* logo = NULL;
        int i,j;


       int m_len = 11;
        int m_L = 4;


        int motif[11][4] = {

                {	6,  2,	8 ,	1},
                { 3,	5,	9 ,	0},
                {	0, 	0,	0, 	17},
                {	5, 	2, 	17, 0},
                { 17, 0, 	0 ,	0},
                {	0,	1,	0 ,	1},
                {	3, 	2 ,	3 ,	9},
                { 4,	7 ,	2 ,	4},
                { 9,	68 ,	1 ,	1},
                { 4,	3 ,	7 ,	3},
                { 6,	33 ,	13 ,	7}
        };


        RUNP(logo = alloc_logo(m_len, m_L));



        double total;
        total = 0;
        for(i = 0; i < m_len;i++){
                for(j = 0; j < m_L;j++){
                        logo->letters[i][j]->c = "ACGT"[j];
                        logo->letters[i][j]->freq = motif[i][j];
                        total += motif[i][j];
                }
        }

        double entropy, sum, e, height;

        e = 1.0 / log(2.0) * ((double)(logo->L-1.0) / (2.0*total));



        for(i = 0; i < m_len;i++){
                sum = 0.0;
                for(j = 0; j < m_L;j++){

                        sum += logo->letters[i][j]->freq;
                }
                ASSERT(sum > 0.0, "empty column!");
                for(j = 0; j < m_L;j++){
                        logo->letters[i][j]->freq /= sum;
                        //fprintf(stdout,"%0.3f ", logo->letters[i][j]->freq);
                }
                entropy = 0.0;
                for(j = 0; j < m_L;j++){
                        if(logo->letters[i][j]->freq){
                                entropy += logo->letters[i][j]->freq *log2(logo->letters[i][j]->freq);
                        }
                }
                entropy = fabs(entropy);


                height =  log2((double) logo->L ) - (entropy + e);
                //fprintf(stdout,"entropy: %f error:%f   height :%f \n",entropy,e, log2(4.0) - (entropy + e));
                for(j = 0; j < m_L;j++){
                        logo->letters[i][j]->scale = logo->letters[i][j]->freq* height * 2.0;
                        //logo->letters[i][j]->scale = round(logo->letters[i][j]->scale* 100.0) / 100.0;
                        //fprintf(stdout,"%0.3f ", logo->letters[i][j]->scale);
                }
                //fprintf(stdout,"\n");

                qsort(logo->letters[i], 4,  sizeof(struct logo_letter*),sort_logo_letter_by_height);

        }
        run_draw_logo(logo,"text.png");

        free_logo(logo);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}
*/
int  run_draw_logo(struct logo_data* logo, char* outname)
{
        int i,j;
        cairo_surface_t *surface;
        cairo_t *cr;
        cairo_text_extents_t extents;
        cairo_font_extents_t font_extents;

        double base_w, base_h;

        surface = cairo_image_surface_create_for_data (image, CAIRO_FORMAT_ARGB32,
                                                       WIDTH, HEIGHT, STRIDE);

        cr = cairo_create (surface);
        cairo_select_font_face (cr, "monospace", 0, CAIRO_FONT_WEIGHT_BOLD);
        cairo_set_font_size (cr, BASE_FONT_SIZE);
        cairo_text_extents (cr, "C", &extents);
        cairo_font_extents(cr, &font_extents);
        base_h =  extents.height;
        base_w  =  extents.width + BASE_FONT_SIZE/6;
        cairo_destroy (cr);

        cairo_surface_destroy (surface);

        fprintf(stdout,"w:%f h:%f\n", base_w ,base_h);

        surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, (int)( base_w * logo->len),(int)( base_h* logo->L));
        fprintf(stdout,"w:%d h:%d\n",        cairo_image_surface_get_width(surface), cairo_image_surface_get_height(surface));

        cr = cairo_create (surface);

        cairo_set_source_rgb (cr, 0., 0., 0.);

        cairo_save (cr);





        fprintf(stdout,"%f w %f h\n", base_w, base_h);
        fprintf(stdout,"%f h (recommended offset)\n",font_extents.height);



        double y = 0.0;
        double x = 0.0;
        for(i = 0;i < logo->len;i++){
                y = 0.0;
                for(j = 0; j <  logo->L;j++){

                        write_letter(cr, x, y, logo->letters[i][j]->letter , logo->letters[i][j]->scale ,1);
                        y += base_h * logo->letters[i][j]->scale ;
                }
                 x+= base_w;
        }

        cairo_surface_write_to_png (surface, outname);

        cairo_destroy (cr);

        cairo_surface_destroy (surface);

        return 0;
}


int write_letter(cairo_t *cr, double x, double y,int c, double scale, int nuc)
{

        cairo_save (cr);

        char tmp[2];
        cairo_surface_t *surface;
        int height;
        double r,g,b;
        //LOG_MSG("Writing Letter");
        surface = cairo_get_group_target (cr);

        //fprintf(stdout,"w:%d h:%d\n",        cairo_image_surface_get_width(surface), cairo_image_surface_get_height(surface));

        height = cairo_image_surface_get_height(surface);

        scale= MACRO_MAX(0.01, scale);
        //fprintf(stdout,"Height: %f\n",y );
        cairo_select_font_face (cr, "monospace", 0, CAIRO_FONT_WEIGHT_BOLD);
        cairo_set_font_size (cr, BASE_FONT_SIZE);

        if(nuc){
                get_rgb_color(nuc_colors[c], &r,&g,&b);
                cairo_set_source_rgb (cr, r,g,b);
                tmp[0] = "ACGT"[c];
                tmp[1] = 0;

        }else{
                get_rgb_color(protein_colors[c], &r,&g,&b);
                cairo_set_source_rgb (cr, r,g,b);
                tmp[0] = "ACDEFGHIKLMNPQRSTVWY"[c];
                tmp[1] = 0;
        }




        cairo_move_to (cr, x,height -y);
        cairo_scale (cr,1.0,scale);

        cairo_show_text (cr, tmp);
        //fprintf(stdout,"Print letter:%s\n",tmp);
        cairo_restore (cr);
        return 0;
}

int sort_logo_letter_by_height(const void *a, const void *b)
{
        struct logo_letter* const *one = a;
        struct logo_letter* const *two = b;

        if((*one)->scale < (*two)->scale){
                return -1;
        }else{
                return 1;
        }
}
