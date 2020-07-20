#include "motif_logo.h"
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

#include "global.h"
#define WIDTH 48
#define HEIGHT 142


#define STRIDE (WIDTH * 4)

#define BASE_FONT_SIZE 60

unsigned char image[STRIDE*HEIGHT];


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




int write_letter(cairo_t *cr, double x, double y,int letter_index, char letter, double scale, int nuc);
int write_letter_gray(cairo_t *cr, double x, double y,int letter_index, char letter, double scale, int nuc);
int  run_draw_logo(struct logo_data* logo, char* outname);
int get_rgb_color(int color, double*r, double *g, double *b);


static int sort_logo_letter_by_height(const void *a, const void *b);
struct logo_data* alloc_logo(int len, int L);
void free_logo(struct logo_data* logo);

int make_logo(int** matrix,int len, int L, char* outname)
{
        struct logo_data* logo = NULL;
        int i,j;

        RUNP(logo = alloc_logo(len, L));



        double total;
        double max;
        total = 0.0;
        max = -1.0;
        if(L ==  ALPHABET_DNA){
                for(i = 0; i < len;i++){
                        total = 0.0;
                        for(j = 0; j < L;j++){
                                logo->letters[i][j]->c = "ACGT"[j];
                                logo->letters[i][j]->freq = matrix[i][j];
                                        fprintf(stdout," %d", matrix[i][j]);
                                total += matrix[i][j];
                        }
                        fprintf(stdout,"\n");
                        max= MACRO_MAX(max, total);

                }
        }else if(L == ALPHABET_PROTEIN){
                for(i = 0; i < len;i++){
                        total = 0.0;
                        for(j = 0; j < L;j++){
                                logo->letters[i][j]->c =  "ACDEFGHIKLMNPQRSTVWY"[j];
                                logo->letters[i][j]->freq = matrix[i][j];
                                fprintf(stdout," %d", matrix[i][j]);
                                total += matrix[i][j];
                        }
                        fprintf(stdout,"\n");
                        max= MACRO_MAX(max, total);

                }
        }else{
                ERROR_MSG("Unknown alphabet! ");
        }

        //fprintf(stdout,"\n");
        double entropy, sum, e, height;

        e = (1.0 / log(2.0)) * ((double)(logo->L-1.0) / (2.0*max));

        //e = 0.0;

        for(i = 0; i < len;i++){
                sum = 0.0;
                for(j = 0; j < L;j++){

                        sum += logo->letters[i][j]->freq;
                }
                ASSERT(sum > 0.0, "empty column!");
                for(j = 0; j < L;j++){
                        logo->letters[i][j]->freq /= sum;
                        //fprintf(stdout,"%0.3f ", logo->letters[i][j]->freq);
                }
                entropy = 0.0;
                for(j = 0; j < L;j++){
                        if(logo->letters[i][j]->freq){
                                entropy += logo->letters[i][j]->freq *log2(logo->letters[i][j]->freq);
                        }
                }
                entropy = fabs(entropy);


                height =  log2((double) logo->L ) - (entropy + e);
                //fprintf(stdout,"pos:%d entropy: %f error:%f   height :%f max:%f L:%d\n",i,entropy,e, log2(4.0) - (entropy + e),max,logo->L );
                if(height < 0.25){
                        for(j = 0; j < L;j++){
                                logo->letters[i][j]->gray = 1;
                                logo->letters[i][j]->scale = logo->letters[i][j]->freq* 2.0 * log2(logo->L);
                        }
                }else{
                        for(j = 0; j < L;j++){
                                logo->letters[i][j]->scale = logo->letters[i][j]->freq* height * log2(logo->L);
                                //logo->letters[i][j]->scale = round(logo->letters[i][j]->scale* 100.0) / 100.0;
                                //fprintf(stdout,"%0.3f ", logo->letters[i][j]->scale);
                        }
                }
                //fprintf(stdout,"\n");

                qsort(logo->letters[i], L,  sizeof(struct logo_letter*),sort_logo_letter_by_height);

        }
        RUN(run_draw_logo(logo,outname));
        free_logo(logo);
        return OK;
ERROR:
        return FAIL;
}

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
                        logo->letters[i][j]->gray = 0;

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

int  run_draw_logo(struct logo_data* logo, char* outname)
{
        int i,j;
        cairo_surface_t *surface;
        cairo_t *cr;
        cairo_text_extents_t extents;
        cairo_font_extents_t font_extents;

        double base_w, base_h;
        double pic_height = 200.0;

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

        //fprintf(stdout,"w:%f h:%f\n", base_w ,base_h);

        //surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, (int)( base_w * logo->len),(int)( base_h* logo->L));


        surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, (int)( base_w * logo->len),(int)(pic_height ));


        //fprintf(stdout,"w:%d h:%d\n",        cairo_image_surface_get_width(surface), cairo_image_surface_get_height(surface));

        cr = cairo_create (surface);

        cairo_set_source_rgb (cr, 0., 0., 0.);

        cairo_save (cr);

        /* set scale modifier */

        double mod = pic_height / (base_h* (double) logo->L);




        //fprintf(stdout,"%f w %f h\n", base_w, base_h);
        //fprintf(stdout,"%f h (recommended offset)\n",font_extents.height);



        double y = 0.0;
        double x = 0.0;
        for(i = 0;i < logo->len;i++){
                y = 0.0;
                for(j = 0; j < logo->L;j++){
                        if(logo->letters[i][j]->gray){
                                if(logo->L == ALPHABET_DNA){
                                        write_letter_gray(cr, x, y, logo->letters[i][j]->letter,logo->letters[i][j]->c , logo->letters[i][j]->scale*mod ,1);
                                }else{
                                        write_letter_gray(cr, x, y, logo->letters[i][j]->letter,logo->letters[i][j]->c , logo->letters[i][j]->scale*mod ,0);
                                }
                        }else{
                                if(logo->L == ALPHABET_DNA){
                                        write_letter(cr, x, y, logo->letters[i][j]->letter ,logo->letters[i][j]->c, logo->letters[i][j]->scale*mod ,1);
                                }else{
                                        write_letter(cr, x, y, logo->letters[i][j]->letter ,logo->letters[i][j]->c, logo->letters[i][j]->scale*mod ,0);
                                }

                        }
                        y += base_h * logo->letters[i][j]->scale*mod ;
                }
                x+= base_w;
        }
        LOG_MSG("Writing to:%s", outname);
        cairo_surface_write_to_png (surface, outname);

        cairo_destroy (cr);

        cairo_surface_destroy (surface);

        return OK;
}

int write_letter_gray(cairo_t *cr, double x, double y,int letter_index, char letter,double scale, int nuc)
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
        r = 192.0 / 255.0;
        g = 192.0 / 255.0;
        b = 192.0 / 255.0;
        cairo_set_source_rgb (cr, r,g,b);
        tmp[0] = letter;
        tmp[1] = 0;





        cairo_move_to (cr, x,height -y);
        cairo_scale (cr,1.0,scale);

        cairo_show_text (cr, tmp);
        //fprintf(stdout,"Print letter:%s\n",tmp);
        cairo_restore (cr);
        return 0;
}

int write_letter(cairo_t *cr, double x, double y,int letter_index, char letter, double scale, int nuc)
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
                get_rgb_color(nuc_colors[letter_index], &r,&g,&b);
                cairo_set_source_rgb (cr, r,g,b);


        }else{
                get_rgb_color(protein_colors[letter_index], &r,&g,&b);
                cairo_set_source_rgb (cr, r,g,b);

        }


        tmp[0] = letter;// "ACGT"[c];
        tmp[1] = 0;

        cairo_move_to (cr, x,height -y);
        cairo_scale (cr,1.0,scale);

        cairo_show_text (cr, tmp);
        //fprintf(stdout,"Print letter:%s\n",tmp);
        cairo_restore (cr);
        return OK;
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



int get_rgb_color(int color, double*r, double *g, double *b)
{


        *r = (double)((color >> 16) & 0xFF) / (double) 0xFF;
        *g = (double)((color >> 8) & 0xFF) / (double) 0xFF;
        *b = (double)(color & 0xFF) / (double) 0xFF;
        return OK;
}
