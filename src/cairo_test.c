/*
 * Copyright © 2003 USC, Information Sciences Institute
 *
 * Permission to use, copy, modify, distribute, and sell this software
 * and its documentation for any purpose is hereby granted without
 * fee, provided that the above copyright notice appear in all copies
 * and that both that copyright notice and this permission notice
 * appear in supporting documentation, and that the name of the
 * University of Southern California not be used in advertising or
 * publicity pertaining to distribution of the software without
 * specific, written prior permission. The University of Southern
 * California makes no representations about the suitability of this
 * software for any purpose.  It is provided "as is" without express
 * or implied warranty.
 *
 * THE UNIVERSITY OF SOUTHERN CALIFORNIA DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE UNIVERSITY OF
 * SOUTHERN CALIFORNIA BE LIABLE FOR ANY SPECIAL, INDIRECT OR
 * CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
 * NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 * Author: Carl D. Worth <cworth@isi.edu>
 */

#include <cairo.h>
#include <stdio.h>
#include <math.h>


#define WIDTH 450
#define HEIGHT 140


#define STRIDE (WIDTH * 4)



unsigned char image[STRIDE*HEIGHT];
int write_letter(cairo_t *cr, double x, double* y,char* c, double scale);


int
main (void)
{
    int i;
    cairo_surface_t *surface;
    cairo_t *cr;
    cairo_text_extents_t extents;
    cairo_font_extents_t font_extents;



    double dx, dy;
    double height;


    surface = cairo_image_surface_create_for_data (image, CAIRO_FORMAT_ARGB32,
						   WIDTH, HEIGHT, STRIDE);


    cr = cairo_create (surface);

    cairo_set_source_rgb (cr, 0., 0., 0.);
    cairo_set_line_width (cr, 2.0);

    cairo_save (cr);


     double base_w, base_h;
        cairo_select_font_face (cr, "monospace", 0, 0);
        cairo_set_font_size (cr, 40);
        cairo_text_extents (cr, "C", &extents);
        base_h =  extents.height;
        base_w  =  extents.width;

        fprintf(stdout,"%f w %f h\n", base_h, base_w);
        base_w =  base_w + 5;
    double y = 0.0;

    write_letter(cr, 0.0, &y,"A", 1.0 );
    write_letter(cr, 0.0,  &y,"C", 1.0 );
    write_letter(cr, 0.0,  &y,"G", 1.0 );
    write_letter(cr, 0.0,  &y,"T", 1.0 );

    y = 0.0;
    write_letter(cr, base_w, &y,"A", 0.5 );
    write_letter(cr, base_w,  &y,"C", 0.5 );
    write_letter(cr, base_w,  &y,"G", 0.01 );// 0.01 is minimum !!!
    write_letter(cr, base_w,  &y,"T", 0.4 );

    cairo_surface_write_to_png (surface, "text.png");

    cairo_destroy (cr);

    cairo_surface_destroy (surface);

    return 0;
}


int write_letter(cairo_t *cr, double x, double* y,char* c, double scale)
{
        cairo_text_extents_t extents;
        double base_w, base_h;
        cairo_save (cr);
        fprintf(stdout,"Height: %f\n",*y );
        cairo_select_font_face (cr, "monospace", 0, 0);
        cairo_set_font_size (cr, 40);
        cairo_text_extents (cr, c, &extents);
        base_h =  extents.height;
        base_w  =  extents.width;

        cairo_set_source_rgb (cr, 1, 1, 0);

        cairo_move_to (cr, x, HEIGHT-*y);
        cairo_scale (cr,1.0,scale);
        cairo_show_text (cr, c);
        *y += base_h*scale +5;

        cairo_restore (cr);
}
