#ifndef MOTIF_LOGO_H
#define MOTIF_LOGO_H

#define WIDTH 48
#define HEIGHT 142


#define STRIDE (WIDTH * 4)

#define BASE_FONT_SIZE 60

unsigned char image[STRIDE*HEIGHT];





struct logo_letter{
        double scale;
        double freq;
        char c;
        int letter;
        int gray;
};


struct logo_data{
        struct logo_letter*** letters;
        int len;
        int L;
};
extern int make_logo(int** matrix,int len, int L, char* outname);


#endif
