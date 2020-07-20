#ifndef MOTIF_LOGO_H
#define MOTIF_LOGO_H



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
