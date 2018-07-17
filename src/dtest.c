#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>

#include "tldevel.h"
#include "distributions.h"


int main (int argc, char *argv[]) 
{
        float emission[100][4];
        float tmp[4];
        float var[4];
        rk_state rndstate;
        float alpha, sum;
        int i,j,c;


        rk_randomseed(&rndstate);

        for(i = 0; i < 1000;i++){
                for(c = 0;c < 4;c++){
                        tmp[c] = 0.0f;
                }
                alpha = (float) (i+1) / 100.0f;
                for(j = 0; j < 100;j++){
                        sum = 0.0f;
                        while(!sum){
                                sum = 0.0f;
                                for(c = 0;c < 4;c++){
                                        emission[j][c] = rk_gamma(&rndstate,alpha , 1.0);
                                        sum += emission[j][c];
                                }
                                for(c = 0;c < 4;c++){
                                        emission[j][c] /= sum;
                                        tmp[c] += emission[j][c];
                                }
                        }
                        
                }
                for(c = 0;c < 4;c++){
                        tmp[c] /= 100.0f;
                        var[c] = 0.0f;
                }
                for(j = 0; j < 100;j++){
                        for(c = 0;c < 4;c++){
                                var[c] += powf(emission[j][c] - tmp[c], 2.0f);
                        }
                }
                for(c = 0;c < 4;c++){
                        var[c] = sqrtf( var[c] / 99.0f );
                }
                fprintf(stdout,"alpha:%f",alpha);
                for(c = 0;c < 4;c++){
                        fprintf(stdout,"\t%f",tmp[c]);
                }
                for(c = 0;c < 4;c++){
                        fprintf(stdout,"\t%f",var[c]);
                }
                fprintf(stdout,"\n");
        }
         
        return EXIT_SUCCESS;
}
