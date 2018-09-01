#ifndef CALIBRATE_H
#define CALIBRATE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"


extern int calibrate(char* model_file,int num_threads, float* lambda, float* beta);

#endif
