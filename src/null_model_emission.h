#ifndef NULL_MODEL_EMISSION_H
#define NULL_MODEL_EMISSION_H



#ifdef NULL_MODEL_EMISSION_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int get_null_model_emissions(double** b, int L);

#undef NULL_MODEL_EMISSION_IMPORT
#undef EXTERN

#endif
