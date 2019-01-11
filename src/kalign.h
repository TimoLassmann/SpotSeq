
#ifndef kalign_header
#define kalign_header

#define SEEK_START 0
#define SEEK_END 2

#define GPO 0.75
#define GPE 0.25
#define TGPE  0.0



#define NODESIZE 16

#define KALIGN_MAX(a, b) (a > b ? a : b)
#define KALIGN_MAX3(a,b,c) KALIGN_MAX(KALIGN_MAX(a,b),c)

#define CHUNK 0x4000

//unsigned char* write_int_to_string(unsigned char* p, int x, int* len, int* pos);
//int  read_int_from_string(unsigned char* p, int* pos);

extern char**  kalign_align(char** sequences, int numseq);

#endif















