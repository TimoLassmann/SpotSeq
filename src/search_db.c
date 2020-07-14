#include "tldevel.h"
#include "tlrng.h"
#include "tlseqio.h"
#include "tlmisc.h"
#include "tlalphabet.h"
#include "tlhdf5wrap.h"

 #include <string.h>


#define SEARCH_DB_IMPORT
#include "search_db.h"

int alloc_hdf_seq_store(struct hdf_seq_store** hs);
int reset_hdf_seq_store(struct hdf_seq_store* h);
int resize_hdf_seq_store(struct hdf_seq_store* h);


int write_hdf_seq_store(struct hdf_seq_store* h, char* filename);
int read_hdf_seq_store_chunk(struct hdf_seq_store** hs, char* filename);

#ifdef SEQDBTEST
int main(int argc, char *argv[])
{
        if(argc == 3){
                build_sequence_database(argv[1], argv[2], 0);

        }else if(argc == 2){
                struct hdf_seq_store*h = NULL;
                int numseq = 0;
                while(1){
                        RUN(read_hdf_seq_store_chunk(&h, argv[1]));
                        if(!h->num_seq){
                                break;
                        }
                        numseq+= h->num_seq;
                }
                free_hdf_seq_store(h);
                LOG_MSG("Read %d", numseq);
        }
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}
#endif

int build_sequence_database(char* filename, char* out,int seed)
{
        struct file_handler* f = NULL;
        struct tl_seq_buffer* sb = NULL;


        struct rng_state* rng = NULL;
        struct alphabet* a = NULL;
        struct len_struct** len_arr = NULL;
        struct hdf_seq_store* hs = NULL;
        //float P_M, P_R;
        int ret;
        int i,j,len;

        int chunk;

        int num_items;
        int pos;
        uint8_t* seq;
        //LOG_MSG("%s",infill);
        if(!my_file_exists(filename)){
                ERROR_MSG("File %s not found");
        }

        if(!seed){
                RUNP(rng = init_rng(42));
        }else{
                RUNP(rng = init_rng(seed));
        }
        pos = 0;


        RUN(alloc_hdf_seq_store(&hs));
        RUN(reset_hdf_seq_store(hs));
        seq = hs->seq;
        RUN(open_fasta_fastq_file(&f, filename, TLSEQIO_READ));
        chunk = 1;
        while(1){
                RUN(read_fasta_fastq_file(f, &sb, 100000));
                if(chunk == 1){
                        if(sb->L == TL_SEQ_BUFFER_DNA){
                                RUN(create_alphabet(&a, rng, TLALPHABET_NOAMBIGUOUS_DNA));
                        }else if(sb->L == TL_SEQ_BUFFER_PROTEIN){
                                RUN(create_alphabet(&a, rng, TLALPHABET_NOAMBIGIOUS_PROTEIN ));
                        }
                }
                if(sb->num_seq == 0){
                        break;
                }

                LOG_MSG("Working on chunk: %d",chunk);
                while(hs->num_seq + sb->num_seq >= hs->alloc_num_seq){
                        RUN(resize_hdf_seq_store(hs));
                }

                for(i = 0; i < sb->num_seq;i++){
                        len = sb->sequences[i]->len;

                        /* Do I have to write? */
                        if(pos + len >= hs->seq_buf_size){

                                LOG_MSG("Need to write:    %d %d",pos + len, hs->seq_buf_size);
                                RUN(write_hdf_seq_store(hs,  out));
                                hs->chunk_number = hs->chunk_number+1;
                                reset_hdf_seq_store(hs);
                                pos = 0;
                        }
                        /* copyname  */

                        memccpy(hs->seq_names[hs->num_seq], sb->sequences[i]->name, 0, TL_SEQ_MAX_NAME_LEN);
                        /* set len */
                        hs->len[hs->num_seq] = pos;
                        /* convert */
                        convert_to_internal(a, (uint8_t*)sb->sequences[i]->seq,len);
                        for(j = 0; j < len;j++){
                                seq[pos] = sb->sequences[i]->seq[j];
                                pos++;
                        }
                        hs->num_seq++;
                        /* set end - this will be repeated in the next sequence   */
                        hs->len[hs->num_seq] = pos;
                }
                //hs->len[hs->num_seq] = pos;

                /*for(i = 0; i < hs->num_seq;i++){
                        fprintf(stdout,"%s\n",hs->seq_names[i]);
                        fprintf(stdout,"%s\n",sb->sequences[i]->name) ;
                        for(j = hs->len[i]; j < hs->len[i+1];j++){
                                fprintf(stdout,"%c",(char) hs->seq[j]);
                        }
                        fprintf(stdout,"\n");
                        fprintf(stdout,"%s\n",sb->sequences[i]->seq);

                        }*/

                //exit(0);
                chunk++;
        }
        /* final write */
        if(hs->num_seq){
                LOG_MSG("Need to write:    %d %d",pos + len, hs->seq_buf_size);
                RUN(write_hdf_seq_store(hs,  out));
                hs->chunk_number = hs->chunk_number+1;
                reset_hdf_seq_store(hs);
                pos = 0;

        }
        RUN(close_seq_file(&f));

        free_hdf_seq_store(hs);
        free_rng(rng);
        free_tl_seq_buffer(sb);
        free_alphabet(a);


        return OK;
ERROR:
        return FAIL;
}


int alloc_hdf_seq_store(struct hdf_seq_store** hs)
{
        struct hdf_seq_store* h = NULL;
        int i,j;
        MMALLOC(h, sizeof(struct hdf_seq_store));
        h->seq = NULL;
        h->len = NULL;
        h->seq_names = NULL;

        h->num_seq = 0;
        h->alloc_num_seq = 0;
        h->seq_buf_size = 0;

        h->seq_buf_size = 1073741824; /*  2** 30 - 1GB  */

        h->alloc_num_seq = 2000000;
        h->chunk_number = 1;
        //MMALLOC(h->seq, sizeof(uint8_t) * h->seq_buf_size);
        RUN(galloc(&h->seq, h->seq_buf_size));
        RUN(galloc(&h->len, (h->alloc_num_seq+1)));
        RUN(galloc(&h->seq_names, h->alloc_num_seq, TL_SEQ_MAX_NAME_LEN));
        for(i = 0; i < h->seq_buf_size;i++){
                h->seq[i] = 0;
        }

        for(i = 0; i < h->alloc_num_seq ;i++){
                h->len[i] = 0;
        }
        h->len[h->alloc_num_seq] = 0;



        for(i = 0; i < h->alloc_num_seq;i++){
                for(j = 0; j < TL_SEQ_MAX_NAME_LEN;j++){
                        h->seq_names[i][j] = 0;
                }
        }

        //MMALLOC(h->seq_names, size)

        *hs = h;
        return OK;
ERROR:
        free_hdf_seq_store(h);
        return FAIL;
}

int reset_hdf_seq_store(struct hdf_seq_store* h)
{
        h->num_seq = 0;
        return OK;
}

int resize_hdf_seq_store(struct hdf_seq_store* h)
{
        int old = 0;
        int i,j;
        old = h->alloc_num_seq;

        h->alloc_num_seq = h->alloc_num_seq + h->alloc_num_seq / 2;
        RUN(galloc(&h->len, (h->alloc_num_seq+1)));
        RUN(galloc(&h->seq_names, h->alloc_num_seq, TL_SEQ_MAX_NAME_LEN));

        for(i = old; i < h->alloc_num_seq ;i++){
                h->len[i] = 0;
        }
        h->len[h->alloc_num_seq] = 0;

        for(i = old; i < h->alloc_num_seq;i++){
                for(j = 0; j < TL_SEQ_MAX_NAME_LEN;j++){
                        h->seq_names[i][j] = 0;
                }
        }


        return OK;
ERROR:
        return FAIL;
}



void free_hdf_seq_store(struct hdf_seq_store* h)
{
        if(h){
                gfree(h->seq);
                gfree(h->len);
                gfree(h->seq_names);
                MFREE(h);
        }
}

int write_hdf_seq_store(struct hdf_seq_store* h, char* filename)
{
        struct hdf5_data* hdf5_data = NULL;
        char buffer[256];

        snprintf(buffer, 256, "/Chunk%d",h->chunk_number);

        RUN(open_hdf5_file(&hdf5_data,filename));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data,"/","NumChunks",h->chunk_number ));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data,buffer,"Chunk Number",h->chunk_number ));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data,buffer,"NumSeq",h->num_seq ));
        //LOG_MSG("p->fit: %p", p->fit);
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,buffer,"Seq", h->seq));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,buffer,"Len", h->len));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,buffer,"Names", h->seq_names));

        close_hdf5_file(&hdf5_data);

        return OK;
ERROR:
        if(hdf5_data){
                close_hdf5_file(&hdf5_data);

        }
        return FAIL;
}

int read_hdf_seq_store_chunk(struct hdf_seq_store** hs, char* filename)
{
        struct hdf_seq_store* h = NULL;
        struct hdf5_data* hdf5_data = NULL;
        char buffer[256];
        int t_chunks;
        if(!*hs){
                RUN(alloc_hdf_seq_store(&h));


        }else{
                h = *hs;
        }

        RUN(open_hdf5_file(&hdf5_data,filename));

        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data,"/","NumChunks",&t_chunks ));
        //LOG_MSG("%d %d",h->chunk_number, t_chunks);
        if(h->chunk_number  <=  t_chunks){
                snprintf(buffer, 256, "/Chunk%d",h->chunk_number);
                LOG_MSG("READING : %s", buffer);
                gfree(h->seq);
                gfree(h->len);
                gfree(h->seq_names);
                RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data,buffer,"NumSeq",&h->num_seq ));
                RUN(HDFWRAP_READ_DATA(hdf5_data ,buffer,"Seq", &h->seq));
                RUN(HDFWRAP_READ_DATA(hdf5_data ,buffer,"Len", &h->len));
                RUN(HDFWRAP_READ_DATA(hdf5_data ,buffer,"Names", &h->seq_names));

        }else{
                h->num_seq =0;
                //LOG_MSG("reading done;");
        }

        close_hdf5_file(&hdf5_data);

        h->chunk_number++;
        *hs = h;

        return OK;
ERROR:
        if(hdf5_data){
                close_hdf5_file(&hdf5_data);

        }
        return FAIL;
}
