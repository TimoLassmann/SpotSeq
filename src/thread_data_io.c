#include "tldevel.h"

#include "randomkit_io.h"

#define THREAD_DATA_IO_IMPORT
#include "thread_data_io.h"

#define BUFFER_LEN 256

struct seqer_thread_data** read_thread_data_to_hdf5(char* filename)
{
        struct seqer_thread_data** td = NULL;
        struct hdf5_data* hdf5_data = NULL;
        char buffer[BUFFER_LEN];
        unsigned int* seeds = NULL;
        int num_threads = 0;
        int max_len = 0;
        int max_K = 0;
        int i;

        RUN(open_hdf5_file(&hdf5_data,filename));

        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data, "/thread_data", "Nthreads",&num_threads));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data, "/thread_data", "Max_L",&max_len));
        RUN(HDFWRAP_READ_ATTRIBUTE(hdf5_data, "/thread_data", "Max_K",&max_K));
        //hdf5_open_file(filename,hdf5_data);

        /*hdf5_open_group("thread_data",hdf5_data);
        hdf5_read_attributes(hdf5_data, hdf5_data->group);
        ASSERT(hdf5_data->num_attr != 0 , "Could not find attributes");
        //print_attributes(hdf5_data);
        for(i = 0; i < hdf5_data->num_attr;i++){
                if(!strcmp("Nthreads", hdf5_data->attr[i]->attr_name)){
                        num_threads = hdf5_data->attr[i]->int_val;
                }
                if(!strcmp("Max_L", hdf5_data->attr[i]->attr_name)){
                        max_len = hdf5_data->attr[i]->int_val;
                }
                if(!strcmp("Max_K", hdf5_data->attr[i]->attr_name)){
                        max_K = hdf5_data->attr[i]->int_val;
                }
                }*/
        ASSERT(num_threads!=0, "No threads???");
        ASSERT(max_len!=0, "No len???");
        ASSERT(max_K!=0, "No states???");

        //hdf5_close_group(hdf5_data);
        //LOG_MSG("MAXLEN: %d",max_len);




        RUN(create_seqer_thread_data(&td,num_threads,  max_len, max_K, NULL));


        RUN(HDFWRAP_READ_DATA(hdf5_data, "/thread_data","seeds", &seeds));

        for(i = 0 ;i < num_threads;i++){
                snprintf(buffer, BUFFER_LEN , "/thread_data/RNG%d",i);
                //LOG_MSG("Trying to create group: %s", buffer);
                RUN(read_RNG_state(hdf5_data, buffer, &td[i]->rndstate));
                td[i]->seed = seeds[i];
        }

        RUN(close_hdf5_file(&hdf5_data));
        gfree(seeds);
        return td;
ERROR:
        return NULL;

}




int write_thread_data_to_hdf5(char* filename,struct seqer_thread_data** td,int num_threads,int max_len,int max_K)
{
        struct hdf5_data* hdf5_data = NULL;
        char buffer[BUFFER_LEN];
        unsigned int* seeds = NULL;
        int i;

        RUN(open_hdf5_file(&hdf5_data,filename));

        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data, "/thread_data", "Nthreads",num_threads));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data, "/thread_data", "Max_L",max_len));
        RUN(HDFWRAP_WRITE_ATTRIBUTE(hdf5_data, "/thread_data", "Max_K",max_K));
        /*
        RUNP(hdf5_data = hdf5_create());

        hdf5_open_file(filename,hdf5_data);


        hdf5_data->num_attr = 0;


        hdf5_add_attribute(hdf5_data, "Nthreads", "",num_threads, 0.0f, HDF5GLUE_INT);
        hdf5_add_attribute(hdf5_data, "Max_L", "",max_len, 0.0f, HDF5GLUE_INT);
        hdf5_add_attribute(hdf5_data, "Max_K", "",max_K, 0.0f, HDF5GLUE_INT);
        hdf5_create_group("thread_data",hdf5_data);
        hdf5_write_attributes(hdf5_data, hdf5_data->group);
        hdf5_data->num_attr = 0;

        hdf5_close_group(hdf5_data);*/


        RUN(galloc(&seeds, num_threads)); // MMALLOC(seeds , sizeof(unsigned int)* num_threads);
        for(i = 0; i < num_threads;i++){
                 seeds[i] = td[i]->seed;
        }


        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,"/thread_data","seeds",seeds));
        gfree(seeds);

        for(i = 0; i < num_threads;i++){
                snprintf(buffer, BUFFER_LEN , "/thread_data/RNG%d",i);
                //LOG_MSG("Trying to create group: %s", buffer);
                RUN(add_RNG_state(hdf5_data, buffer, &td[i]->rndstate));
        }

        RUN(close_hdf5_file(&hdf5_data));

        //hdf5_close_file(hdf5_data);
        //hdf5_free(hdf5_data);
        return OK;

ERROR:
        if(hdf5_data){
                RUN(close_hdf5_file(&hdf5_data));
                //hdf5_close_file(hdf5_data);
                //hdf5_free(hdf5_data);
        }
        return FAIL;
}
