#include "tldevel.h"

#define RANDOMKIT_IO_IMPORT
#include "randomkit_io.h"



int add_RNG_state(struct hdf5_data* hdf5_data, char* group,rk_state* a)
{
        unsigned long* tmp_key = NULL;
        int i;
        RUN(galloc(&tmp_key,624));
        for(i = 0; i < 624;i++){
                tmp_key[i] = a->key[i];
        }


        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"key",tmp_key));
        gfree(tmp_key);

        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"gauss",a->gauss));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"psave",a->psave));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"has_binomial",a->has_binomial));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"has_gauss",a->has_gauss));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"pos",a->pos));
        RUN(HDFWRAP_WRITE_DATA(hdf5_data ,group,"nsave",a->nsave));
        return OK;

ERROR:
        return FAIL;
}

int read_RNG_state(struct hdf5_data* hdf5_data, char* group,rk_state* a)
{
        int i;
        uint64_t* tmp_key = NULL;

        RUN(HDFWRAP_READ_DATA(hdf5_data, group, "key", &tmp_key));

        for(i = 0; i < RK_STATE_LEN;i++){
                a->key[i] = tmp_key[i];
        }
        gfree(tmp_key);
        RUN(HDFWRAP_READ_DATA(hdf5_data ,group,"gauss",&a->gauss));
        RUN(HDFWRAP_READ_DATA(hdf5_data ,group,"psave",&a->psave));
        RUN(HDFWRAP_READ_DATA(hdf5_data ,group,"has_binomial",&a->has_binomial));
        RUN(HDFWRAP_READ_DATA(hdf5_data ,group,"has_gauss",&a->has_gauss));
        RUN(HDFWRAP_READ_DATA(hdf5_data ,group,"pos",&a->pos));
        RUN(HDFWRAP_READ_DATA(hdf5_data ,group,"nsave",&a->nsave));

        return OK;
ERROR:
        return FAIL;
}
