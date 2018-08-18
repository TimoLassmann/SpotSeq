#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"


#include "hdf5_glue.h"
#include <string.h>
#include <stdio.h>


static int hdf5_create_dataset_compact(char* dataset_name,struct hdf5_data* hdf5_data, hid_t type);
static int hdf5_close_dataset(struct hdf5_data* hdf5_data);

static herr_t op_func (hid_t loc_id, const char *name, const H5L_info_t *info, void *operator_data);

int get_group_names(struct hdf5_data* hdf5_data)
{
        //struct hdf5_group_names* grp_names = NULL;
        //int status;
        //struct hdf5_group_names* ptr = NULL;

        MMALLOC(hdf5_data->grp_names, sizeof(struct hdf5_group_names));


        //MMALLOC(hdf5_data->grp_names, sizeof(struct hdf5_group_names));

        hdf5_data->grp_names->names = NULL;
        hdf5_data->grp_names->num_names = 0;
        hdf5_data->grp_names->name_length = 0;
        hdf5_data->grp_names->num_names_mem = 65536;

        RUNP(hdf5_data->grp_names->names = malloc_2d_char(hdf5_data->grp_names->names, hdf5_data->grp_names->num_names_mem , HDF5GLUE_MAX_NAME_LEN,0));
        //hdf5_data->grp_names->names = malloc_2d_char(hdf5_data->grp_names->names, hdf5_data->grp_names->num_names_mem , HDF5GLUE_MAX_NAME_LEN,0);

        if((hdf5_data->status = H5Literate (hdf5_data->file, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, op_func, hdf5_data->grp_names)) < 0) ERROR_MSG("H5Literate failed");

        return OK;
ERROR:

        MFREE(hdf5_data->grp_names);
        return FAIL;
}

herr_t op_func (hid_t loc_id, const char *name, const H5L_info_t *info, void *operator_data)
{
        herr_t status;
        H5O_info_t infobuf;
        struct hdf5_group_names* grp_names = (struct hdf5_group_names*)operator_data;

        status = H5Oget_info_by_name (loc_id, name, &infobuf, H5P_DEFAULT);
        if(status <0){
                ERROR_MSG(" H5Oget_info_by_name failed." );
        }

        switch (infobuf.type){
        case H5O_TYPE_GROUP:
                snprintf(grp_names->names[grp_names->num_names], HDF5GLUE_MAX_NAME_LEN, "%s",name);
                grp_names->num_names++;
                if(grp_names->num_names ==grp_names->num_names_mem){
                        grp_names->num_names_mem = grp_names->num_names_mem << 1;
                        grp_names->names = malloc_2d_char(grp_names->names, grp_names->num_names_mem , HDF5GLUE_MAX_NAME_LEN,0);
                }
                break;
        case H5O_TYPE_DATASET:
                //printf ("  Dataset: %s\n", name);
                break;
        case H5O_TYPE_NAMED_DATATYPE:
                //printf ("  Datatype: %s\n", name);
                break;
        case H5O_TYPE_UNKNOWN:
        case H5O_TYPE_NTYPES:
                break;
        }
        return OK;
ERROR:
        return FAIL;
}

struct hdf5_data* hdf5_create(void)
{
        //int status;
        int i;
        struct hdf5_data* hdf5_data = NULL;


        MMALLOC(hdf5_data, sizeof(struct hdf5_data));

        hdf5_data->attr = NULL;
        hdf5_data->num_attr = 0;
        hdf5_data->num_attr_mem = 100;

        MMALLOC(hdf5_data->attr , sizeof(struct hdf5_attribute*) * hdf5_data->num_attr_mem);


        for(i = 0; i < hdf5_data->num_attr_mem;i++){
                hdf5_data->attr[i] = NULL;
                MMALLOC(hdf5_data->attr[i],  sizeof(struct hdf5_attribute));
        }

        for(i = 0; i < HDF5GLUE_MAX_DIM;i++){
                hdf5_data->dim[i] = 0;
                hdf5_data->chunk_dim[i] = 0;
        }

        hdf5_data->grp_names = NULL;
        hdf5_data->data = 0;

        hdf5_data->fapl = 0;
        hdf5_data->file = 0;
        hdf5_data->group = 0;

        hdf5_data->plist = 0;
        hdf5_data->dataset = 0;

        hdf5_data->attribute_id = 0;
        hdf5_data->attr_dataspace_id = 0;


        hdf5_data->datatype = 0;
        hdf5_data->dataspace = 0;

        hdf5_data->status = 0;

        hdf5_data->rank = 0;

        return hdf5_data;
ERROR:
        if(hdf5_data){
                for(i = 0; i < hdf5_data->num_attr_mem;i++){
                        MFREE(hdf5_data->attr[i]);
                }
                MFREE(hdf5_data->attr);
        }
        MFREE(hdf5_data );

        return NULL;
}


int hdf5_add_attribute(struct hdf5_data* hdf5_data,char* attr_name, char* string, int int_val, float float_val, int type)
{
        int n = hdf5_data->num_attr;

        ASSERT(n != hdf5_data->num_attr_mem,"No space to add more attributes...");
        //if(n == hdf5_data->num_attr_mem){
        //	KSLIB_EXCEPTION(kslFAIL,"No space to add more attributes...");
        //}

        snprintf(hdf5_data->attr[n]->attr_name, HDF5GLUE_MAX_NAME_LEN,"%s", attr_name);

        hdf5_data->attr[n]->type = type;

        if(string){
                snprintf(hdf5_data->attr[n]->string, HDF5GLUE_MAX_CONTENT_LEN,"%s", string);
        }
        hdf5_data->attr[n]->int_val = int_val;
        hdf5_data->attr[n]->float_val = float_val;


        hdf5_data->num_attr++;

        return OK;
ERROR:
        return FAIL;
}


int hdf5_read_attributes(struct hdf5_data* hdf5_data,hid_t target)
{
        H5O_info_t oinfo;
        hid_t atype,atype_mem;
        H5T_class_t type_class;
        int i;
        hdf5_data->status = H5Oget_info(target, &oinfo);
        hdf5_data->num_attr = 0;
        for(i = 0; i < (unsigned)oinfo.num_attrs; i++) {

                if(i == hdf5_data->num_attr_mem){
                        break;
                }

                hdf5_data->attribute_id = H5Aopen_by_idx(target, ".", H5_INDEX_CRT_ORDER, H5_ITER_INC, (hsize_t)i, H5P_DEFAULT, H5P_DEFAULT);
                atype = H5Aget_type(hdf5_data->attribute_id);
                type_class = H5Tget_class(atype);
                H5Aget_name(hdf5_data->attribute_id , HDF5GLUE_MAX_NAME_LEN  , hdf5_data->attr[i]->attr_name);

                if (type_class == H5T_STRING) {
                        hdf5_data->attr[i]->type = HDF5GLUE_CHAR;
                        atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);
                        hdf5_data->status   = H5Aread(hdf5_data->attribute_id, atype_mem, hdf5_data->attr[i]->string  );
                        hdf5_data->status   = H5Tclose(atype_mem);
                }

                if(type_class == H5T_INTEGER){
                        hdf5_data->attr[i]->type = HDF5GLUE_INT;
                        hdf5_data->status = H5Aread(hdf5_data->attribute_id, H5T_NATIVE_INT, &hdf5_data->attr[i]->int_val);
                }
                if(type_class == H5T_FLOAT){
                        hdf5_data->attr[i]->type = HDF5GLUE_FLOAT;
                        hdf5_data->status = H5Aread(hdf5_data->attribute_id, H5T_NATIVE_FLOAT, &hdf5_data->attr[i]->float_val);
                }

                hdf5_data->num_attr++;


                hdf5_data->status   = H5Aclose(hdf5_data->attribute_id);
                hdf5_data->status   = H5Tclose(atype);
        }
        return OK;
}

int hdf5_write_attributes(struct hdf5_data* hdf5_data,hid_t target)
{
        int i;
        hid_t aid;
        hid_t atype;
        hid_t attr;


        int len;

        for(i = 0;i < hdf5_data->num_attr;i++){
                switch (hdf5_data->attr[i]->type) {
                case HDF5GLUE_CHAR:

                        len = (int) strlen(hdf5_data->attr[i]->string) + 1;

                        aid  = H5Screate(H5S_SCALAR);
                        atype = H5Tcopy(H5T_C_S1);
                        H5Tset_size(atype, len);
                        H5Tset_strpad(atype,H5T_STR_NULLTERM);
                        attr = H5Acreate2(target,hdf5_data->attr[i]->attr_name, atype, aid, H5P_DEFAULT, H5P_DEFAULT);


                        hdf5_data->status = H5Awrite(attr, atype, hdf5_data->attr[i]->string);

                        hdf5_data->status = H5Sclose(aid);
                        hdf5_data->status = H5Tclose(atype);
                        hdf5_data->status = H5Aclose(attr);
                        break;
                case HDF5GLUE_INT:

                        aid  = H5Screate(H5S_SCALAR);
                        attr = H5Acreate2(target,hdf5_data->attr[i]->attr_name, H5T_NATIVE_INT, aid,  H5P_DEFAULT, H5P_DEFAULT);

                        hdf5_data->status = H5Awrite(attr, H5T_NATIVE_INT, &hdf5_data->attr[i]->int_val);
                        hdf5_data->status = H5Sclose(aid);
                        hdf5_data->status = H5Aclose(attr);

                        break;
                case HDF5GLUE_FLOAT:
                        aid  = H5Screate(H5S_SCALAR);
                        attr = H5Acreate2(target,hdf5_data->attr[i]->attr_name, H5T_NATIVE_FLOAT, aid,  H5P_DEFAULT, H5P_DEFAULT);

                        hdf5_data->status = H5Awrite(attr, H5T_NATIVE_FLOAT, &hdf5_data->attr[i]->float_val);
                        hdf5_data->status = H5Sclose(aid);
                        hdf5_data->status = H5Aclose(attr);
                        break;
                default:
                        break;
                }
        }
        hdf5_data->num_attr = 0;
        return OK;
}


int hdf5_free(struct hdf5_data* hdf5_data)
{

        int i;
        for(i = 0; i < hdf5_data->num_attr_mem;i++){
                //hdf5_data->attr = NULL;
                MFREE(hdf5_data->attr[i]);//, sizeof(struct hdf5_attribute));

        }
        MFREE(hdf5_data->attr);// , sizeof(struct hdf5_attribute*) * hdf5_data->num_attr_mem);

        if(hdf5_data->grp_names){
                free_2d((void**) hdf5_data->grp_names->names);
                MFREE(hdf5_data->grp_names);
        }

        MFREE(hdf5_data);
        return OK;
}


int hdf5_open_file(char* filename,struct hdf5_data* hdf5_data)
{

        if((hdf5_data->file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT)) == -1)ERROR_MSG("H5Fopen failed");



        //hdf5_read_attributes(hdf5_data, hdf5_data->file);
        return OK;
ERROR:
        return FAIL;
}

int hdf5_open_group(char* groupname,struct hdf5_data* hdf5_data)
{
        if((hdf5_data->group = H5Gopen2(hdf5_data->file, groupname , H5P_DEFAULT)) == -1)ERROR_MSG("H5Gopen2 failed");
        //hdf5_read_attributes(hdf5_data, hdf5_data->group);
        return OK;
ERROR:
        return FAIL;
}


int hdf5_read_dataset(char* dataset_name,struct hdf5_data* hdf5_data)
{
        //int i;
        int type = -1;
        char* m1d_char = NULL;
        int* m1d_int = NULL;

        int* m1d_ushort_int = NULL;

        float* m1d_float = NULL;
        double* m1d_double = NULL;

        char** m2d_char = NULL;
        int** m2d_int = NULL;
        float** m2d_float = NULL;
        float*** m3d_float = NULL;
        float**** m4d_float = NULL;
        double** m2d_double = NULL;

        if((hdf5_data->dataset = H5Dopen(hdf5_data->group, dataset_name,H5P_DEFAULT)) == -1)ERROR_MSG("H5Dopen failed\n");
        //printf ("H5Dopen returns: %d\n", hdf5_data->dataset);
        hdf5_read_attributes(hdf5_data,hdf5_data->dataset);


        /* Get datatype and dataspace handles and then query






           dataset class, order, size, rank and dimensions  */
        hdf5_data->datatype  = H5Dget_type(hdf5_data->dataset );     /* datatype handle */
        //printf ("H5Dget_type returns: %d\n", hdf5_data->datatype);

        //class     = H5Tget_class(hdf5_data->datatype);


        //hid_t native_type = H5Tget_native_type(datatype, H5T_DIR_ASCEND);

        if(H5Tequal(hdf5_data->datatype ,H5T_NATIVE_CHAR) != 0){
                type =HDF5GLUE_CHAR;
        }else if(H5Tequal(hdf5_data->datatype ,H5T_NATIVE_INT) != 0){
                type =HDF5GLUE_INT;
        }else if(H5Tequal(hdf5_data->datatype ,H5T_NATIVE_USHORT) != 0){
                type = HDF5GLUE_USHORT_INT;
        }else if(H5Tequal(hdf5_data->datatype ,H5T_NATIVE_FLOAT) != 0){
                type =HDF5GLUE_FLOAT;
        }else if(H5Tequal(hdf5_data->datatype ,H5T_NATIVE_DOUBLE) != 0){
                type =HDF5GLUE_DOUBLE;
        }

        DPRINTF3( "TYPE is %d", type);


        hdf5_data->dataspace = H5Dget_space(hdf5_data->dataset);
        hdf5_data->rank      = H5Sget_simple_extent_ndims(hdf5_data->dataspace);
        hdf5_data->status  = H5Sget_simple_extent_dims(hdf5_data->dataspace,hdf5_data->dim , NULL);

        hdf5_data->data = NULL;




        if(type == HDF5GLUE_CHAR){
                if(hdf5_data->rank == 1){
                        MMALLOC(m1d_char,sizeof(char)* hdf5_data->dim[0]);
                        hdf5_data->status = H5Dread(hdf5_data->dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m1d_char[0]);
                        hdf5_data->data = m1d_char;
                }
                if(hdf5_data->rank == 2){
                        m2d_char = malloc_2d_char(m2d_char,(int)hdf5_data->dim[0],(int)hdf5_data->dim[1],0);
                        hdf5_data->status = H5Dread(hdf5_data->dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m2d_char[0][0]);
                        hdf5_data->data = m2d_char;
                }
        }

        if(type == HDF5GLUE_INT){
                if(hdf5_data->rank == 1){
                        MMALLOC(m1d_int,sizeof(int)* hdf5_data->dim[0]);
                        hdf5_data->status = H5Dread(hdf5_data->dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m1d_int[0]);

                        hdf5_data->data = m1d_int;
                }
                if(hdf5_data->rank == 2){


                        m2d_int = malloc_2d_int(m2d_int,(int)hdf5_data->dim[0],(int)hdf5_data->dim[1],0);
                        hdf5_data->status = H5Dread(hdf5_data->dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m2d_int[0][0]);

                        hdf5_data->data = m2d_int;

                }
        }

        if(type == HDF5GLUE_USHORT_INT){
                if(hdf5_data->rank == 1){
                        //fprintf(stdout,"%llu allocing\n", sizeof(unsigned short int)* hdf5_data->dim[0]);
                        MMALLOC(m1d_ushort_int,sizeof(unsigned short int)* hdf5_data->dim[0]);
                        //fprintf(stdout,"Done\n");//, sizeof(unsigned short int)* hdf5_data->dim[0]);
                        hdf5_data->status = H5Dread(hdf5_data->dataset, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m1d_ushort_int[0]);

                        hdf5_data->data = m1d_ushort_int;
                }

        }

        if(type == HDF5GLUE_FLOAT){
                if(hdf5_data->rank == 1){
                        MMALLOC(m1d_float,sizeof(float)* hdf5_data->dim[0]);
                        hdf5_data->status = H5Dread(hdf5_data->dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m1d_float[0]);

                        hdf5_data->data = m1d_float;
                }
                if(hdf5_data->rank == 2){


                        m2d_float = malloc_2d_float(m2d_float,(int)hdf5_data->dim[0],(int)hdf5_data->dim[1],0.0f);
                        hdf5_data->status = H5Dread(hdf5_data->dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m2d_float[0][0]);

                        hdf5_data->data = m2d_float;

                }

                if(hdf5_data->rank == 3){


                        m3d_float = malloc_3d_float((int)hdf5_data->dim[0],(int)hdf5_data->dim[1],(int)hdf5_data->dim[2],0.0f);
                        hdf5_data->status = H5Dread(hdf5_data->dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m3d_float[0][0][0]);

                        hdf5_data->data = m3d_float;

                }

                if(hdf5_data->rank == 4){


                        m4d_float = malloc_4d_float((int)hdf5_data->dim[0],(int)hdf5_data->dim[1],(int)hdf5_data->dim[2],(int)hdf5_data->dim[3],0.0f);
                        hdf5_data->status = H5Dread(hdf5_data->dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m4d_float[0][0][0][0]);

                        hdf5_data->data = m4d_float;

                }


        }

        if(type == HDF5GLUE_DOUBLE){
                if(hdf5_data->rank == 1){
                        MMALLOC(m1d_double,sizeof(double)* hdf5_data->dim[0]);
                        hdf5_data->status = H5Dread(hdf5_data->dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m1d_double[0]);

                        hdf5_data->data = m1d_double;
                }
                if(hdf5_data->rank == 2){
                        m2d_double = malloc_2d_double(m2d_double,(int)hdf5_data->dim[0],(int)hdf5_data->dim[1],0.0);
                        hdf5_data->status = H5Dread(hdf5_data->dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m2d_double[0][0]);

                        hdf5_data->data = m2d_double;
                }
        }




        if((hdf5_data->status = H5Tclose(hdf5_data->datatype)) < 0) ERROR_MSG("H5Tclose failed");
        if((hdf5_data->status = H5Dclose(hdf5_data->dataset)) < 0) ERROR_MSG("H5Dclose failed");

        return OK;
ERROR:

        MFREE(m1d_char);

        MFREE(m1d_double);

        MFREE(m1d_float);

        MFREE(m1d_int);

        free_2d((void**)m2d_char);

        free_2d((void**)m2d_int);

        free_2d((void**)m2d_float);

        free_2d((void**)m2d_double);

        free_3d((void***)m3d_float);

        free_4d((void****)m4d_float);


        return FAIL;
}




int hdf5_create_file(char* filename,struct hdf5_data* hdf5_data)
{
        snprintf(hdf5_data->file_name , HDF5GLUE_MAX_NAME_LEN,"%s",filename);

        //if((hdf5_data->fapl = H5Pcreate (H5P_FILE_ACCESS)) == -1)ERROR_MSG("H5Pcreate failed\n");


        //if((hdf5_data->status = H5Pset_libver_bounds (hdf5_data->fapl , H5F_LIBVER_LATEST, H5F_LIBVER_LATEST)) < 0)ERROR_MSG("H5Pset_libver_bounds failed\n");

        if((hdf5_data->file = H5Fcreate(hdf5_data->file_name, H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT)) < 0)  ERROR_MSG("H5Fcreate failed: %s\n",hdf5_data->file_name);



        //if(hdf5_data->num_attr){
        //	if((hdf5_data->status = hdf5_write_attribute(hdf5_data, hdf5_data->file)) != kslOK)KSLIB_EXCEPTION(kslFAIL,"hdf5_write_attribute failed\n");
        //}




        return OK;
ERROR:
        return FAIL;
}

int hdf5_close_group(struct hdf5_data* hdf5_data)
{
        if((hdf5_data->status = H5Gclose(hdf5_data->group)) < 0) ERROR_MSG("H5Gclose failed");

        return OK;
ERROR:
        return FAIL;
}

int hdf5_close_file(struct hdf5_data* hdf5_data)
{
        if(hdf5_data->fapl){
                if((hdf5_data->status = H5Pclose(hdf5_data->fapl)) < 0) ERROR_MSG("H5Pclose failed");
        }
        if(hdf5_data->file){
                if((hdf5_data->status = H5Fclose(hdf5_data->file)) < 0) ERROR_MSG("H5Fclose failed");
        }
        return OK;
ERROR:
        return FAIL;
}


int hdf5_create_group(char* groupname,struct hdf5_data* hdf5_data)
{
        snprintf(hdf5_data->group_name , HDF5GLUE_MAX_NAME_LEN,"%s",groupname);

        hdf5_data->status = H5Eset_auto(hdf5_data->status ,NULL, NULL);
        hdf5_data->status = H5Gget_objinfo (hdf5_data->file, hdf5_data->group_name, 0, NULL);

        if (hdf5_data->status == 0){
                //printf ("The group %s exists.\n",hdf5_data->group_name );
                if((hdf5_data->group = H5Gopen2(hdf5_data->file,  hdf5_data->group_name , H5P_DEFAULT)) == -1)ERROR_MSG("H5Gopen2 failed\n");
        }else{
                //printf ("The group %s either does NOT exist\n or some other error occurred.\n",hdf5_data->group_name);
                if((hdf5_data->group = H5Gcreate (hdf5_data->file , hdf5_data->group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0)ERROR_MSG("H5Gcreate failed\n");
        }

        return OK;
ERROR:
        return FAIL;
}

/* this function should create an empty data group that is expandable.*/

int hdf5_create_dataset_chunked_expandable(char* dataset_name,struct hdf5_data* hdf5_data, hid_t type)
{

        hsize_t max_dim[HDF5GLUE_MAX_DIM];
        int i;
        for(i = 0; i< hdf5_data->rank;i++){
                hdf5_data->chunk_dim[i] = 256;

                max_dim[i] = H5S_UNLIMITED;
                //set dim to zero...
                hdf5_data->dim[i] = 0;

        }
        snprintf(hdf5_data->dataset_name , HDF5GLUE_MAX_NAME_LEN,"%s",dataset_name);
        if((hdf5_data->dataspace = H5Screate_simple(hdf5_data->rank,  hdf5_data->dim , max_dim)) < 0)ERROR_MSG("H5Screate_simple failed.");

        if((hdf5_data->datatype = H5Tcopy(type)) < 0) ERROR_MSG("H5Tcopy failed");

        if((hdf5_data->status = H5Tset_order(hdf5_data->datatype, H5T_ORDER_LE)) < 0) ERROR_MSG("H5Tset_order failed.");

        if((hdf5_data->plist = H5Pcreate (H5P_DATASET_CREATE)) < 0) ERROR_MSG("H5Pcreate failed.");

        if((hdf5_data->status = H5Pset_shuffle (hdf5_data->plist )) < 0 )ERROR_MSG("H5Pset_shuffle failed.");

        if((hdf5_data->status = H5Pset_deflate (hdf5_data->plist, 9)) < 0 )ERROR_MSG("H5Pset_deflate failed.");
        if((hdf5_data->status = H5Pset_chunk (hdf5_data->plist, hdf5_data->rank,  hdf5_data->chunk_dim)) < 0 )ERROR_MSG("H5Pset_chunk failed.");


        if((hdf5_data->dataset = H5Dcreate(hdf5_data->group, hdf5_data->dataset_name, hdf5_data->datatype, hdf5_data->dataspace,    H5P_DEFAULT, hdf5_data->plist , H5P_DEFAULT)) < 0 )ERROR_MSG("H5Dcreate failed");
        //hdf5_data->plist
        if(hdf5_data->num_attr){

                RUN(hdf5_write_attributes(hdf5_data, hdf5_data->dataset));
        }


        hdf5_close_dataset(hdf5_data);
        hdf5_close_group(hdf5_data);

        return OK;
ERROR:
        return FAIL;


}

int hdf5_create_dataset_chunked(char* dataset_name,struct hdf5_data* hdf5_data, hid_t type)
{
        int i;


        for(i = 0; i< hdf5_data->rank;i++){
                if(hdf5_data->chunk_dim[i] >  hdf5_data->dim[i]){
                        ERROR_MSG("chunk dimenson exceeds dataset dimension:%d (rank) %d %d \n", i,hdf5_data->chunk_dim[i], hdf5_data->dim[i] );
                }
                if(hdf5_data->chunk_dim[i] == 0){
                        hdf5_data->chunk_dim[i] = hdf5_data->dim[i];
                }
        }


        snprintf(hdf5_data->dataset_name , HDF5GLUE_MAX_NAME_LEN,"%s",dataset_name);
        if((hdf5_data->dataspace = H5Screate_simple(hdf5_data->rank,  hdf5_data->dim , NULL)) < 0)ERROR_MSG("H5Screate_simple failed.");

        if((hdf5_data->datatype = H5Tcopy(type)) < 0) ERROR_MSG("H5Tcopy failed");

        if((hdf5_data->status = H5Tset_order(hdf5_data->datatype, H5T_ORDER_LE)) < 0) ERROR_MSG("H5Tset_order failed.");

        if((hdf5_data->plist = H5Pcreate (H5P_DATASET_CREATE)) < 0) ERROR_MSG("H5Pcreate failed.");

        if((hdf5_data->status = H5Pset_shuffle (hdf5_data->plist )) < 0 )ERROR_MSG("H5Pset_shuffle failed.");

        if((hdf5_data->status = H5Pset_deflate (hdf5_data->plist, 2)) < 0 )ERROR_MSG("H5Pset_deflate failed.");
        if((hdf5_data->status = H5Pset_chunk (hdf5_data->plist, hdf5_data->rank,  hdf5_data->chunk_dim)) < 0 )ERROR_MSG("H5Pset_chunk failed.");


        if((hdf5_data->dataset = H5Dcreate(hdf5_data->group, hdf5_data->dataset_name, hdf5_data->datatype, hdf5_data->dataspace,    H5P_DEFAULT, hdf5_data->plist , H5P_DEFAULT)) < 0 )ERROR_MSG("H5Dcreate failed");
        //hdf5_data->plist
        if(hdf5_data->num_attr){

                RUN(hdf5_write_attributes(hdf5_data, hdf5_data->dataset));
        }

        return OK;
ERROR:
        return FAIL;
}



static int hdf5_create_dataset_compact(char* dataset_name,struct hdf5_data* hdf5_data, hid_t type)
{
        int i;


        for(i = 0; i< hdf5_data->rank;i++){
                if(hdf5_data->chunk_dim[i] >  hdf5_data->dim[i]){
                        ERROR_MSG("chunk dimenson exceeds dataset dimension:%d (rank) %d %d \n", i,hdf5_data->chunk_dim[i], hdf5_data->dim[i] );
                }
                if(hdf5_data->chunk_dim[i] == 0){
                        hdf5_data->chunk_dim[i] = hdf5_data->dim[i];
                }
        }


        snprintf(hdf5_data->dataset_name , HDF5GLUE_MAX_NAME_LEN,"%s",dataset_name);
        if((hdf5_data->dataspace = H5Screate_simple(hdf5_data->rank,  hdf5_data->dim , NULL)) < 0)ERROR_MSG("H5Screate_simple failed\n");

        if((hdf5_data->datatype = H5Tcopy(type)) < 0) ERROR_MSG("H5Tcopy failed");

        if((hdf5_data->status = H5Tset_order(hdf5_data->datatype, H5T_ORDER_LE)) < 0) ERROR_MSG("H5Tset_order failed");

        if((hdf5_data->plist = H5Pcreate (H5P_DATASET_CREATE)) < 0) ERROR_MSG("H5Pcreate failed");
        //status = H5Pset_layout (dcpl, H5D_COMPACT);
        if((hdf5_data->status = H5Pset_layout (hdf5_data->plist, H5D_COMPACT )) < 0 )ERROR_MSG("H5Pset_shuffle failed");

        //if((hdf5_data->status = H5Pset_deflate (hdf5_data->plist, 9)) < 0 )KSLIB_EXCEPTION(kslFAIL,"H5Pset_deflate failed");
        //if((hdf5_data->status = H5Pset_chunk (hdf5_data->plist, hdf5_data->rank,  hdf5_data->chunk_dim)) < 0 )KSLIB_EXCEPTION(kslFAIL,"H5Pset_chunk failed");



        if((hdf5_data->dataset = H5Dcreate(hdf5_data->group, hdf5_data->dataset_name, hdf5_data->datatype, hdf5_data->dataspace,    H5P_DEFAULT, hdf5_data->plist , H5P_DEFAULT)) < 0 )ERROR_MSG("H5Dcreate failed");
        //hdf5_data->plist
        if(hdf5_data->num_attr){
                RUN(hdf5_write_attributes(hdf5_data, hdf5_data->dataset));
        }

        return OK;
ERROR:
        return FAIL;
}



static int hdf5_close_dataset(struct hdf5_data* hdf5_data)
{
        if((hdf5_data->status = H5Sclose(hdf5_data->dataspace)) < 0) ERROR_MSG("H5Sclose failed");
        if((hdf5_data->status = H5Pclose(hdf5_data->plist)) < 0) ERROR_MSG("H5Pclose failed");
        if((hdf5_data->status = H5Tclose(hdf5_data->datatype)) < 0) ERROR_MSG("H5Tclose failed");
        if((hdf5_data->status = H5Dclose(hdf5_data->dataset)) < 0) ERROR_MSG("H5Dclose failed");
        return OK;
ERROR:
        return FAIL;
}



int hdf5_write(char* dataset_name,void* data,struct hdf5_data* hdf5_data)
{


        RUN(hdf5_create_dataset_chunked(dataset_name,hdf5_data, hdf5_data->native_type));

        if((hdf5_data->status  = H5Dwrite(hdf5_data->dataset,hdf5_data->native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)) < 0) ERROR_MSG("H5Dwrite failed");

        hdf5_close_dataset(hdf5_data);
        return OK;
ERROR:
        return FAIL;
}

int hdf5_write_1D_ushortint_expandable(char* dataset_name,unsigned short int* data,int len,struct hdf5_data* hdf5_data)
{
        //open dataset...
        hid_t memspace;
        int i;
        hsize_t old_dim[HDF5GLUE_MAX_DIM];

        hsize_t new_dim[HDF5GLUE_MAX_DIM];

        hsize_t size[HDF5GLUE_MAX_DIM];
        hsize_t offset[HDF5GLUE_MAX_DIM];
        //hsize_t   maxdims[1] = {H5S_UNLIMITED};


        if((hdf5_data->dataset = H5Dopen(hdf5_data->group, dataset_name,H5P_DEFAULT)) == -1)ERROR_MSG("H5Dopen failed\n");
        //dset = H5Dopen (hdf5_data->file, dataset_name, H5P_DEFAULT);

        //int hdf5_open_group(char* groupname,struct hdf5_data* hdf5_data)

        //if((hdf5_data->dataset = H5Dopen(hdf5_data->group, dataset_name,H5P_DEFAULT)) == -1)ERROR_MSG("H5Dopen failed\n");

        hdf5_data->dataspace = H5Dget_space(hdf5_data->dataset);
        hdf5_data->rank     = H5Sget_simple_extent_ndims(hdf5_data->dataspace);
        hdf5_data->status  = H5Sget_simple_extent_dims(hdf5_data->dataspace,old_dim , NULL);

        new_dim[0] = len;


        fprintf(stdout,"rank %d\n", hdf5_data->rank);
        for(i = 0; i < hdf5_data->rank;i++){
                fprintf(stdout,"dim%d = %lu\n",i,(unsigned long)(old_dim[i]));



        }
        if((hdf5_data->status = H5Sclose(hdf5_data->dataspace)) < 0) ERROR_MSG("H5Dclose failed");


        size[0] = old_dim[0] + new_dim[0] ;
        //now bigger..
        hdf5_data->status = H5Dset_extent (hdf5_data->dataset, size);

        //selected
        hdf5_data->dataspace = H5Dget_space(hdf5_data->dataset);
        //select region..
        offset[0] = old_dim[0] ;

        //fprintf(stdout,"%d %d %d \n", old_dim[0], new_dim[0] ,  size[0]  );


        hdf5_data->status = H5Sselect_hyperslab(hdf5_data->dataspace, H5S_SELECT_SET, offset, NULL,     new_dim, NULL);


        memspace = H5Screate_simple (1, new_dim, NULL);





        if((hdf5_data->status  = H5Dwrite(hdf5_data->dataset, H5T_NATIVE_USHORT, memspace, hdf5_data->dataspace , H5P_DEFAULT,  &data[0]))< 0) ERROR_MSG("H5Dwrite failed");

        //status = H5Dclose (dataset);

        //status = H5Sclose (filespace);
        //status = H5Fclose (file);



        //hdf5_data-

        /* Define memory space */
        //memspace = H5Screate_simple (RANK, dimsext, NULL);

        /* Write the data to the extended portion of dataset  */
        //status = H5Dwrite (dataset, H5T_NATIVE_INT, memspace, filespace,			   H5P_DEFAULT, dataext);


        //hdf5_close_dataset(hdf5_data);
        if((hdf5_data->status = H5Dclose(hdf5_data->dataset)) < 0) ERROR_MSG("H5Dclose failed");
        if((hdf5_data->status = H5Sclose(memspace)) < 0) ERROR_MSG("H5Sclose failed");
        if((hdf5_data->status = H5Sclose(hdf5_data->dataspace)) < 0) ERROR_MSG("H5Sclose failed");
        //if((hdf5_data->status = H5Pclose(hdf5_data->plist)) < 0) ERROR_MSG("H5Pclose failed");
        //if((hdf5_data->status = H5Tclose(hdf5_data->datatype)) < 0) ERROR_MSG("H5Tclose failed");



        return OK;
ERROR:
        return FAIL;
        //close dataset
}




//attribute_id = H5Acreate2 (dataset_id, "Units", H5T_STD_I32BE, dataspace_id,   H5P_DEFAULT, H5P_DEFAULT);

/* Write the attribute data. */
//status = H5Awrite(attribute_id, H5T_NATIVE_INT, attr_data);

/* Close the attribute. */
//status = H5Aclose(attribute_id);

#ifdef ITEST

int print_attributes(struct hdf5_data* hdf5_data);

int expand_test(void);

int main (int argc,char * argv[])
{
        int dim1, dim2,dim3,dim4;
        int  i,j,c,f;
        fprintf(stdout,"Running libhdf5glue sanity tests\n");



        char* string = NULL;

        char** m2d_char = NULL;
        int** m2d_int = NULL;
        float** m2d_float = NULL;
        float*** m3d_float = NULL;

        float**** m4d_float = NULL;

        double** m2d_double = NULL;



        MMALLOC(string, sizeof(char) * 10);

        for(i = 0;i < 10;i++){
                string[i] = 0;
        }

        snprintf(string , 10,"%s","DONKEY");





        DPRINTF1("Testing double 2D ");

        dim1 = 3;
        dim2 = 15;
        m2d_char = malloc_2d_char(m2d_char,dim1,dim2,0);
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        m2d_char[i][j] = i+j + 65;
                }
        }

        dim1 = 10;
        dim2 = 10;
        m2d_int = malloc_2d_int(m2d_int,dim1,dim2,0);
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        m2d_int[i][j] = i+j;
                }
        }

        dim1 = 5;
        dim2 = 5;
        m2d_float = malloc_2d_float(m2d_float,dim1,dim2,0.0f);
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        m2d_float[i][j] = i+j + 0.1;
                }
        }


        dim1 = 5;
        dim2 = 5;
        m2d_double = malloc_2d_double(m2d_double,dim1,dim2,0.0);
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        m2d_double[i][j] = i+j + 0.1;
                }
        }

        dim1 = 3;

        dim2 = 5;
        dim3 = 5;

        m3d_float = malloc_3d_float(dim1,dim2,dim3,0.0f);
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        for(c = 0; c < dim3;c++){
                                m3d_float[i][j][c] = (float)i/10.0+  (float)j +  (float)c;
                        }
                }
        }

        dim1 = 2;
        dim2 = 3;
        dim3 = 5;
        dim4 = 5;


        m4d_float = malloc_4d_float(dim1,dim2,dim3,dim4,0.0f);
        for(i =0;i < dim1;i++){
                for(j = 0; j < dim2;j++){
                        for(c = 0; c < dim3;c++){
                                for(f = 0; f < dim4;f++){
                                        m4d_float[i][j][c][f] =(float)i/1000.0 + (float)j/10.0+  (float)c +  (float)f;
                                }
                        }
                }
        }





        struct hdf5_data* hdf5_data = hdf5_create();


        hdf5_add_attribute(hdf5_data, "KALIGN", "READY", 0, 0.0f, HDF5GLUE_CHAR);

        RUN(hdf5_create_file("libhdf5gluetest.h5",hdf5_data));

        hdf5_write_attributes(hdf5_data, hdf5_data->file);

        hdf5_data->num_attr = 0;

        hdf5_add_attribute(hdf5_data, "GROUP1_ATTRIBUTE", "", 4, 0.0f, HDF5GLUE_INT);

        hdf5_create_group("group1",hdf5_data);

        hdf5_write_attributes(hdf5_data, hdf5_data->group);



        hdf5_data->rank = 1;
        hdf5_data->dim[0] = (int)strlen (string) + 1;
        hdf5_data->dim[1] = -1;
        hdf5_data->chunk_dim[0] =  (int)strlen (string) + 1;
        hdf5_data->chunk_dim[1] = -1;

        hdf5_add_attribute(hdf5_data, "DonkeyDATA_ATTRIBUTE1", "FLUP", 4, 0.0f, HDF5GLUE_CHAR);
        hdf5_add_attribute(hdf5_data, "DonkeyDATA_ATTRIBUTE2", "", 4, 0.0f, HDF5GLUE_INT);
        hdf5_add_attribute(hdf5_data, "DonkeyDATA_ATTRIBUTE3", "", 4, 1.0f, HDF5GLUE_FLOAT);
        hdf5_add_attribute(hdf5_data, "DonkeyDATA_ATTRIBUTE4", "", 4, 42.355, HDF5GLUE_FLOAT);
        hdf5_data->native_type = H5T_NATIVE_CHAR;
        hdf5_write("DonkeyDATA",&string[0],hdf5_data );



        DPRINTF1("Write 2D char array");
        hdf5_data->rank = 2;
        hdf5_data->dim[0] = 3;
        hdf5_data->dim[1] = 15;
        hdf5_data->chunk_dim[0] = 3;
        hdf5_data->chunk_dim[1] = 15;

        hdf5_data->native_type = H5T_NATIVE_CHAR;
        hdf5_write("char2D",&m2d_char[0][0], hdf5_data);
        //hdf5_write_2D_char("char2D",m2d_char, hdf5_data);
        DPRINTF1("Write 2D int array");

        hdf5_data->rank = 2;
        hdf5_data->dim[0] = 10;
        hdf5_data->dim[1] = 10;
        hdf5_data->chunk_dim[0] = 10;
        hdf5_data->chunk_dim[1] = 10;
        hdf5_data->native_type = H5T_NATIVE_INT;
        hdf5_write("int2D",&m2d_int[0][0], hdf5_data);


        hdf5_close_group(hdf5_data);

        hdf5_create_group("group2",hdf5_data);


        hdf5_data->rank = 2;
        hdf5_data->dim[0] = 5;
        hdf5_data->dim[1] = 5;
        hdf5_data->chunk_dim[0] = 5;
        hdf5_data->chunk_dim[1] = 5;
        hdf5_data->native_type = H5T_NATIVE_FLOAT;
        hdf5_write("float2d",&m2d_float[0][0], hdf5_data);

        hdf5_data->rank = 2;
        hdf5_data->dim[0] = 5;
        hdf5_data->dim[1] = 5;
        hdf5_data->chunk_dim[0] = 5;
        hdf5_data->chunk_dim[1] = 5;
        hdf5_data->native_type = H5T_NATIVE_DOUBLE ;
        hdf5_write("double2D",&m2d_double[0][0], hdf5_data);
        hdf5_close_group(hdf5_data);



        hdf5_add_attribute(hdf5_data, "THREEDEEE", "INDEED", 4, 0.0f, HDF5GLUE_CHAR);
        hdf5_create_group("group3(D)",hdf5_data);

        hdf5_write_attributes(hdf5_data, hdf5_data->group);

        hdf5_data->rank = 3;
        hdf5_data->dim[0] = 3;
        hdf5_data->dim[1] = 5;
        hdf5_data->dim[2] = 5;
        hdf5_data->chunk_dim[0] = 3;
        hdf5_data->chunk_dim[1] = 5;
        hdf5_data->chunk_dim[2] = 5;
        hdf5_data->native_type = H5T_NATIVE_FLOAT;
        hdf5_write("float3d",&m3d_float[0][0][0], hdf5_data);
        hdf5_close_group(hdf5_data);


        //hdf5_data->num_attr = 0;
        hdf5_add_attribute(hdf5_data, "FOURDEEE", "INDEED", 4, 0.0f, HDF5GLUE_CHAR);
        hdf5_create_group("group4(D)",hdf5_data);
        hdf5_write_attributes(hdf5_data, hdf5_data->group);

        //hdf5_data->num_attr = 0;

        hdf5_data->rank = 4;
        hdf5_data->dim[0] = 2;
        hdf5_data->dim[1] = 3;
        hdf5_data->dim[2] = 5;
        hdf5_data->dim[3] = 5;
        hdf5_data->chunk_dim[0] = 1;
        hdf5_data->chunk_dim[1] = 1;
        hdf5_data->chunk_dim[2] = 1;
        hdf5_data->chunk_dim[3] = 1;
        hdf5_data->native_type = H5T_NATIVE_FLOAT;
        hdf5_write("float4d",&m4d_float[0][0][0][0], hdf5_data);

        hdf5_close_group(hdf5_data);

        hdf5_create_group("group4(D)",hdf5_data);
        hdf5_write_attributes(hdf5_data, hdf5_data->group);

        //hdf5_data->num_attr = 0;

        hdf5_data->rank = 4;
        hdf5_data->dim[0] = 2;
        hdf5_data->dim[1] = 3;
        hdf5_data->dim[2] = 5;
        hdf5_data->dim[3] = 5;
        hdf5_data->chunk_dim[0] = 1;
        hdf5_data->chunk_dim[1] = 1;
        hdf5_data->chunk_dim[2] = 1;
        hdf5_data->chunk_dim[3] = 1;
        hdf5_data->native_type = H5T_NATIVE_FLOAT;

        hdf5_write("float4d_no2",&m4d_float[0][0][0][0], hdf5_data);


        hdf5_close_group(hdf5_data);




        hdf5_close_file(hdf5_data);


        hdf5_free(hdf5_data);


        MFREE(string);

        free_2d((void**) m2d_char);
        free_2d((void**) m2d_int);
        free_2d((void**) m2d_float);
        free_2d((void**) m2d_double);
        free_3d((void***)m3d_float);
        free_4d((void****)m4d_float);

        hdf5_data = NULL;
        hdf5_data = hdf5_create();

        hdf5_open_file("libhdf5gluetest.h5",hdf5_data);
        hdf5_read_attributes(hdf5_data,hdf5_data->file);
        print_attributes(hdf5_data);



        get_group_names(hdf5_data);
        fprintf(stdout,"Groups:\n");
        for(i = 0; i < hdf5_data->grp_names->num_names;i++){
                fprintf(stdout,"%d %s\n",i,hdf5_data->grp_names->names[i]);
        }

        //free_hdf5_group_names(grp_names);


        hdf5_open_group("group1",hdf5_data);
        hdf5_read_attributes(hdf5_data, hdf5_data->group);

        print_attributes(hdf5_data);
        hdf5_read_dataset("DonkeyDATA",hdf5_data);
        print_attributes(hdf5_data);
        if(hdf5_data->data){
                char* string  =hdf5_data->data;
                fprintf(stdout,"%s\n", string);
                MFREE(string);
        }

        hdf5_read_dataset("char2D",hdf5_data);
        print_attributes(hdf5_data);
        if(hdf5_data->data){
                char** m =hdf5_data->data;
                for(i = 0;i < hdf5_data->dim[0];i++){
                        for(j = 0;j < hdf5_data->dim[1];j++){
                                fprintf(stdout," %c",(char) m[i][j]);
                        }
                        fprintf(stdout,"\n");
                }
                fprintf(stdout,"\n");
                free_2d((void**) m);
        }

        hdf5_read_dataset("int2D",hdf5_data);
        print_attributes(hdf5_data);
        if(hdf5_data->data){
                if(hdf5_data->rank == 2){

                        int** int_ptr =hdf5_data->data;
                        for(i = 0;i < hdf5_data->dim[0];i++){
                                for(j = 0;j < hdf5_data->dim[1];j++){
                                        fprintf(stdout," %d", int_ptr[i][j]);
                                }
                                fprintf(stdout,"\n");
                        }
                        fprintf(stdout,"\n");
                        free_2d((void**) int_ptr);
                }
        }


        hdf5_close_group(hdf5_data);


        hdf5_open_group("group2",hdf5_data);

        hdf5_read_attributes(hdf5_data, hdf5_data->group);
        print_attributes(hdf5_data);
        hdf5_read_dataset("float2d",hdf5_data);
        print_attributes(hdf5_data);
        if(hdf5_data->data){
                if(hdf5_data->rank == 2){

                        float** ptr =hdf5_data->data;
                        for(i = 0;i < hdf5_data->dim[0];i++){
                                for(j = 0;j < hdf5_data->dim[1];j++){
                                        fprintf(stdout," %0.1f", ptr[i][j]);
                                }
                                fprintf(stdout,"\n");
                        }
                        fprintf(stdout,"\n");
                        free_2d((void**) ptr);
                }
        }

        hdf5_read_dataset("double2D",hdf5_data);

        print_attributes(hdf5_data);
        if(hdf5_data->data){
                if(hdf5_data->rank == 2){

                        double** ptr =hdf5_data->data;
                        for(i = 0;i < hdf5_data->dim[0];i++){
                                for(j = 0;j < hdf5_data->dim[1];j++){
                                        fprintf(stdout," %0.1f", ptr[i][j]);
                                }
                                fprintf(stdout,"\n");
                        }
                        fprintf(stdout,"\n");
                        free_2d((void**) ptr);
                }
        }

        hdf5_close_group(hdf5_data);

        hdf5_open_group("group3(D)",hdf5_data);
        hdf5_read_attributes(hdf5_data, hdf5_data->group);

        print_attributes(hdf5_data);
        hdf5_read_dataset("float3d",hdf5_data);
        print_attributes(hdf5_data);
        if(hdf5_data->data){
                if(hdf5_data->rank == 3){

                        float*** ptr =hdf5_data->data;
                        for(i = 0;i < hdf5_data->dim[0];i++){
                                for(j = 0;j < hdf5_data->dim[1];j++){
                                        for(c = 0;c < hdf5_data->dim[2];c++){
                                                fprintf(stdout," %0.1f", ptr[i][j][c]);
                                        }
                                        fprintf(stdout,"\n");
                                }
                                fprintf(stdout,"\n");
                                fprintf(stdout,"\n");
                        }
                        fprintf(stdout,"\n");


                        free_3d((void***) ptr);
                }
        }

        hdf5_close_group(hdf5_data);



        hdf5_open_group("group4(D)",hdf5_data);
        hdf5_read_attributes(hdf5_data, hdf5_data->group);

        print_attributes(hdf5_data);
        hdf5_read_dataset("float4d",hdf5_data);
        print_attributes(hdf5_data);
        if(hdf5_data->data){
                if(hdf5_data->rank == 4){

                        float**** ptr =hdf5_data->data;
                        for(i = 0;i < hdf5_data->dim[0];i++){
                                for(j = 0;j < hdf5_data->dim[1];j++){
                                        for(c = 0;c < hdf5_data->dim[2];c++){
                                                for(f = 0;f < hdf5_data->dim[3];f++){
                                                        fprintf(stdout," %0.3f", ptr[i][j][c][f]);
                                                }
                                                fprintf(stdout,"\n");
                                        }
                                        fprintf(stdout,"\n");
                                }
                                fprintf(stdout,"\n");
                                fprintf(stdout,"\n");
                        }
                        fprintf(stdout,"\n");
                        fprintf(stdout,"\n");
                        free_4d((void****) ptr);
                }
        }

        hdf5_close_group(hdf5_data);








        hdf5_close_file(hdf5_data);
        hdf5_free(hdf5_data);


        expand_test();

        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}


int expand_test(void)
{
        struct hdf5_data* hdf5_data = hdf5_create();


        if(!my_file_exists("libhdf5gluetest_exp.h5")){
                RUN(hdf5_create_file("libhdf5gluetest_exp.h5",hdf5_data));
                hdf5_create_group("group2",hdf5_data);

                hid_t native_type =H5T_NATIVE_USHORT;


                hdf5_data->rank = 1;
                hdf5_data->dim[0] = 0;

                hdf5_create_dataset_chunked_expandable("EXPand1",hdf5_data,native_type);
        }else{
                hdf5_open_file("libhdf5gluetest_exp.h5",hdf5_data);
                //hdf5_open_group("group2",hdf5_data);
        }
        hdf5_open_group("group2", hdf5_data);


        unsigned short int* array = NULL;

        MMALLOC(array, sizeof(unsigned short int) * 5);
        int i;

        for(i =0;i < 5;i++){
                array[i] = i;
        }

        hdf5_data->rank = 1;
        hdf5_data->dim[0] = 5;

        hdf5_write_1D_ushortint_expandable("EXPand1",array,5,hdf5_data);


        for(i =0;i < 5;i++){
                array[i] = i+10;
        }

        //hdf5_free(hdf5_data);
        //hdf5_data = NULL;

        //hdf5_data = hdf5_create();

        //hdf5_open_file("libhdf5gluetest_exp.h5",hdf5_data);
        //hdf5_open_group("group2",hdf5_data);

        hdf5_data->rank = 1;
        hdf5_data->dim[0] = 5;

        hdf5_write_1D_ushortint_expandable("EXPand1",array,5,hdf5_data);




        hdf5_read_dataset("EXPand1",hdf5_data);
        print_attributes(hdf5_data);
        if(hdf5_data->data){
                if(hdf5_data->rank == 1){

                        unsigned short int* ptr =hdf5_data->data;
                        for(i = 0;i < hdf5_data->dim[0];i++){
                                //for(j = 0;j < hdf5_data->dim[1];j++){
                                fprintf(stdout," %d", ptr[i]);
                                //}
                                //fprintf(stdout,"\n");
                        }
                        fprintf(stdout,"\n");
                        free((void*) ptr);
                }
        }



        hdf5_close_file(hdf5_data);
        hdf5_free(hdf5_data);

        MFREE(array);
        return OK;
ERROR:
        return FAIL;
}





#endif


int print_attributes(struct hdf5_data* hdf5_data)
{
        int i;
        if(hdf5_data->num_attr){
                for(i = 0; i < hdf5_data->num_attr;i++){
                        switch (hdf5_data->attr[i]->type) {
                        case HDF5GLUE_CHAR:
                                fprintf(stdout,"%s => %s\n", hdf5_data->attr[i]->attr_name, hdf5_data->attr[i]->string);
                                break;
                        case HDF5GLUE_INT:
                                fprintf(stdout,"%s => %d\n", hdf5_data->attr[i]->attr_name, hdf5_data->attr[i]->int_val);
                                break;
                        case HDF5GLUE_FLOAT:
                                fprintf(stdout,"%s => %f\n", hdf5_data->attr[i]->attr_name, hdf5_data->attr[i]->float_val);
                                break;
                        default:
                                break;
                        }
                }
        }
        return OK;
}
