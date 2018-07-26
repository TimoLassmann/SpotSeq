
#ifndef hdf5_glue_header

#define hdf5_glue_header

#define HDF5GLUE_MAX_NAME_LEN 200
#define HDF5GLUE_MAX_CONTENT_LEN 1000
#define HDF5GLUE_MAX_DIM 5

#define HDF5GLUE_CHAR 1
#define HDF5GLUE_INT 2
#define HDF5GLUE_USHORT_INT 3

#define HDF5GLUE_FLOAT 4
#define HDF5GLUE_DOUBLE 5

#include "hdf5.h"


struct hdf5_attribute{
        char attr_name[HDF5GLUE_MAX_NAME_LEN];
        char string[HDF5GLUE_MAX_CONTENT_LEN];
        int int_val;
        float float_val;
        int type;
};

struct hdf5_group_names{
        char** names;
        int num_names;
        int name_length;
        int num_names_mem;
};

struct hdf5_data{
        char dataset_name[HDF5GLUE_MAX_NAME_LEN];
        char group_name[HDF5GLUE_MAX_NAME_LEN];
        char file_name[HDF5GLUE_MAX_NAME_LEN];
        char tmp_name[HDF5GLUE_MAX_NAME_LEN];
	
        hsize_t dim[HDF5GLUE_MAX_DIM];
        hsize_t chunk_dim[HDF5GLUE_MAX_DIM];
        struct hdf5_group_names* grp_names;
        struct hdf5_attribute** attr;
        void* data;
        int num_attr;
        int num_attr_mem;
        int rank;
	
        hid_t fapl;
        hid_t	 file;
        hid_t group;
	
        hid_t plist;
        hid_t dataset;
	
        hid_t attribute_id;
        hid_t attr_dataspace_id;
	
        hid_t datatype;
        hid_t dataspace;
	
        herr_t status;
};



extern struct hdf5_data* hdf5_create(void);
extern int hdf5_free(struct hdf5_data* hdf5_data);

extern int hdf5_create_file(char* filename,struct hdf5_data* hdf5_data);
extern int hdf5_create_group(char* groupname,struct hdf5_data* hdf5_data);

extern int hdf5_close_group(struct hdf5_data* hdf5_data);
extern int hdf5_close_file(struct hdf5_data* hdf5_data);

extern int hdf5_add_attribute(struct hdf5_data* hdf5_data,char* attr_name, char* string, int int_val, float float_val, int type);

extern int hdf5_write_attributes(struct hdf5_data* hdf5_data,hid_t target);
extern int hdf5_read_attributes(struct hdf5_data* hdf5_data,hid_t target);


extern int print_attributes(struct hdf5_data* hdf5_data);


extern int hdf5_write_1D_char(char* dataset_name, char* data,struct hdf5_data* hdf5_data);
extern int hdf5_write_1D_int(char* dataset_name, int* data,struct hdf5_data* hdf5_data);
extern int hdf5_write_1D_ushortint(char* dataset_name,unsigned short int* data,struct hdf5_data* hdf5_data);
extern int hdf5_write_1D_float(char* dataset_name, float* data,struct hdf5_data* hdf5_data);
extern int hdf5_write_1D_double(char* dataset_name, double* data,struct hdf5_data* hdf5_data);

extern int hdf5_write_2D_char(char* dataset_name, char** data,struct hdf5_data* hdf5_data);
extern int hdf5_write_2D_int(char* dataset_name,int** data,struct hdf5_data* hdf5_data);
extern int hdf5_write_2D_float(char* dataset_name,float** data,struct hdf5_data* hdf5_data);
extern int hdf5_write_2D_double(char* dataset_name, double** data,struct hdf5_data* hdf5_data);

extern int hdf5_write_3D_float(char* dataset_name,float*** data,struct hdf5_data* hdf5_data);

extern int hdf5_write_4D_float(char* dataset_name,float**** data,struct hdf5_data* hdf5_data);

extern int hdf5_open_file(char* filename,struct hdf5_data* hdf5_data);
extern int hdf5_open_group(char* groupname,struct hdf5_data* hdf5_data);
extern int hdf5_read_dataset(char* dataset_name,struct hdf5_data* hdf5_data);

extern int get_group_names(struct hdf5_data* hdf5_data);

extern int hdf5_create_dataset_chunked_expandable(char* dataset_name,struct hdf5_data* hdf5_data, hid_t type);
extern int hdf5_create_dataset_chunked(char* dataset_name,struct hdf5_data* hdf5_data, hid_t type);
extern int hdf5_write_1D_ushortint_expandable(char* dataset_name,unsigned short int* data,int len, struct hdf5_data* hdf5_data);

#ifndef STRINGIFY
#define STRINGIFY(x) #x

#endif

#ifndef TOSTRING
#define TOSTRING(x) STRINGIFY(x)
#endif


#ifndef AT
#define AT __FILE__ " line " TOSTRING(__LINE__)
#endif

#define RUN_HDF5(EXP) do {                                        \
                if((EXP) < 0){                                    \
                        ERROR_MSG_HDF5(TOSTRING(EXP),"failed.");  \
                }                                                 \
        }while (0)



#define ERROR_MSG_HDF5(...) do {                \
                error(AT, __VA_ARGS__ );        \
                H5Eprint2(H5E_DEFAULT, stdout); \
                goto ERROR;                     \
        }while (0)

#endif
