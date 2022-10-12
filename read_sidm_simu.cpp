#include <iostream>
#include <string>
#include <H5Cpp.h>
#include <H5Apublic.h>
#include <stdlib.h>

#include "read_sidm_simu.h"

using namespace std;

float read_float_attr(string path, const char *attr_name){

    float result = -1;

    int status = -1;

    hid_t file_id, attr;

    file_id = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    attr = H5Aopen_by_name(file_id, ".", attr_name, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Aread(attr, H5T_NATIVE_FLOAT, &result);

    status = H5Aclose(attr);

    status = H5Fclose(file_id);

    return result;
}

int read_int_attr(string path, const char *attr_name){

    int result = -1;

    int status = -1;

    hid_t file_id, attr;

    file_id = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    attr = H5Aopen_by_name(file_id, ".", attr_name, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Aread(attr, H5T_NATIVE_INT, &result);

    status = H5Aclose(attr);

    status = H5Fclose(file_id);

    return result;
}

float *read_dataset(string path, string dset_name){
    hid_t file_id, dataset_type, dataset_id, dspace_id;
    hid_t mem_type, group_id;
    hsize_t nobjects;
    float *result = NULL;

    file_id = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset_id = H5Dopen(file_id, dset_name.c_str(), H5P_DEFAULT);
    dataset_type = H5Dget_type(dataset_id);
    
    cout << "dtype: " << dataset_type << endl;    
    cout << "H5T_NATIVE_FLOAT: " << H5T_NATIVE_FLOAT << endl;    
    
    if (H5Tequal(H5T_NATIVE_FLOAT, dataset_type)){

        dspace_id = H5Dget_space(dataset_id);
        H5Sget_simple_extent_dims(dspace_id, &nobjects, NULL);
        H5Sclose(dspace_id);

        result = (float *) malloc(nobjects * sizeof(H5T_NATIVE_FLOAT));

        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, result);

    }

    H5Tclose(dataset_type);

    H5Dclose(dataset_id);

    H5Fclose(file_id);

    return result;
}
