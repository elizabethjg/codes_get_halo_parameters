#ifndef READ_SIDM_SIMU_H
#define READ_SIDM_SIMU_H

using namespace std;

float read_float_attr(string path, const char *attr_name);

int read_int_attr(string path, const char *attr_name);

float *read_dataset(string path, string dset_name);

#endif
