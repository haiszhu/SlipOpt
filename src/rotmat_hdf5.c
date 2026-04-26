#include <stdio.h>
#include <stdint.h>
#include <string.h>
#ifdef HAVE_HDF5_CACHE
#include <hdf5.h>
#endif

#ifndef SLIPOPT_ROOT
#define SLIPOPT_ROOT "."
#endif

static void build_path(int p, char *buf, size_t nbuf) {
  (void)snprintf(buf, nbuf, "%s/data/rotmat_p%03d.h5", SLIPOPT_ROOT, p);
}

static void build_alt_path(int p, char *buf, size_t nbuf) {
  (void)snprintf(buf, nbuf, "data/rotmat_p%03d.h5", p);
}

#ifndef HAVE_HDF5_CACHE
int rotmat_h5_read_p(int p, int np, int nu, double *rall) {
  (void)p; (void)np; (void)nu; (void)rall;
  return 0;
}

int rotmat_h5_write_p(int p, int np, int nu, const double *rall) {
  (void)p; (void)np; (void)nu; (void)rall;
  return 0;
}
#else
int rotmat_h5_read_p(int p, int np, int nu, double *rall) {
  char path[256], path2[256];
  hid_t fid = -1, dset = -1, space = -1;
  hid_t err_stack = H5E_DEFAULT;
  H5E_auto2_t old_func = NULL;
  void *old_client_data = NULL;
  hsize_t dims[3] = {0, 0, 0};
  int ndims;
  herr_t st;

  build_path(p, path, sizeof(path));
  build_alt_path(p, path2, sizeof(path2));

  /* Missing cache file is expected on first call; silence HDF5 stderr noise. */
  H5Eget_auto2(err_stack, &old_func, &old_client_data);
  H5Eset_auto2(err_stack, NULL, NULL);
  fid = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (fid < 0) fid = H5Fopen(path2, H5F_ACC_RDONLY, H5P_DEFAULT);
  H5Eset_auto2(err_stack, old_func, old_client_data);
  if (fid < 0) return 0;

  dset = H5Dopen2(fid, "/rall", H5P_DEFAULT);
  if (dset < 0) goto fail;

  space = H5Dget_space(dset);
  if (space < 0) goto fail;

  ndims = H5Sget_simple_extent_ndims(space);
  if (ndims != 3) goto fail;

  st = H5Sget_simple_extent_dims(space, dims, NULL);
  if (st < 0) goto fail;
  if ((int)dims[0] != np || (int)dims[1] != np || (int)dims[2] != nu) goto fail;

  st = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rall);
  if (st < 0) goto fail;

  H5Sclose(space);
  H5Dclose(dset);
  H5Fclose(fid);
  return 1;

fail:
  if (space >= 0) H5Sclose(space);
  if (dset >= 0) H5Dclose(dset);
  if (fid >= 0) H5Fclose(fid);
  return 0;
}

int rotmat_h5_write_p(int p, int np, int nu, const double *rall) {
  char path[256], path2[256];
  hid_t fid = -1, dset = -1, space = -1;
  hsize_t dims[3];
  herr_t st;

  build_path(p, path, sizeof(path));
  build_alt_path(p, path2, sizeof(path2));

  dims[0] = (hsize_t)np;
  dims[1] = (hsize_t)np;
  dims[2] = (hsize_t)nu;

  fid = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (fid < 0) fid = H5Fcreate(path2, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (fid < 0) return 0;

  space = H5Screate_simple(3, dims, NULL);
  if (space < 0) goto fail;

  dset = H5Dcreate2(fid, "/rall", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dset < 0) goto fail;

  st = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rall);
  if (st < 0) goto fail;

  H5Dclose(dset);
  H5Sclose(space);
  H5Fclose(fid);
  return 1;

fail:
  if (dset >= 0) H5Dclose(dset);
  if (space >= 0) H5Sclose(space);
  if (fid >= 0) H5Fclose(fid);
  return 0;
}
#endif
