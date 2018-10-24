/***
 ***  FILE COMPARE
 ***
 *** Simple program that compares the contents of a binary snapshot file to its
 *** corresponding HDF5 counterpart.
 ***
 *** Author: Mark Cheeseman, Pawsey Supercomputing Center
 *** Date  : October 19, 2018
 ***
 ***===========================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>

int main( int argc, char *argv[] ){

     int NX, NY, NZ;
     FILE *fp_in;

     char binary_file[50], hdf5_file[50], dsetname[10];
     int n, count;
     float *binary_data, *hdf5_data, comparison_ratio;

     hid_t fpout, dset;
     herr_t status;

 /**
  ** Set the dataset dimensions 
  **-------------------------------------------------------------------*/
     NX = 160; 
     NY = 184;
     NZ = 160;

 /**
  ** Parse command line arguments
  **-------------------------------------------------------------------*/
     if ( argc != 3 ) {
        printf( "usage:\n\n" );
        printf( "      h5dump <binary filename> <dataset name>\n\n" );
        return 0;
     }
     sprintf( binary_file, "%s", argv[1] );
     sprintf( dsetname, "%s", argv[2] );

 /**
  ** Check if the binary file exists.  If so, open it and read its
  ** contents. 
  **-------------------------------------------------------------------*/
     fp_in = fopen( binary_file, "r" );
     if ( fp_in==NULL ) {
        printf( "\nERROR: could not find a binary datafile called %s\n\n", binary_file );
        return 1;
     }

     binary_data = (float *) malloc( NY*NX*NZ*sizeof(float) );
     fread( binary_data, sizeof(float), NY*NX*NZ, fp_in );
     fclose( fp_in );

 /**
  ** Check if the corresponding HDF5 file exists.  If so, open it and 
  ** read its contents. 
  **-------------------------------------------------------------------*/
     sprintf( hdf5_file, "%s_new.h5", binary_file ); 
     fpout = H5Fopen( hdf5_file, H5F_ACC_RDONLY, H5P_DEFAULT );
     if ( fpout<0 ) {
        printf( "\nERROR: could not find the hdf5 file called %s\n\n", hdf5_file );
        return 1;
     }

     dset = H5Dopen2( fpout, dsetname, H5P_DEFAULT );
     if ( dset<0 ) {
        printf( "\nERROR: could not find dsetname %s in hdf5 file\n\n", dsetname );
        return 1;
     }

     hdf5_data = (float *) malloc( NY*NX*NZ*sizeof(float) );
     status = H5Dread( dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                       hdf5_data );

     status = H5Dclose( dset );
     status = H5Fclose( fpout );

 /**
  ** iDetermine the amount of misfit between the contents of the binary
  ** and HDF5 data files. 
  **-------------------------------------------------------------------*/
     count = NX*NY*NZ;
     for ( n=0; n<NX*NY*NZ; n++ ) {
         if ( abs( binary_data[n]-hdf5_data[n])>0.0000001 ) {
            count--;
         }
     }

     comparison_ratio = ((float )(count) / (float ) (NX*NY*NZ)) * 100.0;
     printf( "%4.1f percent fit determined\n\n", comparison_ratio );
  
 /**
  ** Clean-up 
  **-------------------------------------------------------------------*/
     free( binary_data );
     free( hdf5_data );
     binary_data = NULL;
     hdf5_data = NULL;

     return 0;
}

