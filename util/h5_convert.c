#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>

int main( int argc, char *argv[] ){

     int NX, NY, NZ;
     FILE *fp_in;

     char filename[50], hdf5_filename[50], dsetname[10];
     int i, j, k;
     float *a;

     hid_t fpout, dset, dspace;
     hsize_t dims[3];
     herr_t status;

 /**
  ** Read in the binary file's name from the command line 
  **-------------------------------------------------------------------*/
     if ( (argc == 1)||(argc>3) ) {
        printf( "usage:\n\n" );
        printf( "      h5convert <binary filename> <dataset name>\n\n" );
        return 0;
     }
     if ( argc>1 ) {
        sprintf( filename, "%s", argv[1] );
     }
     if ( argc==3 ) {
        sprintf( dsetname, "%s", argv[2] );
     } else {
        sprintf( dsetname, "data" );
     }
     printf( "\nInput binary file: %s\n", filename );

 /**
  ** Set the dataset dimensions 
  **-------------------------------------------------------------------*/
     NX = 160; 
     NY = 184;
     NZ = 160;

 /**
  ** Open the original binary file and read its contents 
  **-------------------------------------------------------------------*/
     a = (float *) malloc( NY*NX*NZ*sizeof(float) );

     fp_in = fopen( filename, "r" );
     if ( fp_in==NULL ) {
        printf( "\nERROR: could not find a binary datafile called %s\n\n", filename );
        return 1;
     }
     fread( a, sizeof(float), NY*NX*NZ, fp_in );	
     fclose( fp_in );

 /**
  ** Create a new HDF5 file and dataset 
  **-------------------------------------------------*/
     sprintf( hdf5_filename, "%s.h5", filename ); 
     printf( "Output HDF5 file : %s\n\n", hdf5_filename );
     fpout = H5Fcreate( hdf5_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

 /**
  ** Define the 3D dataset to be written into the new HDF5 file. Remember that:
  **   Y -> vertical direction into the ground
  **   X -> horizontal direction along the ground (parallel to your p.o.v) 
  **   Z -> horizontal direction along the ground (away from you)
  **---------------------------------------------------------------------------*/
     dims[0] = NY;
     dims[1] = NX;
     dims[2] = NZ;
     dspace = H5Screate_simple( 3, dims, NULL );

     dset = H5Dcreate( fpout, dsetname, H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

 /**
  ** Writes data into the new HDF5 file 
  **------------------------------------------------------------------------*/
     status = H5Dwrite( dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, a );

     status = H5Dclose( dset );
     status = H5Sclose( dspace );
     status = H5Fclose( fpout );

     free( a );
     a = NULL;
}

