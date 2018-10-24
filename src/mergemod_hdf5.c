/*------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS3D.
 * 
 * IFOS3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * IFOS3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IFOS3D. See file COPYING and/or 
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
--------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   merge model files written by the different processes to 
 *   a single HDF5 file                                 
 *   last update 15/10/18   M. Cheeseman (Pawsey)
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include <hdf5.h>

void mergemod_hdf5(char modfile[STRING_SIZE], char dsetname[10], float ***val){

     extern int POS[4], NXG, NYG, NZG, NX, NY, NZ;

     char file[STRING_SIZE], filename[STRING_SIZE];
     int i,j,k,ind;
     float *buf;

     hid_t fpout, dset, dspace, mspace, plist;
     hsize_t dims[3], offset[3], count[3], block[3], stride[3];
     herr_t status;

 /**
  ** Strip off the halo rows and columns from the input data array.
  **---------------------------------------------------------------*/
     buf = (float *) malloc( sizeof(float)*NX*NY*NZ );

     ind = 0;
     for ( k=0; k<NZ; k++ ) {
     for ( i=0; i<NX; i++ ) {
     for ( j=0; j<NY; j++ ) {
         buf[ind] = val[j+1][i+1][k+1];
         ind++;
     }}}

 /**
  ** Create a new HDF5 and enable parallel I/O on it.
  **-------------------------------------------------*/
     sprintf( filename, "%s_new.h5", modfile ); 

     plist = H5Pcreate( H5P_FILE_ACCESS );
     H5Pset_fapl_mpio( plist, MPI_COMM_WORLD, MPI_INFO_NULL );
     fpout = H5Fcreate( filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist );
     status = H5Pclose( plist );

 /**
  ** Define the 3D dataset to be written into the new HDF5 file. Remember that:
  **   Y -> vertical direction into the ground
  **   X -> horizontal direction along the ground (parallel to your p.o.v) 
  **   Z -> horizontal direction along the ground (away from you)
  **---------------------------------------------------------------------------*/
     dims[0] = NYG;
     dims[1] = NXG;
     dims[2] = NZG;
     dspace = H5Screate_simple( 3, dims, NULL );

     block[0] = NY;
     block[1] = NX;
     block[2] = NZ;
     mspace = H5Screate_simple( 3, block, NULL );

     plist = H5Pcreate( H5P_DATASET_CREATE );
     status = H5Pset_chunk( plist, 3, block );
     dset = H5Dcreate( fpout, dsetname, H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, plist, H5P_DEFAULT );
     status = H5Pclose( plist );
     status = H5Sclose( dspace );

 /**
  ** Each MPI task then writes its data into the new HDF5 file concurrently
  **------------------------------------------------------------------------*/
     for ( k=0; k<3; k++ ) {
         stride[k] = 1;
         count[k] = 1;
     }

     offset[0] = POS[2]*NY;
     offset[1] = POS[1]*NX;
     offset[2] = POS[3]*NZ;
     dspace = H5Dget_space( dset );
     status = H5Sselect_hyperslab( dspace, H5S_SELECT_SET, offset, stride, count, block );
       
     plist = H5Pcreate( H5P_DATASET_XFER );
     status = H5Pset_dxpl_mpio( plist, H5FD_MPIO_COLLECTIVE );     
     status = H5Dwrite( dset, H5T_NATIVE_FLOAT, mspace, dspace, plist, buf );

 /**
  ** Clean up
  **----------*/
     status = H5Dclose( dset );
     status = H5Sclose( dspace );
     status = H5Sclose( mspace );
     status = H5Pclose( plist );
     status = H5Fclose( fpout );

     free( buf );
     buf = NULL;
}

