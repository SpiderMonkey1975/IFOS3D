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
 *   Write seismograms to disk
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "segy.h"
#include <hdf5.h>

/* ****************************  UNDER CONSTRUCTION !!!  ***************************** */

void outseis_hdf5( hid_t gid, char dsetname[10], float **data, int **recpos_local, int ntr_glob, 
                   int ntr_local, int ns ) {

     extern int MYID;
     int n, ierr;

     hid_t dspace, dset, plist, mspace;
     herr_t status;
     hsize_t dims[2], block[2], stride[2], count[2], offset[2];
	
 /**
  ** Define a 2D dataset to be written into the new HDF5 file with:
  **   dim[0] -> ID of the receiver data (1 to ntr_global)
  **   dim[1] -> # of samples (ns) 
  **---------------------------------------------------------------------------*/
     dims[0] = ntr_glob;
     dims[1] = ns;
     dspace = H5Screate_simple( 2, dims, NULL );

     block[0] = 1;
     block[1] = ns;
     mspace = H5Screate_simple( 2, block, NULL );

     plist = H5Pcreate( H5P_DATASET_CREATE );
     status = H5Pset_chunk( plist, 2, block );
     dset = H5Dcreate( gid, dsetname, H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, plist, H5P_DEFAULT );
     if ( dset<0 ) {
        n = 0;
        if (MYID==0) { printf( "\n\nERROR: could not create HDF5 dataset %s\n\n", dsetname ); }
        ierr = MPI_Abort( MPI_COMM_WORLD, n );
     }
     status = H5Pclose( plist );

 /**
  ** Each MPI task then writes its ireceiver data at locations it possesses
  ** independently.
  **------------------------------------------------------------------------*/
     for ( n=0; n<2; n++ ) {
         stride[n] = 1;
         count[n] = 1;
     }
     offset[1] = 0;

     plist = H5Pcreate( H5P_DATASET_XFER );
     status = H5Pset_dxpl_mpio( plist, H5FD_MPIO_INDEPENDENT );

     for ( n=1; n<=ntr_local; n++ ) {

         offset[0] = recpos_local[4][n];
         if ( offset[0]==ntr_glob ) { offset[0] = ntr_glob-1; }

         status = H5Sselect_hyperslab( dspace, H5S_SELECT_SET, offset, stride, count, block );
         status = H5Dwrite( dset, H5T_NATIVE_FLOAT, mspace, dspace, plist, &data[n][1] );
         if ( status<0 ) {
            printf( "\nMPI Task %d has offset of [%d,%d]\n", MYID, offset[0], offset[1] );
         }

     }

     status = H5Dclose( dset );
     status = H5Sclose( mspace );
     status = H5Sclose( dspace );
     status = H5Pclose( plist );
     
//	for (tracl=1; tracl<=ntr; tracl++) {
//	   tr.data[j]=section[tracl][1];
//	   fwrite(&tr.data[1],4,ns,fpdata);
//	}

 }

