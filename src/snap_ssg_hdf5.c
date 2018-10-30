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
 *   Write 3D snapshot for current timestep to disk                                   
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include <unistd.h>
#include <hdf5.h>

void snap_hdf5( FILE *fp, int nt, int nsnap, int type, float ***vx, float ***vy, float ***vz, 
                float ***sxx, float ***syy, float ***szz, float ***u, float ***pi, int idx, int idy, int idz, 
                int nx1, int ny1, int nz1, int nx2, int ny2, int nz2 ) {

	/* 
	different types:
	type=1 : values in vx, vy, and vz
	type=2 : -(sxx+syy+szz) (pressure field)
	type=3 : divergence of vx, vy and vz (energy of compressional waves)
	         and curl of vx, vy and vz (energy of shear waves)
	type=4 : both particle velocities (type=1) and energy (type=3)
	*/
	
        hid_t fid, plist, dspace, mspace, dset[3];
        hsize_t block[3], offset[3], stride[3], dims[3], count[3];
        herr_t status;

	char filename[STRING_SIZE], ext[8], wm[1];
	int i,j,k;
	float tmp, a=0.0, ***amp, dh24x, dh24y, dh24z, vyx, vxy, vxx, vyy, vzx, vyz, vxz, vzy, vzz;

	extern float DX, DY, DZ, DT;
	extern char SNAP_FILE[STRING_SIZE];
	extern int MYID, SNAP_PLANE, LOG, POS[4], NX, NY, NZ, NXG, NYG, NZG;

     /**
      ** Create or open the appropriate HDF5 file 
      **/

        sprintf( filename, "%s.h5" );

        plist = H5Pcreate( H5P_FILE_ACCESS );
        H5Pset_fapl_mpio( plist, MPI_COMM_WORLD, MPI_INFO_NULL );

        if ( access( filename, F_OK ) != -1 ) {
           fid = H5Fopen( filename, H5F_ACC_RDWR, plist );
        } else {
           fid = H5Fcreate( filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist );
        }

        if ( fid<0 ) {
           if (MYID==0) { printf( "\n\nERROR: could not create the snapfile: %s\n\n", filename ); }
           k = 0;
           i = MPI_Abort( MPI_COMM_WORLD, k );
        }
        status = H5Pclose( plist );

     /**
      ** Create the appropriate dataset 
      **/

        dims[0] = NYG;
        dims[1] = NXG;
        dims[2] = NZG;
        dspace = H5Screate_simple( 3, dims, NULL );

        block[0] = NY;
        block[1] = NX;
        block[2] = NZ;
        mspace = H5Screate_simple( 3, block, NULL );

        for ( k=0; k<3; k++ ) {
            stride[k] = 1;
            count[k] = 1;
        }

        offset[0] = POS[2]*NY;
        offset[1] = POS[1]*NX;
        offset[2] = POS[3]*NZ;

        plist = H5Pcreate( H5P_DATASET_CREATE );
        status = H5Pset_chunk( plist, 3, block );

	switch(type){
	  case 1 : // particle velocities
               dset[0] = H5Dcreate( fid, "vx", H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, plist, H5P_DEFAULT );
               dset[1] = H5Dcreate( fid, "vy", H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, plist, H5P_DEFAULT );
               dset[2] = H5Dcreate( fid, "vz", H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, plist, H5P_DEFAULT );
               status = H5Pclose( plist );
               status = H5Sclose( dspace );

               plist = H5Pcreate( H5P_DATASET_XFER );
               status = H5Pset_dxpl_mpio( plist, H5FD_MPIO_COLLECTIVE );

               dspace = H5Dget_space( dset[0] );
               status = H5Sselect_hyperslab( dspace, H5S_SELECT_SET, offset, stride, count, block );
               status = H5Dwrite( dset[0], H5T_NATIVE_FLOAT, mspace, dspace, plist, vx );
               status = H5Sclose( dspace );
               status = H5Dclose( dset[0] );

               dspace = H5Dget_space( dset[1] );
               status = H5Sselect_hyperslab( dspace, H5S_SELECT_SET, offset, stride, count, block );
               status = H5Dwrite( dset[1], H5T_NATIVE_FLOAT, mspace, dspace, plist, vy );
               status = H5Sclose( dspace );
               status = H5Dclose( dset[1] );

               dspace = H5Dget_space( dset[2] );
               status = H5Sselect_hyperslab( dspace, H5S_SELECT_SET, offset, stride, count, block );
               status = H5Dwrite( dset[2], H5T_NATIVE_FLOAT, mspace, dspace, plist, vz );
               status = H5Sclose( dspace );
               status = H5Dclose( dset[2] );
	       break;

	  case 2 : // pressure
               dset[0] = H5Dcreate( fid, "pressure", H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, plist, H5P_DEFAULT );
               status = H5Pclose( plist );
               status = H5Sclose( dspace );

               amp = f3tensor(0,NY+1,0,NX+1,0,NZ+1);
	       for (k=ny1;k<=ny2;k+=idy){
	       for (i=nx1;i<=nx2;i+=idx){
	       for (j=nz1;j<=nz2;j+=idz){
		    amp[j][i][k] = -sxx[j][i][k]-syy[j][i][k]-szz[j][i][k];
	       }}}
	
               plist = H5Pcreate( H5P_DATASET_XFER );
               status = H5Pset_dxpl_mpio( plist, H5FD_MPIO_COLLECTIVE );
               dspace = H5Dget_space( dset[0] );
               status = H5Sselect_hyperslab( dspace, H5S_SELECT_SET, offset, stride, count, block );
               status = H5Dwrite( dset[0], H5T_NATIVE_FLOAT, mspace, dspace, plist, amp );
               status = H5Sclose( dspace );               
               status = H5Dclose( dset[0] );               
               
               free_f3tensor(amp,0,NY+1,0,NX+1,0,NZ+1);
	       break;

	  case 4 :
               if ( MYID==0 ) {
                  printf( "\n\ntype option 4 no longer supported with HDF5 enabled.  Make 2 seperate calls instead: option 1 then option 2\n\n" );
                  i = MPI_Abort( MPI_COMM_WORLD, k );
               }
               break;

	  case 3 : // curl of the velocity field according to Dougherty and Stephen (PAGEOPH, 1988) 
               dset[0] = H5Dcreate( fid, "divergence", H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, plist, H5P_DEFAULT );
               dset[1] = H5Dcreate( fid, "curl", H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, plist, H5P_DEFAULT );
               status = H5Pclose( plist );
               status = H5Sclose( dspace );
		
               plist = H5Pcreate( H5P_DATASET_XFER );
               status = H5Pset_dxpl_mpio( plist, H5FD_MPIO_COLLECTIVE );

	       dh24x = 1.0/DX;
	       dh24y = 1.0/DY;
	       dh24z = 1.0/DZ;
		
               amp = f3tensor(0,NY+1,0,NX+1,0,NZ+1);
	       for (k=ny1;k<=ny2;k+=idy){
	       for (i=nx1;i<=nx2;i+=idx){
	       for (j=nz1;j<=nz2;j+=idz){
 		    vxy = (vx[j+1][i][k]-vx[j][i][k])*(dh24y);
		    vxz = (vx[j][i][k+1]-vx[j][i][k])*(dh24z);
		    vyx = (vy[j][i+1][k]-vy[j][i][k])*(dh24x);
		    vyz = (vy[j][i][k+1]-vy[j][i][k])*(dh24z);
		    vzx = (vz[j][i+1][k]-vz[j][i][k])*(dh24x);
		    vzy = (vz[j+1][i][k]-vz[j][i][k])*(dh24y);
					
		/*amp= absolute value of curl(v)), without sqrt!!!*/
		    a = ((vzy-vyz)*(vzy-vyz)+(vxz-vzx)*(vxz-vzx)+(vyx-vxy)*(vyx-vxy)); 
					
		/*note that "y" is used for the vertical coordinate*/
		    switch( SNAP_PLANE ){
			case 1 : /* energy without sign */
				amp[j][i][k]=sqrt((u[j][i][k])*a);
				break;
			case 2 : /* energy with sign true for x-y-plane */
                                tmp = (float ) fsign((vxy-vyx)); 
				amp[j][i][k] = tmp*sqrt((u[j][i][k])*a);
				break;
			case 3 : /* energy with sign true for x-z-plane */
                                tmp = (float ) fsign((vxz-vzx));
				amp[j][i][k] = tmp*sqrt((u[j][i][k])*a);
				break;
			case 4 :/* energy with sign true for y-z-plane */
                                tmp = (float ) fsign((vzy-vyz));
				amp[j][i][k] = tmp*sqrt((u[j][i][k])*a);
				break;
			case 5 : /*custom force*/ /*not yet working properly*/
                                tmp = (float ) fsign((vxz-vzx));
				amp[j][i][k] = tmp*sqrt((u[j][i][k])*a);
				break;
		    }
	       }}}

               dspace = H5Dget_space( dset[0] );
               status = H5Sselect_hyperslab( dspace, H5S_SELECT_SET, offset, stride, count, block );
               status = H5Dwrite( dset[0], H5T_NATIVE_FLOAT, mspace, dspace, plist, amp );
               status = H5Sclose( dspace );               
               status = H5Dclose( dset[0] );               

		/* output of the divergence of the velocity field according to Dougherty and
		                  Stephen (PAGEOPH, 1988) */
	       for (k=ny1;k<=ny2;k+=idy){
	       for (i=nx1;i<=nx2;i+=idx){
	       for (j=nz1;j<=nz2;j+=idz){
					
		    vxx = (vx[j][i][k]-vx[j][i-1][k])*(dh24x);
		    vyy = (vy[j][i][k]-vy[j-1][i][k])*(dh24y);
		    vzz = (vz[j][i][k]-vz[j][i][k-1])*(dh24z);
					
		    a = (vxx+vyy+vzz);
					
		    switch(SNAP_PLANE){
			case 1 : /* energy without sign */
				/* Ep with Ep=pi*amp*amp */
				amp[j][i][k] = sqrt((pi[j][i][k])*a*a);
				break;
			case 2 : /* single force in x */
				/*sign of div(v) * Ep with Ep=pi*amp*amp */
                                tmp = (float ) fsign(a);
				amp[j][i][k] = tmp*sqrt((pi[j][i][k])*a*a);
				break;
			case 3 : /* single force in y */
				/*sign of div(v) * Ep with Ep=pi*amp*amp */
                                tmp = (float ) fsign(a);
				amp[j][i][k] = tmp*sqrt((pi[j][i][k])*a*a);
				break;
			case 4 : /* single force in z */
				/*sign of div(v) * Ep with Ep=pi*amp*amp */
                                tmp = (float ) fsign(a);
				amp[j][i][k] = tmp*sqrt((pi[j][i][k])*a*a);
				break;
			}
	       }}}

               dspace = H5Dget_space( dset[1] );
               status = H5Sselect_hyperslab( dspace, H5S_SELECT_SET, offset, stride, count, block );
               status = H5Dwrite( dset[1], H5T_NATIVE_FLOAT, mspace, dspace, plist, amp );
               status = H5Sclose( dspace );
               status = H5Dclose( dset[1] );               
               free_f3tensor(amp,0,NY+1,0,NX+1,0,NZ+1);
               break;

        } 
        status = H5Sclose( mspace );
        status = H5Pclose( plist );
        status = H5Fclose( fid );

}


