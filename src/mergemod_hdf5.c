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

void mergemod_hdf5(char modfile[STRING_SIZE]){

     extern int MYID, NPROCX, NPROCY, NPROCZ, NXG, NYG, NZG, NX, NY, NZ, NPROC, IDX, IDY, IDZ;
     extern FILE *FP;

     char file[STRING_SIZE], filename[STRING_SIZE];
     FILE *fp[NPROCY_MAX][NPROCX_MAX][NPROCZ_MAX];
     int i, j, k, ip, jp, kp;
     float *a;

     hid_t fpout, dset, dspace, mspace;
     hsize_t dims[3], offset[3], count[3], block[3], stride[3];
     herr_t status;

     if ((NPROCX>NPROCX_MAX)||(NPROCY>NPROCY_MAX)||(NPROCZ>NPROCZ_MAX))
	err(" mergemod_hdf5.c: constant expression NPROC?_MAX < NPROC? ");
    
     sprintf( filename, "%s.h5", modfile ); 
     fprintf( FP, " writing merged model file to  %s \n", filename );

     fpout = H5Fcreate( filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

     for (kp=0;kp<=NPROCZ-1; kp++) {
     for (ip=0;ip<=NPROCX-1; ip++) {
     for (jp=0;jp<=NPROCY-1; jp++) {
      	sprintf(file,"%s.%i.%i.%i",modfile,ip,jp,kp);
      	fp[jp][ip][kp]=fopen(file,"r");
      	if (fp[jp][ip][kp]==NULL) err("mergemod_hdf5.c: can't read modfile !"); 
     }}}

     fprintf(FP," Copying...");

     dims[0] = NYG;
     dims[1] = NXG;
     dims[2] = NZG;
     dspace = H5Screate_simple( 3, dims, NULL );
     dset = H5Dcreate2( fpout, "data", H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

     a = (float *) malloc( NY*sizeof(float) );

     for ( k=0; k<3; k++ ) {
         stride[k] = 1;
         block[k] = 1;
        count[k] = 1;
     }
     count[0] = NY;
     mspace = H5Screate_simple(3, count, NULL);

     for (kp=0;kp<=NPROCZ-1; kp++) {
     for (k=1;k<=NZ;k+=IDZ) {
     for (ip=0;ip<=NPROCX-1; ip++) {
     for (i=1;i<=NX;i+=IDX) {
     for (jp=0;jp<=NPROCY-1; jp++) {

        fread( a, sizeof(float), NY, fp[jp][ip][kp] );	

        offset[0] = jp*NY;
        offset[1] = ip*NX;
        offset[2] = kp*NZ;
        status = H5Sselect_hyperslab( dspace, H5S_SELECT_SET, offset, stride, count, block );

        status = H5Dwrite( dset, H5T_NATIVE_FLOAT, mspace, dspace, H5P_DEFAULT, a );

     }}}}}
     status = H5Sclose( mspace );

     free( a );
     a = NULL;

     fprintf(FP," ... finished. \n");

     for (kp=0;kp<=NPROCZ-1; kp++) {
     for (ip=0;ip<=NPROCX-1; ip++) {
     for (jp=0;jp<=NPROCY-1; jp++) {
     	fclose(fp[jp][ip][kp]);
     }}}
      	
     status = H5Dclose( dset );
     status = H5Sclose( dspace );
     status = H5Fclose( fpout );

}

