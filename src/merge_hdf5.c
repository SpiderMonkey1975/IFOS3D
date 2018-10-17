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
 *   merge snapshots files written by the different processes to 
 *   a single file                                 
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include <hdf5.h>
#include <unistd.h>

void merge_hdf5( int nsnap, int type ) {

     extern char SNAP_FILE[STRING_SIZE];
     extern int NXG, NYG, SNAP_FORMAT, NPROCX, NPROCY, NPROCZ, NX, NY, NZ, IDX, IDY, IDZ, NXG, NYG, NZG;
     extern FILE *FP;
	
     char file[STRING_SIZE], mfile[STRING_SIZE], outfile[STRING_SIZE], ext[10], dsetname[7];
     FILE *fp[NPROCY_MAX][NPROCX_MAX][NPROCZ_MAX];
     int i, j, k, ip, jp, kp, n;
     float *a;

     hid_t fpout, dset, dspace, mspace;
     hsize_t dims2[3], offset2[3], count2[3], block2[3], stride2[3];
     hsize_t dims[4], offset[4], count[4], block[4], stride[4];
     herr_t status;

     if ((NPROCX>NPROCX_MAX)||(NPROCY>NPROCY_MAX)||(NPROCZ>NPROCZ_MAX))
	err(" merge_hdf5.c: constant expression NPROC?_MAX < NPROC? ");

 /**
  ** Set the name of the output HDF5 file
  **/

     sprintf(mfile,"%s%s",SNAP_FILE,ext);
     fprintf(FP," (files: %s.??? ).\n",mfile);
	
     sprintf(ext,".h5");
     sprintf(outfile,"%s%s",SNAP_FILE,ext);
     fprintf(FP,"\n writing merged snapshot file to  %s \n",outfile);

 /**
  ** Create or open the HDF5 output file
  **/

     if ( access(outfile, F_OK) ) {
        fpout = H5Fopen( outfile, H5F_ACC_RDWR, H5P_DEFAULT );
     } else {
        fpout = H5Fcreate( outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
     }

 /**
  ** Add the appropriate dataset
  **/

     if ( nsnap>0 ) {
        dims[0] = nsnap;
        dims[1] = NZG;
        dims[2] = NXG;
        dims[3] = NYG;
        dspace = H5Screate_simple( 4, dims, NULL );
     } else {
        dims2[0] = NZG;
        dims2[1] = NXG;
        dims2[2] = NYG;
        dspace = H5Screate_simple( 3, dims2, NULL );
     }

     switch (type){
            case 1: 
               fprintf(FP," x-component of particle velocity");
               sprintf(dsetname,"x");
               break;
            case 2: 
               fprintf(FP," y-component of particle velocity");
               sprintf(dsetname,"y");
               break;
            case 3: 
               fprintf(FP," z-component of particle velocity");
               sprintf(dsetname,"z");
               break;
            case 4: 
               fprintf(FP," P-wave energyfield");
               sprintf(dsetname,"div");
               break;
            case 5: 
	       fprintf(FP," S-wave energyfield");
               sprintf(dsetname,"rot");
	       break;
      	    case 6: 
	       fprintf(FP," pressure");
               sprintf(dsetname,"p");
	       break;
	    case 7:
	       fprintf(FP," Gradient lambda");
               sprintf(dsetname,"grad1");
	       break;
   	    case 8:
	       fprintf(FP," Gradient mu");
               sprintf(dsetname,"grad2");
	       break;    
	    case 9:
	       fprintf(FP," Gradient rho");
               sprintf(dsetname,"grad3");
	       break;        
      	    default: 
	       err(" merge: nonexistant dataset specified! ");
	       break;
     }
     dset = H5Dcreate2( fpout, dsetname, H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

 /**
  ** Combine the contents of all the raw binary input files into the HDF5 file
  **/

     fprintf(FP," Combining snapshot files: %s.??? ",mfile);
	
     for (kp=0;kp<=NPROCZ-1; kp++) {
     for (ip=0;ip<=NPROCX-1; ip++) {
     for (jp=0;jp<=NPROCY-1; jp++) {

      	 sprintf(file,"%s.%i%i%i",mfile,ip,jp,kp);
      	 fp[jp][ip][kp]=fopen(file,"r");
      	 if (fp[jp][ip][kp]==NULL) err("merge: can't read snapfile !"); 

     }}}

     a = (float *) malloc( NZ*sizeof(float) );

     if ( nsnap>0 ) {
        for ( n=0; n<4; n++ ) {
            stride[n] = 1;
            block[n] = 1;
            count[n] = 1;
        }
        count[1] = NZ; 
        mspace = H5Screate_simple(4, count, NULL);

        for (n=0;n<=nsnap; n++) {/*for gradient calculation?*/
        for (jp=0;jp<=NPROCY-1; jp++) {
        for (j=1;j<=NY;j+=IDY) {
        for (ip=0;ip<=NPROCX-1; ip++) {
        for (i=1;i<=NX;i+=IDX) {
        for (kp=0;kp<=NPROCZ-1; kp++) {

            fread( a, sizeof(float), NZ, fp[jp][ip][kp] );

            offset[0] = n;
            offset[1] = kp*NZ;
            offset[2] = jp*NY;
            offset[3] = ip*NX;
            status = H5Sselect_hyperslab( dspace, H5S_SELECT_SET, offset, stride, count, block );

            status = H5Dwrite( dset, H5T_NATIVE_FLOAT, mspace, dspace, H5P_DEFAULT, a );     
            status = H5Sclose( mspace );

        }}}}}}

     } else {
        for ( n=0; n<3; n++ ) {
            stride2[n] = 1;
            block2[n] = 1;
            count2[n] = 1;
        }
        count2[0] = NZ; 
        mspace = H5Screate_simple(3, count2, NULL);

        for (jp=0;jp<=NPROCY-1; jp++) {
        for (j=1;j<=NY;j+=IDY) {
        for (ip=0;ip<=NPROCX-1; ip++) {
        for (i=1;i<=NX;i+=IDX) {
        for (kp=0;kp<=NPROCZ-1; kp++) {

            fread( a, sizeof(float), NZ, fp[jp][ip][kp] );

            offset2[0] = kp*NZ;
            offset2[1] = jp*NY;
            offset2[2] = ip*NX;
            status = H5Sselect_hyperslab( dspace, H5S_SELECT_SET, offset2, stride2, count2, block2 );

            status = H5Dwrite( dset, H5T_NATIVE_FLOAT, mspace, dspace, H5P_DEFAULT, a );     
            status = H5Sclose( mspace );

        }}}}}
     }

     free( a );
     a = NULL;

     for (kp=0;kp<=NPROCZ-1; kp++) {
     for (ip=0;ip<=NPROCX-1; ip++) {
     for (jp=0;jp<=NPROCY-1; jp++) {
         fclose(fp[jp][ip][kp]);
     }}}

 /**
  ** Close the new HDF5 output file
  **/

     status = H5Dclose( dset );
     status = H5Sclose( dspace );
     status = H5Fclose( fpout );

}

