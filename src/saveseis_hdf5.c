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
 *   write seismograms to files using parallel HDF5 
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include <unistd.h>
#include <hdf5.h>

void saveseis_hdf5( float **sectionvx, float **sectionvy, float **sectionvz, float **sectionp, float **sectioncurl, 
                    float **sectiondiv, int **recpos, int **recpos_loc, int ntr_glob, int ntr, float **srcpos, 
                    int ishot, int ns, int obs, int iteration, int nshot ) { 
		
     extern int SEISMO, MYID, RUN_MULTIPLE_SHOTS;	
     extern char  SEIS_FILE[STRING_SIZE], SEIS_OBS_FILE[STRING_SIZE];

     char filename[STRING_SIZE], grpname[10], dsetname[10];
     float *buf2,*buf1, **srcpos1=NULL;
     int nsrc=1,i,k,m,l,nt;

     hid_t gid, plist, dset, dspace, mspace, fid;
     herr_t status;
     hsize_t dims[3], offset[3], stride[3], count[3], block[3];

  /* 
   * Check if the appropriate HDF5 seismogram data file already exists
   *-------------------------------------------------------------------*/		
     if (obs==0) sprintf( filename, "%s_it%d.h5", SEIS_FILE, iteration );
     if (obs==1) sprintf( filename, "%s_it%d.h5", SEIS_OBS_FILE, iteration );

     if ( (access( filename, F_OK ) == -1) && (MYID==0) ) {
        fid = H5Fcreate( filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
        for ( i=1; i<=nshot; i++ ) {
            sprintf( grpname, "Shot_%d", i );
            gid = H5Gcreate2(fid, grpname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
            status = H5Gclose( gid );
        }
        status = H5Fclose( fid );
     }
     i = MPI_Barrier( MPI_COMM_WORLD );

  /* 
   * All MPI tasks now open the HDF5 seismogram data file  
   *-------------------------------------------------------------------*/		
     plist = H5Pcreate( H5P_FILE_ACCESS );
     H5Pset_fapl_mpio( plist, MPI_COMM_WORLD, MPI_INFO_NULL );
     fid = H5Fopen( filename, H5F_ACC_RDWR, plist );
     status = H5Pclose( plist );

     sprintf( grpname, "Shot_%d", ishot );
     gid = H5Gopen2( fid, grpname, H5P_DEFAULT );

  /* 
   * Check if any receiver data is available
   *------------------------------------------------------------------*/		
//     if (ntr<1) {
	
     srcpos1 = fmatrix(1,7,1,1);
     for (nt=1;nt<=7;nt++) srcpos1[nt][1]=srcpos[nt][ishot];

  /* 
   * Write out seismogram data 
   *------------------------------------------------------------------*/		
     switch (SEISMO){
	case 1 : /* particle velocities only */
	     if (MYID==0 ) { printf("\n\n writing velocity data seismogramtraces to %s \n",filename); }
             sprintf( dsetname, "vx" );
             outseis_hdf5( gid, dsetname, sectionvx, recpos_loc, ntr_glob, ntr, ns );
             sprintf( dsetname, "vy" );
             outseis_hdf5( gid, dsetname, sectionvy, recpos_loc, ntr_glob, ntr, ns );
             sprintf( dsetname, "vz" );
             outseis_hdf5( gid, dsetname, sectionvz, recpos_loc, ntr_glob, ntr, ns );
//	     outseis_hdf5( gid, 1, sectionvx, recpos, recpos_loc, ntr_glob, ntr, srcpos1, nsrc, ns );
//	     outseis_hdf5( gid, 2, sectionvy, recpos, recpos_loc, ntr_glob, ntr, srcpos1, nsrc, ns );
//	     outseis_hdf5( gid, 3, sectionvz, recpos, recpos_loc, ntr_glob, ntr, srcpos1, nsrc, ns );
	     break;

	case 2 : /* pressure only */
	     if (MYID==0 ) { 
                printf("\n\n writing seismogramtraces of pressure data to %s \n",filename);
             }
//	     outseis_hdf5( gid, 0, sectionp, recpos,recpos_loc, ntr_glob, ntr, srcpos1, nsrc, ns );
	     break;

	case 3 : /* curl and div only */
	     if (MYID==0 ) { 
                printf("\n\n writing seismogramtraces of curl and div data to %s \n",filename);
             }
//	     outseis_hdf5( gid, 0,  sectiondiv, recpos, recpos_loc, ntr_glob, ntr, srcpos1, nsrc, ns );
//	     outseis_hdf5( gid, 0, sectioncurl, recpos, recpos_loc, ntr_glob, ntr, srcpos1, nsrc, ns );	
	     break;	

	case 4 : /* everything */
	     if (MYID==0 ) { 
                printf("\n\n ERROR: option 4 in saveseis_hdf5 not a valid value. Make multiple calls to saveseis_hdf5 with options 1,2 and 3 instead.\n\n");
                k = 0;
                i = MPI_Abort( MPI_COMM_WORLD, k );
             }
	     break;
     }

     status = H5Gclose( gid );
     status = H5Fclose( fid );

     free_matrix(srcpos1,1,7,1,1);


//     status = H5Fclose( fid );
//   }	
/*	MPI_Barrier(MPI_COMM_WORLD);
	if(MYID==0){
		sprintf(outfile,"%s_vx_it%d.%s.shot%d",seisfile,iteration,file_ext,ishot);
		fout=fopen(outfile,"w");
		buf2=vector(1,ns);
		buf1=vector(1,240/sizeof(float));	
			for(i=0;i<=NPROC-1;i++){
				k=0;
				sprintf(infile,"%s_vx_it%d.%s.shot%d.%d",seisfile,iteration,file_ext,ishot,i);
				fin=fopen(infile,"r");
				if(fin!=0){
					fseek(fin,0,SEEK_END);
					k=ftell(fin);
					if(k%(240+4*ns)==0)l=k/(240+4*ns);
					else{ fprintf(FP,"attention: Error in seisfile");
						l=0;
					}
					fseek(fin,0,SEEK_SET);
					
					for(m=1;m<=l;m++){
					fread(buf1,sizeof(float),240/sizeof(float),fin);
					fread(buf2,sizeof(float),ns,fin);
							
					fwrite(buf1,sizeof(float),240/sizeof(float),fout);
					fwrite(buf2,sizeof(float),ns,fout);
					}
					fclose(fin);
					
				}
			}
		fclose(fout);
		free_vector(buf2,1,240/sizeof(float));
		free_vector(buf1,1,ns);
	}

	if(MYID==0){
		sprintf(outfile,"%s_vy_it%d.%s.shot%d",seisfile,iteration,file_ext,ishot);
		fout=fopen(outfile,"w");
		buf2=vector(1,ns);
		buf1=vector(1,240/sizeof(float));
		k=0;
		
		for(i=0;i<=NPROC-1;i++){
			sprintf(infile,"%s_vy_it%d.%s.shot%d.%d",seisfile,iteration,file_ext,ishot,i);
			fin=fopen(infile,"r");
			if(fin!=0){
				fseek(fin,0,SEEK_END);
				k=ftell(fin);
				if(k%(240+4*ns)==0)l=k/(240+4*ns);
				else{ fprintf(FP,"attention: Error in seisfile");
					l=0;
				}
				fseek(fin,0,SEEK_SET);
				
				for(m=1;m<=l;m++){
				fread(buf1,sizeof(float),240/sizeof(float),fin);
				fread(buf2,sizeof(float),ns,fin);
						
				fwrite(buf1,sizeof(float),240/sizeof(float),fout);
				fwrite(buf2,sizeof(float),ns,fout);
				}
				fclose(fin);
			}
		}
		fclose(fout);
		free_vector(buf2,1,240/sizeof(float));
		free_vector(buf1,1,ns);
	}
	
	if(MYID==0){
		sprintf(outfile,"%s_vz_it%d.%s.shot%d",seisfile,iteration,file_ext,ishot);
		fout=fopen(outfile,"w");
		buf2=vector(1,ns);
		buf1=vector(1,240/sizeof(float));
		k=0;
		for(i=0;i<=NPROC-1;i++){
			sprintf(infile,"%s_vz_it%d.%s.shot%d.%d",seisfile,iteration,file_ext,ishot,i);
			fin=fopen(infile,"r");
			if(fin!=0){
				fseek(fin,0,SEEK_END);
				k=ftell(fin);
				if(k%(240+4*ns)==0)l=k/(240+4*ns);
				else{ fprintf(FP,"attention: Error in seisfile");
					l=0;
				}
				fseek(fin,0,SEEK_SET);
				for(m=1;m<=l;m++){
				fread(buf1,sizeof(float),240/sizeof(float),fin);
				fread(buf2,sizeof(float),ns,fin);
						
				fwrite(buf1,sizeof(float),240/sizeof(float),fout);
				fwrite(buf2,sizeof(float),ns,fout);
				}
				fclose(fin);
				
			}
		}
		fclose(fout);
		free_vector(buf2,1,240/sizeof(float));
		free_vector(buf1,1,ns);
	}
	free_matrix(srcpos1,1,7,1,1);
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(iteration>0&&ntr>0&&METHOD&&!obs){
		sprintf(infile,"%s_vx_it%d.%s.shot%d.%d",seisfile,iteration,file_ext,ishot,MYID);
		remove(infile);
		sprintf(infile,"%s_vy_it%d.%s.shot%d.%d",seisfile,iteration,file_ext,ishot,MYID);
		remove(infile);
		sprintf(infile,"%s_vz_it%d.%s.shot%d.%d",seisfile,iteration,file_ext,ishot,MYID);
		remove(infile);
	} */
}
