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

/*------------------------------------------------------------
 * output of model parameters vp,vs and rho to MOD_OUT_FILE;
 * S. Dunkl 2013 
 -------------------------------------------------------------*/



#include "fd.h"

void outmod(float ***rho, float ***pi,float ***u, int iteration){

	extern int POS[4];
	extern char  MOD_OUT_FILE[STRING_SIZE];
	FILE *fpmod1, *fpmod2, *fpmod3;
	extern int MYID, NX, NY, NZ;

	int i,j,k;
	char modfile1[STRING_SIZE],modfile2[STRING_SIZE],modfile3[STRING_SIZE], modfile4[STRING_SIZE];
	char modfile5[STRING_SIZE],modfile6[STRING_SIZE],hdf5_prefix[STRING_SIZE];
	float ***vp, ***vs;
#ifdef HDF5
        char dsetname[10];
#endif

        vs = f3tensor(0,NY+1,0,NX+1,0,NZ+1);
        vp = f3tensor(0,NY+1,0,NX+1,0,NZ+1);

        for (k=1;k<=NZ;k++){
        for (i=1;i<=NX;i++){
        for (j=1;j<=NY;j++){
            vp[j][i][k] = sqrt(pi[j][i][k]/rho[j][i][k]);
            vs[j][i][k] = sqrt(u[j][i][k]/rho[j][i][k]);
        }}}


#ifdef HDF5
        sprintf( hdf5_prefix, "%s_iteration%d", MOD_OUT_FILE, iteration );

        sprintf( dsetname, "vp" );
	mergemod_hdf5( hdf5_prefix, dsetname, vp );
        sprintf( dsetname, "vs" );
	mergemod_hdf5( hdf5_prefix, dsetname, vs );
        sprintf( dsetname, "rho" );
	mergemod_hdf5( hdf5_prefix, dsetname, rho );
#endif
	
	sprintf(modfile1,"%s.vp_it%d.%i.%i.%i",MOD_OUT_FILE,iteration,POS[1],POS[2],POS[3]);
	sprintf(modfile2,"%s.vs_it%d.%i.%i.%i",MOD_OUT_FILE,iteration,POS[1],POS[2],POS[3]);
	sprintf(modfile3,"%s.rho_it%d.%i.%i.%i",MOD_OUT_FILE,iteration,POS[1],POS[2],POS[3]);
	sprintf(modfile4,"%s.vp_it%d",MOD_OUT_FILE,iteration);
	sprintf(modfile5,"%s.vs_it%d",MOD_OUT_FILE,iteration);
	sprintf(modfile6,"%s.rho_it%d",MOD_OUT_FILE,iteration);

	fpmod1=fopen(modfile1,"w");
	fpmod2=fopen(modfile2,"w");
	fpmod3=fopen(modfile3,"w");
	
	for (k=1;k<=NZ;k++){
        for (i=1;i<=NX;i++){	
	for (j=1;j<=NY;j++){
	    fwrite(&vp[j][i][k], sizeof(float), 1,fpmod1);
	    fwrite(&vs[j][i][k], sizeof(float), 1,fpmod2);
	    fwrite(&rho[j][i][k], sizeof(float), 1,fpmod3);
	}}}

	fclose(fpmod1);
	fclose(fpmod2);
	fclose(fpmod3);
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(MYID==0){
		mergemod(modfile4,3);
		mergemod(modfile5,3);
		mergemod(modfile6,3);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	remove(modfile1);
	remove(modfile2);
	remove(modfile3);

        free_f3tensor(vs,0,NY+1,0,NX+1,0,NZ+1);
        free_f3tensor(vp,0,NY+1,0,NX+1,0,NZ+1);
	
}
