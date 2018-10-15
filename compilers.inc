
ifeq ($(ALLINEA), 1)
  PROFILE=1
endif   
ifeq ($(CRAYPAT), 1)
  PROFILE=1
endif   

##-------------------------------------------------------
## CRAY-specific compiler settings (default)
##-------------------------------------------------------

CXX=CC
CC=cc

CFLAGS=-hPIC -O3 -h fp3
CXXFLAGS=-hPIC -O3 -h gnu -h fp3  

ifeq ($(DEBUG),1)
  CFLAGS=-hPIC -g -h msglevel_3 -O0 -h bounds
  CXXFLAGS=-hPIC -g -h msglevel_3 -O0 -h bounds
endif
ifeq ($(PROFILE), 1)
  CFLAGS+= -G2
  CXXFLAGS+= -G2
endif

##-------------------------------------------------------
## INTEL-specific compiler settings 
##-------------------------------------------------------

ifeq ($(INTEL_COMPILER),1)

ifeq ($(NONCRAY),1)
CC=mpiicc
CXX=mpiicpc
endif

CFLAGS=-fPIC -fast -no-ipo
CXXFLAGS=-fPIC -fast -no-ipo

ifeq ($(DEBUG),1)
  CFLAGS=-fPIC -g -w2 -O0 -check=stack -traceback
  CXXFLAGS=-fPIC -g -w2 -O0 -check=stack -traceback
endif

ifeq ($(PROFILE), 1)
  CFLAGS+= -debug extended
  CXXFLAGS+= -debug extended
endif
endif

##-------------------------------------------------------
## GNU-specific compiler settings
##-------------------------------------------------------

ifeq ($(GNU_COMPILER),1)

ifeq ($(NONCRAY),1)
CC=mpicc
CXX=mpicxx
endif

CFLAGS=-fPIC -O3 -mfma -mavx  
CXXFLAGS=-fPIC -O3 -mfma -mavx -std=c++11

ifeq ($(DEBUG),1)
  CFLAGS=-fPIC -g -Wall -O0 -check=stack -traceback
  CXXFLAGS=-fPIC -g -Wall -O0 -check=stack -traceback
endif

ifeq ($(PROFILE), 1)
  CFLAGS+= -g
  CXXFLAGS+= -g
endif
endif

