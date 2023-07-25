cmake_minimum_required(VERSION 3.0)
PROJECT(master Fortran C)
SET(CMAKE_Fortran_COMPILER mpifort)#nvfortran
SET(CMAKE_C_COMPILER nvc)

SET(CMAKE_C_FLAGS "-O0 -g  -m64 -DPOINTER_64 -DINTEGER_IS_INT -DREAL_8 -DREAL_BIGGER_THAN_INTERGER -DLITTLE_ENDIAN  -DLINUX ")
#SET(CMAKE_Fortran_FLAGS " -g -traceback -DNDEBUG  -fPIC  -m64 -r8  -Ktrap=fp -Minfo=all -fast -Minline -acc=gpu -Mcuda -Mlarge_arrays -mcmodel=medium")
#SET(CMAKE_Fortran_FLAGS "   -fPIC  -m64 -r8  -Ktrap=fp -Minfo=all -fast -Minline -acc=gpu -Mcuda -Mlarge_arrays -mcmodel=medium -ta=tesla")
SET(CMAKE_Fortran_FLAGS " -g   -fPIC  -m64 -r8  -Ktrap=fp -Minfo=all -fast -Minline -acc=gpu -Mcuda -Mlarge_arrays -mcmodel=medium -ta=nvidia  -Mchkptr")

#set mpi
SET(MPI_HOME =/opt/nvidia/hpc_sdk/Linux_x86_64/22.9/comm_libs/openmpi/openmpi-3.1.5)
SET(MPI_INCLUDE_PATH  ${MPI_HOME}/lib)
#SET(MPI_LIBRARIES ${MPI_HOME}/lib/libmpi_usempif08.so ${MPI_HOME}/lib/libmpi.so)# ${MPI_HOME}/lib/libmpi.so)
#set for eccodes
SET(ECCODES_HOME /public/home/nwpyewu/app/nv/eccodes)
SET(GRIBAPI_INCLUDE_PATH ${ECCODES_HOME}/include)
SET(GRIBLIB ${ECCODES_HOME}/lib64/libeccodes_f90.so ${ECCODES_HOME}/lib64/libeccodes.so)

INCLUDE_DIRECTORIES( interface ${MPI_INCLUDE_PATH} ${GRIBAPI_INCLUDE_PATH}) 
SET(WAM_HAVE_UNWAM FALSE)
SET(HAVE_SINGLE_PRECISION FALSE)

#for GPU 
ADD_DEFINITIONS(-DWAM_PHYS_GPU)


### Set third party libraries
if(WAM_HAVE_UNWAM)
    SET(NORMAL_FLASGS "${CMAKE_Fortran_FLAGS} -DWAM_HAVE_UNWAM")
endif()

if( HAVE_SINGLE_PRECISION )
  set( prec sp )
  SET(NORMAL_FLASGS "${CMAKE_Fortran_FLAGS} -DPARKIND1_SINGLE_NEMO")
else()
  set( prec dp )
endif()


### Sources
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D__FILENAME__='\"$(notdir $(abspath $<))\"'")

FILE(GLOB unwam_srcs src/ecwam/unwam/*.F90)
FILE(GLOB ecwam_srcs src/ecwam/*.F90 src/util/*.F90 src/module/*.F90 src/util/*.c)
FILE(GLOB ecwam_exe_srcs src/programs/*.F90)

if(WAM_HAVE_UNWAM)
    LIST(APPEND ecwam_srcs ${unwam_srcs})
endif()

ADD_LIBRARY(ecwam STATIC ${ecwam_srcs})

#add exe files
foreach(tool ${ecwam_exe_srcs})
    STRING(REGEX REPLACE ".+/(.+)\\..*" "\\1" toolname ${tool})
    ADD_EXECUTABLE(${toolname} ${tool})
    TARGET_LINK_LIBRARIES(${toolname} ecwam  ${MPI_LIBRARIES} ${GRIBLIB} )
endforeach()
