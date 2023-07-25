# ecwam的cmake编译方式
## ecWAM的cmake版本，代码同开源ecwam一致，但是不用ecbuild和fiat等辅助库，可生产ecwam库文件和可执行文件
## 编译器采用的时ifort和icc,运行前需要再CMakeList.txt中设置mpi库和eccodes库
###设置mpi库
SET(MPI_INCLUDE_PATH  /vol7/app_share/mpich/intel/3.3.2/include)
SET(MPI_LIBRARIES /vol7/app_share/mpich/intel/3.3.2/lib/libmpifort.so /vol7/app_share/mpich/intel/3.3.2/lib/libmpi.so)
###设置eccodes库
SET(GRIBAPI_INCLUDE_PATH /vol7/app_share/eccodes-2.12.0/include)
SET(GRIBLIB /vol7/app_share/eccodes-2.12.0/lib/libeccodes_f90.so /vol7/app_share/eccodes-2.12.0/lib/libeccodes.so)

杨锦辉 20230413
