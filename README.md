ecwam_gpu
## Introduction
GPU version of ecWAM(https://github.com/ecmwf-ifs/ecwam). Both the propagation and physical part are rewrited by OpenAcc directives.

## Requirements
1. nvhpc22.1
2. CMake
3. eccodes
The project are reorgnized by cmake, and some packages such as ecbuild and fiat are no longer needed.

## Building ecwam_acc

1. SET(GRIBAPI_INCLUDE_PATH pathtoeccodes/include)
2. SET(GRIBLIB pathtoeccodes/lib/libeccodes_f90.so pathtoeccodes/lib/libeccodes.so)
3. make
## Running ecwam_acc
If some  wave initial condition files such as  BLS\*, LAW\* are not prepared, the corresponding field in code are automatically setted by default value.
For more  conditions see https://github.com/ecmwf-ifs/ecwam.


