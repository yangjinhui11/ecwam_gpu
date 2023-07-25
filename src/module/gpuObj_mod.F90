module GPUOBJ_MOD
    USE YOWFIELD_MOD, ONLY : FREQUENCY_FIELD, ENVIRONMENT_FIELD,FORCING_FIELDS_FIELD, &
 &  WAVE2OCEAN_FIELD, INTGT_PARAM_FIELDS_FIELD, SOURCE_CONTRIBS_FIELD
USE YOWDRVTYPE  , ONLY : WVGRIDGLO


implicit none
SAVE
! Objects to store fields for wamintgr
TYPE(FREQUENCY_FIELD) :: WVPRPT_FIELD
TYPE(ENVIRONMENT_FIELD) :: WVENVI_FIELD
TYPE(FORCING_FIELDS_FIELD) :: FF_NOW_FIELD
TYPE(FORCING_FIELDS_FIELD) :: FF_NEXT_FIELD
TYPE(WAVE2OCEAN_FIELD) :: WAM2NEMO_FIELD
TYPE(INTGT_PARAM_FIELDS_FIELD) :: INTFLDS_FIELD
TYPE(SOURCE_CONTRIBS_FIELD) :: SRC_CONTRIBS

INTEGER:: COMM_SHARED
!Object for propag_wam
!TYPE(FIELD_3D_WRAPPER):: FL1_FIELD,FL3_FIELD
!TYPE(FIELD_2D_WRAPPER):: WAVNUM_FIELD,CGROUP_FIELD,OMOSNH2KD_FIELD

contains

  subroutine get_device_info(IU06)
    USE openacc

    integer,intent(in):: IU06
    integer:: num_threads, num_gpus
    num_gpus = acc_get_num_devices(acc_device_nvidia)
    WRITE(IU06,*)"Total devices ",num_gpus,"used"
  endsubroutine
  subroutine UPDATE_TOHOST(WVPRPT_FIELD,WVENVI_FIELD,FF_NOW_FIELD,WAM2NEMO_FIELD,INTFLDS_FIELD,SRC_CONTRIBS)
    TYPE(FREQUENCY_FIELD),OPTIONAL:: WVPRPT_FIELD
    TYPE(ENVIRONMENT_FIELD),OPTIONAL :: WVENVI_FIELD
    TYPE(FORCING_FIELDS_FIELD),OPTIONAL :: FF_NOW_FIELD
  
    TYPE(WAVE2OCEAN_FIELD),OPTIONAL :: WAM2NEMO_FIELD
    TYPE(INTGT_PARAM_FIELDS_FIELD),OPTIONAL :: INTFLDS_FIELD
    TYPE(SOURCE_CONTRIBS_FIELD),OPTIONAL :: SRC_CONTRIBS
    IF(PRESENT(FF_NOW_FIELD))THEN
        CALL FF_NOW_FIELD%ENSURE_HOST(LAIRD=.TRUE.,LWDWAVE=.TRUE.,LCICOVER=.TRUE.,LWSWAVE=.TRUE., &
            &LWSTAR=.TRUE., LUFRIC=.TRUE., LTAUW=.TRUE., LTAUWDIR=.TRUE., &
            &LZ0M=.TRUE., LZ0B=.TRUE., LCHRNCK=.TRUE., LCITHICK=.TRUE.)
    ENDIF
    IF(PRESENT(INTFLDS_FIELD))THEN
        CALL INTFLDS_FIELD%ENSURE_HOST(LWSEMEAN=.TRUE.,LWSFMEAN=.TRUE.,LUSTOKES=.TRUE.,LVSTOKES=.TRUE.,&
        &LSTRNMS=.TRUE.,LTAUXD=.TRUE.,LTAUYD=.TRUE.,LTAUOCXD=.TRUE., &
        &LTAUOCYD=.TRUE.,LTAUOC=.TRUE.,LPHIOCD=.TRUE.,LPHIEPS=.TRUE., &
        &LPHIAW=.TRUE.)
    ENDIF
    IF(PRESENT(SRC_CONTRIBS))THEN
        CALL SRC_CONTRIBS%ENSURE_HOST(LFL1=.TRUE.,LXLLWS=.TRUE.,LMIJ=.TRUE.)
    ENDIF
  endsubroutine
  subroutine init_devices(myrank,IU06)
    use openacc
    use cudafor
    use mpi
    integer,intent(in):: IU06,MYRANK
    integer:: num_gpus,ierr,nprocsh,iprocsh,nlen
    integer(kind=cuda_count_kind):: heapsize
    character(20):: hostname
    !create a shared communicator for all ranks in the node
    CALL MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,COMM_SHARED,ierr)
    ! get the total number of processors in node
    CALL MPI_Comm_size(COMM_SHARED,nprocsh,ierr)
    ! get the index of the local processor in the local node
    CALL MPI_Comm_rank(COMM_SHARED,iprocsh,ierr)
    CALL MPI_GET_PROCESSOR_NAME(hostname,nlen,ierr)

    num_gpus = acc_get_num_devices(acc_device_nvidia)
    ! set the gpu device number based on local rank
    call acc_set_device_num(mod(iprocsh,num_gpus),acc_device_nvidia)
    !call acc_set_device_num(mod(MYRANK-1,4)+1,acc_device_nvidia)
    WRITE(IU06,*)"Local size",nprocsh
    WRITE(IU06,*)"Global rank ",MYRANK,"Node ",hostname(:nlen)," local rank ",iprocsh
    WRITE(IU06,*)" Devices local node ",num_gpus,"used, current ID is",acc_get_device_num(acc_device_nvidia)
    ! set heap size

    heapsize=16_8*1024_8*1024_8*10
    ierr=cudaDeviceSetLimit(cudaLimitMallocHeapSize,heapsize)
  endsubroutine
  
  subroutine get_mem_info(IU06)
    use openacc
    use cudafor
    use accel_lib
    integer,intent(in):: IU06
    integer :: devNum,devType,freeMem,totalMem,freeMem1,allocMem
    INTEGER(C_SIZE_T) :: free, total
    character(14):: devName,devVendor
    integer(C_SIZE_T),save :: last_free=0

    devType = acc_get_device_type()
    devNum = acc_get_device_num(devType)
    !freeMem = acc_get_property(devNum,acc_device_current,acc_property_free_memory)
    !totalMem = acc_get_property(devNum,acc_device_current,acc_property_memory)
    call acc_get_property_string(devNum,acc_device_current,acc_property_name,devName)
    call acc_get_property_string(devNum,acc_device_current,acc_property_vendor,devVendor)
    
    !WRITE(*,'(I2,2A10,3(A10,F5.2)')devNum,trim(devName),trim(devVendor),"Alloced mem(M)",acc_bytesalloc()/1024.0/1024.0," Free mem(M),",freeMem/1024./1024.,"Total mem",totalMem/1024./1024.
    freeMem1 = cudamemgetinfo(free,total)
    
    WRITE(IU06,'((1X,I2),(1X,A12),2(1X,A12,F9.3),(1X,A12,F10.3))')devNum,trim(devName),"Alloc mem(M)",acc_bytesalloc()/1024.0/1024.0," Free mem(M),",free/1024./1024.,"Inc mem",(last_free-free)/1024./1024.
    !call acc_present_dump
    last_free=free
  endsubroutine
endmodule GPUOBJ_MOD
