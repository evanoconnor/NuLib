!-*-f90-*-
subroutine nulibtable_reader(filename)
  
  use nulibtable
  use hdf5
  implicit none

  character(*) :: filename
  
  !H5 stuff
  integer :: error,rank,cerror
  integer(HID_T) :: file_id,dset_id,dspace_id
  integer(HSIZE_T) :: dims1(1), dims2(2), dims3(3), dims4(4), dims5(5)!, etc....

  !local
  real*8, allocatable :: nulibtable_temp(:,:,:,:,:)
  real*8 :: timestamp
  integer :: startindex,endindex
  integer :: i,j

  cerror = 0

  !open HDF5 file, given filename                                                                                             
  call h5open_f(error)
  if (error.ne.0) then
     stop "Error reading in nulib file"
  endif
  
  call h5fopen_f(trim(adjustl(filename)), &
       H5F_ACC_RDONLY_F,file_id,error)
  if (error.ne.0) then
     write(*,*) trim(adjustl(filename))
     stop "Error reading in nulib table"
  endif

  !write scalars (rank=1, dims1(1) = 1)
  rank = 1
  dims1(1) = 1

  !first lets read number of species and number of groups and number of rho,temp, and ye points, also timestamp
  call h5dopen_f(file_id, "number_species",dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nulibtable_number_species, dims1, error)
  call h5dclose_f(dset_id, error)
  cerror = cerror + error

  call h5dopen_f(file_id, "timestamp",dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, timestamp, dims1, error)
  call h5dclose_f(dset_id, error)
  cerror = cerror + error  

  call h5dopen_f(file_id, "number_groups",dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nulibtable_number_groups, dims1, error)
  call h5dclose_f(dset_id, error)
  cerror = cerror + error

  allocate(nulibtable_energies(nulibtable_number_groups))
  allocate(nulibtable_inv_energies(nulibtable_number_groups))
  allocate(nulibtable_ewidths(nulibtable_number_groups))
  allocate(nulibtable_ebottom(nulibtable_number_groups))
  allocate(nulibtable_etop(nulibtable_number_groups))

  call h5dopen_f(file_id, "nrho",dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nulibtable_nrho, dims1, error)
  call h5dclose_f(dset_id, error)
  cerror = cerror + error

  call h5dopen_f(file_id, "ntemp",dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nulibtable_ntemp, dims1, error)
  call h5dclose_f(dset_id, error)
  cerror = cerror + error

  call h5dopen_f(file_id, "nye",dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nulibtable_nye, dims1, error)
  call h5dclose_f(dset_id, error)
  cerror = cerror + error

  !lets also read neutrino energies, bin bottoms,tops and widths
  rank = 1
  dims1(1) = nulibtable_number_groups
  call h5dopen_f(file_id, "neutrino_energies",dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,nulibtable_energies, dims1, error)
  call h5dclose_f(dset_id, error)
  cerror = cerror + error      
  
  !calculate inverse energies
  nulibtable_inv_energies = 1.0d0/nulibtable_energies


  call h5dopen_f(file_id, "bin_widths", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,nulibtable_ewidths, dims1, error)
  call h5dclose_f(dset_id, error)
  cerror = cerror + error   

  call h5dopen_f(file_id, "bin_bottom", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,nulibtable_ebottom, dims1, error)
  call h5dclose_f(dset_id, error)
  cerror = cerror + error   

  call h5dopen_f(file_id, "bin_top", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,nulibtable_etop, dims1, error)
  call h5dclose_f(dset_id, error)
  cerror = cerror + error   

  allocate(nulibtable_logrho(nulibtable_nrho))
  allocate(nulibtable_logtemp(nulibtable_ntemp))
  allocate(nulibtable_ye(nulibtable_nye))
    
  rank = 1
  dims1(1) = nulibtable_nrho
  call h5dopen_f(file_id, "rho_points", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,nulibtable_logrho, dims1, error)
  call h5dclose_f(dset_id, error)
  cerror = cerror + error     

  rank = 1
  dims1(1) = nulibtable_ntemp
  call h5dopen_f(file_id, "temp_points", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,nulibtable_logtemp, dims1, error)
  call h5dclose_f(dset_id, error)
  cerror = cerror + error  

  rank = 1
  dims1(1) = nulibtable_nye
  call h5dopen_f(file_id, "ye_points", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,nulibtable_ye, dims1, error)
  call h5dclose_f(dset_id, error)
  cerror = cerror + error   

  nulibtable_logrho = log10(nulibtable_logrho)
  nulibtable_logrho_min = nulibtable_logrho(1)
  nulibtable_logrho_max = nulibtable_logrho(nulibtable_nrho)
  nulibtable_logtemp = log10(nulibtable_logtemp) 
  nulibtable_logtemp_min = nulibtable_logtemp(1)
  nulibtable_logtemp_max = nulibtable_logtemp(nulibtable_ntemp)
  nulibtable_ye_min = nulibtable_ye(1)
  nulibtable_ye_max = nulibtable_ye(nulibtable_nye)
  
  !now read three tables
  rank = 5
  dims5(1) = nulibtable_nrho
  dims5(2) = nulibtable_ntemp
  dims5(3) = nulibtable_nye
  dims5(4) = nulibtable_number_species  
  dims5(5) = nulibtable_number_groups  


  allocate(nulibtable_temp(nulibtable_nrho,nulibtable_ntemp, &
       nulibtable_nye,nulibtable_number_species,nulibtable_number_groups))  
  allocate(nulibtable_emissivities(nulibtable_nrho,nulibtable_ntemp, &
       nulibtable_nye,nulibtable_number_species*nulibtable_number_groups))
  allocate(nulibtable_absopacity(nulibtable_nrho,nulibtable_ntemp, &
       nulibtable_nye,nulibtable_number_species*nulibtable_number_groups))
  allocate(nulibtable_scatopacity(nulibtable_nrho,nulibtable_ntemp, &
       nulibtable_nye,nulibtable_number_species*nulibtable_number_groups))

  call h5dopen_f(file_id, "emissivities", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,nulibtable_temp, dims5, error)
  call h5dclose_f(dset_id, error)
  cerror = cerror + error   

  do i=1,nulibtable_number_species
     startindex = (i-1)*nulibtable_number_groups+1
     endindex = startindex + nulibtable_number_groups - 1
     do j = 1,nulibtable_number_groups
        nulibtable_emissivities(:,:,:,startindex+j-1) = log10(max(1.0d-30,&
             nulibtable_temp(:,:,:,i,j)*nulibtable_ewidths(j)))
     enddo
  enddo

  call h5dopen_f(file_id, "absorption_opacity",dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,nulibtable_temp, dims5, error)
  call h5dclose_f(dset_id, error)
  cerror = cerror + error   

  do i=1,nulibtable_number_species
     startindex = (i-1)*nulibtable_number_groups+1
     endindex = startindex + nulibtable_number_groups - 1
     nulibtable_absopacity(:,:,:,startindex:endindex) = log10(nulibtable_temp(:,:,:,i,:))
  enddo

  call h5dopen_f(file_id, "scattering_opacity", dset_id, error)
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,nulibtable_temp, dims5, error)
  call h5dclose_f(dset_id, error)
  cerror = cerror + error   

  do i=1,nulibtable_number_species
     startindex = (i-1)*nulibtable_number_groups+1
     endindex = startindex + nulibtable_number_groups - 1
     nulibtable_scatopacity(:,:,:,startindex:endindex) = log10(nulibtable_temp(:,:,:,i,:))
  enddo

  nulibtable_number_easvariables = 3

  !must close h5 files, check for error
  if (cerror.ne.0) then
     write (*,*) "We have errors on reading HDF5 file", cerror
     stop
  endif  

  call h5fclose_f(file_id,error)
  call h5close_f(error)

end subroutine nulibtable_reader
