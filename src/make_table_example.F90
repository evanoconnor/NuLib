!-*-f90-*-
program make_table_example
  
  use nulib
  use inputparser
#if NUCLEI_HEMPEL
  use nuclei_hempel
#endif
  implicit none
#ifdef __MPI__
  include 'mpif.h'
#endif
  !many people use different number of species, this is to denote how they are devided up.
  ! mytable_neutrino_scheme = 1 (three output species)
  ! species #1: electron neutrino             #2 electron antineutrino
  !         #3: muon+tau neutrino+antineutrino
  ! neutrino_scheme = 2 (four output species)
  ! species #1: electron neutrino             #2 electron antineutrino
  !         #3: muon+tau neutrino             #4 mu and tau antineutrino
  ! neutrino_scheme = 3 (six output species)
  ! species #1: electron neutrino             #2 electron antineutrino
  !         #3: muon neutrino                 #4 mu antineutrino
  !         #5: tau neutrino                  #6 tau antineutrino
  integer :: mytable_neutrino_scheme = 1

  !number of species for nulib to calculation interactions for, must
  !be six currently, average is done via above parameter
  integer :: mytable_number_species = 6

  !number of species to output, must match comments above
  integer :: number_output_species = 3

  !number of energy groups
  integer :: mytable_number_groups = 18

  !NuLib parameters file (weak rates and EOS)
  character*200 :: parameters_filename = "./parameters"


  !final table parameters
  integer :: final_table_size_ye, final_table_size_rho, final_table_size_temp
  real*8  :: min_logrho,max_logrho
  real*8  :: min_logtemp,max_logtemp
  real*8  :: min_ye,max_ye

  character*512 :: finaltable_filename
  real*8, allocatable,dimension(:) :: table_rho
  real*8, allocatable,dimension(:) :: table_temp
  real*8, allocatable,dimension(:) :: table_ye
  real*8, allocatable,dimension(:,:,:,:,:) :: table_emission 
  real*8, allocatable,dimension(:,:,:,:,:) :: table_absopacity
  real*8, allocatable,dimension(:,:,:,:,:) :: table_scatopacity 
  real*8, allocatable,dimension(:,:,:,:,:) :: table_delta 


  !final Itable parameters
  integer :: final_Itable_size_temp, final_Itable_size_eta, final_Itable_size_inE
  real*8  :: Imin_logtemp,Imax_logtemp
  real*8  :: Imin_logeta,Imax_logeta
  real*8, allocatable,dimension(:) :: Itable_temp
  real*8, allocatable,dimension(:) :: Itable_eta
  real*8, allocatable,dimension(:) :: Itable_inE
  real*8, allocatable,dimension(:,:,:,:,:) :: Itable_Phi0
  real*8, allocatable,dimension(:,:,:,:,:) :: Itable_Phi1   
  real*8, allocatable,dimension(:,:,:,:,:,:) :: epannihiltable_Phi0 !for ep-annihilation kernels, need both production and destruction kernels
  real*8, allocatable,dimension(:,:,:,:,:,:) :: epannihiltable_Phi1 !for ep-annihilation kernels, need both production and destruction kernels


  !versioning
  real*8 :: timestamp
  character(8) :: date
  integer :: values(8)
  character(100) :: outdir,base,vnum,srho,stemp,sye,sng,sns,sItemp,sIeta

  !local variables to help in making tables
  integer :: irho,itemp,iye,ns,ng
  integer :: ieta,iinE
  real*8, allocatable,dimension(:,:) :: local_emissivity
  real*8, allocatable,dimension(:,:) :: local_absopacity
  real*8, allocatable,dimension(:,:) :: local_scatopacity
  real*8, allocatable,dimension(:,:) :: local_delta
  real*8, allocatable,dimension(:,:) :: local_Phi0
  real*8, allocatable,dimension(:,:) :: local_Phi1
  real*8, allocatable,dimension(:,:,:) :: local_Phi0_epannihil 
  real*8, allocatable,dimension(:,:,:) :: local_Phi1_epannihil 
  real*8, allocatable,dimension(:) :: eos_variables
  real*8 :: matter_prs,matter_ent,matter_cs2,matter_dedt,matter_dpderho,matter_dpdrhoe
  integer :: keytemp,keyerr
  real*8 :: precision = 1.0d-10
  integer :: i
  real*8 dxfac,mindx
  logical :: doing_inelastic, doing_epannihil

#ifdef __MPI__  
  !MPI variables
  integer :: mpirank, nprocs, ierror
  integer :: table_size,irho_save,remainder_size,counter,nmin,recvcount
  integer, allocatable, dimension(:) :: sendcounts,displs
  real*8, allocatable,dimension(:) :: table_rho_subset
  integer :: mpi_final_table_size_rho
  real*8 :: mpitime1,mpitime2
  real*8, allocatable,dimension(:,:,:,:,:) :: table_emission_node
  real*8, allocatable,dimension(:,:,:,:,:) :: table_absopacity_node
  real*8, allocatable,dimension(:,:,:,:,:) :: table_scatopacity_node
  real*8, allocatable,dimension(:,:,:,:,:) :: table_delta_node
  real*8, allocatable,dimension(:) :: Itable_temp_subset
  real*8, allocatable,dimension(:,:,:,:,:) :: Itable_Phi0_node
  real*8, allocatable,dimension(:,:,:,:,:) :: Itable_Phi1_node
  integer :: mpi_final_Itable_size_temp

  
  !MPI initialization
  call mpi_init(ierror)
  call mpi_comm_rank(mpi_comm_world, mpirank, ierror)
  call mpi_comm_size(mpi_comm_world, nprocs, ierror)
  mpitime1=mpi_wtime()
#endif

  call input_parser(parameters_filename)
  !this sets up many coefficients and creates the energy grid (one
  !zone + log spacing) see nulib.F90:initialize_nulib
  call initialize_nulib(mytable_neutrino_scheme,mytable_number_species,mytable_number_groups)
  call read_eos_table(eos_filename) !read in EOS table & set reference mass
#if NUCLEI_HEMPEL
  call set_up_Hempel !set's up EOS for nuclear abundances
#endif

  outdir="./"
  base="NuLib"
  vnum="1.0"
  
  adhoc_nux_factor = 0.0d0 !increase for adhoc nux heating (also set
                           !add_nux_absorption_on_n_and_p to true)
  !set up table
  final_table_size_ye = 51
  final_table_size_rho = 82
  final_table_size_temp = 65
  
  final_Itable_size_temp = 65
  final_Itable_size_eta = 61
  final_Itable_size_inE = mytable_number_groups

  min_ye = 0.035d0
  max_ye = 0.55d0
  min_logrho = 6.0d0
  max_logrho = 15.5d0
  min_logtemp = log10(0.05d0)
  max_logtemp = log10(150.0d0)
  Imin_logtemp = log10(0.05d0)
  Imax_logtemp = log10(150.0d0)
  Imin_logeta = log10(0.1d0)
  Imax_logeta = log10(100.0d0)

  !set up energies bins
  do_integrated_BB_and_emissivity = .false.
  mindx = 2.0d0
  bin_bottom(1) = 0.0d0 !MeV
  bin_bottom(2) = 2.0d0 !MeV
  bin_bottom(3) = bin_bottom(2)+mindx
  bin_bottom(number_groups) = 250.0d0

#ifdef __MPI__
  !set up mpi arrays
  counter = 0
  table_size = final_table_size_rho
  remainder_size = mod(table_size,nprocs)
  nmin = int(table_size/nprocs)
  allocate(table_rho_subset(nmin+1))
  allocate(sendcounts(0:nprocs-1))
  allocate(displs(0:nprocs-1))
  table_rho_subset = 0.0d0
  sendcounts = 0
  displs = 0

  !calculate dimensions for mpi_scatterv
  do i=0,nprocs-1
     if(i.lt.remainder_size)then
        sendcounts(i) = nmin + 1
     else
        sendcounts(i) = nmin
     end if
     displs(i) = counter
     counter = counter + sendcounts(i)
     if(i.eq.mpirank)then
        recvcount = sendcounts(mpirank)
     end if
  end do
#endif
  
  call nulib_series2(number_groups-1,bin_bottom(2),bin_bottom(number_groups),mindx,dxfac)
  do i=4,number_groups
     bin_bottom(i) = bin_bottom(i-1)+(bin_bottom(i-1)-bin_bottom(i-2))*dxfac
  enddo
  
  !calculate bin widths & energies from the bottom of the bin & energy at top on bin
  do i=1,number_groups-1
     energies(i) = (bin_bottom(i)+bin_bottom(i+1))/2.0d0
     bin_widths(i) = bin_bottom(i+1)-bin_bottom(i)
     bin_top(i) = bin_bottom(i+1)
  enddo
  energies(number_groups) = bin_bottom(number_groups)+bin_widths(number_groups-1)*dxfac/2.0d0
  bin_widths(number_groups) = 2.0*(energies(number_groups)-bin_bottom(number_groups))
  bin_top(number_groups) = bin_bottom(number_groups)+bin_widths(number_groups)


  allocate(table_ye(final_table_size_ye))
  allocate(table_rho(final_table_size_rho))
  allocate(table_temp(final_table_size_temp))

  allocate(table_emission(final_table_size_rho,final_table_size_temp, &
       final_table_size_ye,number_output_species,mytable_number_groups))
  allocate(table_absopacity(final_table_size_rho,final_table_size_temp, &
       final_table_size_ye,number_output_species,mytable_number_groups))
  allocate(table_scatopacity(final_table_size_rho,final_table_size_temp, &
       final_table_size_ye,number_output_species,mytable_number_groups))
  allocate(table_delta(final_table_size_rho,final_table_size_temp, &
       final_table_size_ye,number_output_species,mytable_number_groups))

#ifdef __MPI__
  !mpi node tables
  allocate(table_emission_node(final_table_size_rho,final_table_size_temp, &
       final_table_size_ye,number_output_species,mytable_number_groups))
  allocate(table_absopacity_node(final_table_size_rho,final_table_size_temp, &
       final_table_size_ye,number_output_species,mytable_number_groups))
  allocate(table_scatopacity_node(final_table_size_rho,final_table_size_temp, &
       final_table_size_ye,number_output_species,mytable_number_groups))
  allocate(table_delta_node(final_table_size_rho,final_table_size_temp, &
       final_table_size_ye,number_output_species,mytable_number_groups))

  table_emission_node = 0.0d0
#endif
  
  table_emission = 0.0d0

  do iye=1,final_table_size_ye
     table_ye(iye) = min_ye+dble(iye-1)/dble(final_table_size_ye-1)*(max_ye-min_ye)
  enddo

  do irho=1,final_table_size_rho
     table_rho(irho) =  &
          10.0d0**(min_logrho+dble(irho-1)/dble(final_table_size_rho-1)*(max_logrho-min_logrho))
  enddo

  do itemp=1,final_table_size_temp
     table_temp(itemp) = &
          10.0d0**(min_logtemp+dble(itemp-1)/dble(final_table_size_temp-1)*(max_logtemp-min_logtemp))
  enddo

#ifdef __MPI__
  !mpi_scatterv sends portions of table_rho to different nodes
  call mpi_scatterv(table_rho,sendcounts,displs,mpi_double,table_rho_subset,&
       recvcount,mpi_double,0,mpi_comm_world,ierror)
  mpi_final_table_size_rho = recvcount
#endif

  !$OMP PARALLEL DO PRIVATE(itemp,iye,local_emissivity,local_absopacity,local_scatopacity, &
  !$OMP ns,ng,eos_variables,keytemp,keyerr,matter_prs,matter_ent,matter_cs2,matter_dedt, &
  !$OMP matter_dpderho,matter_dpdrhoe,local_delta) 
#ifdef __MPI__
  do irho=1,mpi_final_table_size_rho
#else
  do irho=1,final_table_size_rho
#endif
     !must do declarations here for openmp
     allocate(local_emissivity(number_output_species,mytable_number_groups))
     allocate(local_absopacity(number_output_species,mytable_number_groups))
     allocate(local_scatopacity(number_output_species,mytable_number_groups))
     allocate(local_delta(number_output_species,mytable_number_groups))
     allocate(eos_variables(total_eos_variables))
#ifdef __MPI__
     write(*,*) "Rho:", 100.0*dble(displs(mpirank)+irho-1)/dble(final_table_size_rho),"%"
#else
     write(*,*) "Rho:", 100.0*dble(irho-1)/dble(final_table_size_rho),"%"
#endif
     do itemp=1,final_table_size_temp
        do iye=1,final_table_size_ye

           eos_variables = 0.0d0
#ifdef __MPI__           
           eos_variables(rhoindex) = table_rho_subset(irho)
#else
           eos_variables(rhoindex) = table_rho(irho)
#endif
           eos_variables(tempindex) = table_temp(itemp)
           eos_variables(yeindex) = table_ye(iye)

           !! EOS stuff
           call set_eos_variables(eos_variables)

           !calculate the rho,temp,ye
           call single_point_return_all(eos_variables, &
                local_emissivity,local_absopacity,local_scatopacity, local_delta, &
                mytable_neutrino_scheme)
           
           !check that the number is not NaN or Inf (.gt.1.0d300)
           do ns=1,number_output_species
              do ng=1,mytable_number_groups
                 if (local_emissivity(ns,ng).ne.local_emissivity(ns,ng)) then
                    write(*,"(a,1P3E18.9,i6,i6)") "We have a NaN in emissivity", &
                         eos_variables(rhoindex),eos_variables(tempindex),eos_variables(yeindex),ns,ng
                    stop
                 endif
                 if (local_absopacity(ns,ng).ne.local_absopacity(ns,ng)) then
                    write(*,"(a,1P3E18.9,i6,i6)") "We have a NaN in abs", &
                         eos_variables(rhoindex),eos_variables(tempindex),eos_variables(yeindex),ns,ng
                    stop
                 endif
                 if (local_scatopacity(ns,ng).ne.local_scatopacity(ns,ng)) then
                    write(*,"(a,1P3E18.9,i6,i6)") "We have a NaN in scat", &
                         eos_variables(rhoindex),eos_variables(tempindex),eos_variables(yeindex),ns,ng
                    stop
                 endif
                 if (.not. do_transport_opacities) then
                    if (local_delta(ns,ng).ne.local_delta(ns,ng)) then
                        write(*,"(a,1P3E18.9,i6,i6)") "We have a NaN in scat delta", &
                            eos_variables(rhoindex),eos_variables(tempindex),eos_variables(yeindex),ns,ng
                        stop
                    endif
                 endif
                 
                 if (log10(local_emissivity(ns,ng)).ge.300.0d0) then
                    write(*,"(a,1P4E18.9,i6,i6)") "We have a Inf in emissivity", &
                         local_emissivity(ns,ng),eos_variables(rhoindex), &
                         eos_variables(tempindex),eos_variables(yeindex),ns,ng
                    stop
                 endif
                 if (log10(local_absopacity(ns,ng)).ge.300.0d0) then
                    write(*,"(a,1P4E18.9,i6,i6)") "We have a Inf in abs", &
                         local_absopacity(ns,ng),eos_variables(rhoindex), &
                         eos_variables(tempindex),eos_variables(yeindex),ns,ng
                    stop
                 endif
                 if (log10(local_scatopacity(ns,ng)).ge.300.0d0) then
                    write(*,"(a,1P4E18.9,i6,i6)") "We have a Inf in scat", &
                         local_scatopacity(ns,ng),eos_variables(rhoindex), &
                         eos_variables(tempindex),eos_variables(yeindex),ns,ng
                    stop
                 endif
                 if (.not. do_transport_opacities) then
                    if (local_delta(ns,ng).gt.1.0d0) then
                        write(*,"(a,1P4E18.9,i6,i6)") "delta > 1", &
                            local_delta(ns,ng),eos_variables(rhoindex), &
                            eos_variables(tempindex),eos_variables(yeindex),ns,ng
                        stop
                    endif
                    if (local_delta(ns,ng).lt.-1.0d0) then
                        write(*,"(a,1P4E18.9,i6,i6)") "delta < -1", &
                            local_delta(ns,ng),eos_variables(rhoindex), &
                            eos_variables(tempindex),eos_variables(yeindex),ns,ng
                        stop
                    endif
                 endif
              enddo !do ng=1,mytable_number_groups
           enddo !do ns=1,number_output_species

           !set global table
           do ns=1,number_output_species
              do ng=1,mytable_number_groups
#ifdef __MPI__
                 table_emission_node(displs(mpirank)+irho,itemp,iye,ns,ng) = local_emissivity(ns,ng) !ergs/cm^3/s/MeV/srad
                 table_absopacity_node(displs(mpirank)+irho,itemp,iye,ns,ng) = local_absopacity(ns,ng) !cm^-1
                 table_scatopacity_node(displs(mpirank)+irho,itemp,iye,ns,ng) = local_scatopacity(ns,ng) !cm^-1
                 table_delta_node(displs(mpirank)+irho,itemp,iye,ns,ng) = local_delta(ns,ng) !dimensionless
#else
                 table_emission(irho,itemp,iye,ns,ng) = local_emissivity(ns,ng) !ergs/cm^3/s/MeV/srad
                 table_absopacity(irho,itemp,iye,ns,ng) = local_absopacity(ns,ng) !cm^-1
                 table_scatopacity(irho,itemp,iye,ns,ng) = local_scatopacity(ns,ng) !cm^-1
                 table_delta(irho,itemp,iye,ns,ng) = local_delta(ns,ng) !dimensionless
#endif
              enddo !do ns=1,number_output_species
           enddo !do ng=1,mytable_number_groups
           
        enddo!do iye=1,final_table_size_ye
     enddo!do itemp=1,final_table_size_temp

     deallocate(local_emissivity)
     deallocate(local_absopacity)
     deallocate(local_scatopacity)
     deallocate(local_delta)
     deallocate(eos_variables)
  enddo!do irho=1,final_table_size_rho
  !$OMP END PARALLEL DO! end do

#ifdef __MPI__
  call mpi_barrier(mpi_comm_world, ierror)
  if(mpirank.eq.0)write(*,*) "Finished Opacity Table" 

  table_size = final_table_size_rho*final_table_size_temp* &
       final_table_size_ye*number_output_species*mytable_number_groups
  call mpi_reduce(table_emission_node,table_emission,&
       table_size,mpi_double,mpi_sum,0,mpi_comm_world,ierror)
  call mpi_reduce(table_absopacity_node,table_absopacity,&
       table_size,mpi_double,mpi_sum,0,mpi_comm_world,ierror)
  call mpi_reduce(table_scatopacity_node,table_scatopacity,&
       table_size,mpi_double,mpi_sum,0,mpi_comm_world,ierror)
  call mpi_reduce(table_delta_node,table_delta,&
       table_size,mpi_double,mpi_sum,0,mpi_comm_world,ierror)



  !begin inelastic, timing for mpi purposes
  call mpi_barrier(mpi_comm_world, ierror)
  mpitime2 = mpi_wtime()
  if(mpirank.eq.0)write(*,*) "Total time = ",mpitime2-mpitime1
  mpitime1 = mpi_wtime()

#else
  write(*,*) "Finished Opacity Table" 
#endif

  !now generate the inelastic electron scattering table.  This is
  !stored as a function of T,matter_eta, and ingoing electron energy
  !instead of the standard rho,t,ye.  This is to save space.  Note, we
  !store the zeroth and first moment of the scattering kernal for each
  !point, outgoing energy, and neutrino species.  So far we neglect
  !inelastic positron scattering.
  if (add_nue_Iscattering_electrons.or.add_anue_Iscattering_electrons.or. &
       add_numu_Iscattering_electrons.or.add_anumu_Iscattering_electrons.or. &
       add_nutau_Iscattering_electrons.or.add_anutau_Iscattering_electrons) then

     doing_inelastic = .true.
  else
     doing_inelastic = .false.
  endif

  if (add_nue_kernel_epannihil.or.add_anue_kernel_epannihil.or. &
       add_numu_kernel_epannihil.or.add_anumu_kernel_epannihil.or. &
       add_nutau_kernel_epannihil.or.add_anutau_kernel_epannihil) then

     doing_epannihil = .true.
  else
     doing_epannihil = .false.
  endif

  if (doing_inelastic.or.doing_epannihil) then

     write(*,*) "Making Inelastic Table Opacity Table" 

     if (final_Itable_size_inE.ne.mytable_number_groups) then
        stop "make_table_example: inelastic scattering table not square in energy &
             we assume the same energy sturcture for inelastic scattering"
     endif

#ifdef __MPI__
     !set up mpi variables for inelastic scattering
     deallocate(sendcounts)
     deallocate(displs)
     counter = 0
     table_size = final_Itable_size_temp
     remainder_size = mod(table_size,nprocs)
     nmin = int(table_size/nprocs)
     allocate(Itable_temp_subset(nmin+1))
     allocate(sendcounts(0:nprocs-1))
     allocate(displs(0:nprocs-1))
     Itable_temp_subset = 0.0d0
     sendcounts = 0
     recvcount = 0
     displs = 0

     !calculate dimensions for mpi_scatterv
     do i=0,nprocs-1
        if(i.lt.remainder_size)then
           sendcounts(i) = nmin + 1
        else
           sendcounts(i) = nmin
        end if
        displs(i) = counter
        counter = counter + sendcounts(i)
        if(i.eq.mpirank)then
           recvcount = sendcounts(mpirank)
        end if
     end do
#endif

     allocate(Itable_temp(final_Itable_size_temp))
     allocate(Itable_eta(final_Itable_size_eta))
     allocate(Itable_inE(final_Itable_size_inE))

     allocate(Itable_Phi0(final_Itable_size_temp,final_Itable_size_eta, &
          final_Itable_size_inE,number_output_species,mytable_number_groups))
     allocate(Itable_Phi1(final_Itable_size_temp,final_Itable_size_eta, &
          final_Itable_size_inE,number_output_species,mytable_number_groups))
     allocate(epannihiltable_Phi0(final_Itable_size_temp,final_Itable_size_eta, &
          final_Itable_size_inE,number_output_species,mytable_number_groups,2))
     allocate(epannihiltable_Phi1(final_Itable_size_temp,final_Itable_size_eta, &
          final_Itable_size_inE,number_output_species,mytable_number_groups,2))

#ifdef __MPI__
     !mpi node tables for inelastic
     allocate(Itable_Phi0_node(final_Itable_size_temp,final_Itable_size_eta, &
          final_Itable_size_inE,number_output_species,mytable_number_groups))
     allocate(Itable_Phi1_node(final_Itable_size_temp,final_Itable_size_eta, &
          final_Itable_size_inE,number_output_species,mytable_number_groups))
#endif


     do itemp=1,final_Itable_size_temp
        Itable_temp(itemp) = &
             10.0d0**(Imin_logtemp+dble(itemp-1)/dble(final_Itable_size_temp-1)*(Imax_logtemp-Imin_logtemp))
     enddo
     do ieta=1,final_Itable_size_eta
        Itable_eta(ieta) = &
             10.0d0**(Imin_logeta+dble(ieta-1)/dble(final_Itable_size_eta-1)*(Imax_logeta-Imin_logeta))
     enddo

#ifdef __MPI__
     !mpi_scatterv sends portions of Itable_temp to different nodes
     call mpi_scatterv(Itable_temp,sendcounts,displs,mpi_double,Itable_temp_subset,&
          recvcount,mpi_double,0,mpi_comm_world,ierror)
     mpi_final_Itable_size_temp = recvcount
#endif
     !$OMP PARALLEL DO PRIVATE(local_Phi0,local_Phi1,local_Phi0_epannihil,local_Phi1_epannihil,ieta,iinE,ns,ng)
     !loop over temp,eta,inE of table, do each point
#ifdef __MPI__
     do itemp=1,mpi_final_Itable_size_temp
#else
     do itemp=1,final_Itable_size_temp
#endif
        !must do declarations here for openmp
        allocate(local_Phi0(number_output_species,mytable_number_groups))
        allocate(local_Phi1(number_output_species,mytable_number_groups))
        allocate(local_Phi0_epannihil(number_output_species,mytable_number_groups,2))
        allocate(local_Phi1_epannihil(number_output_species,mytable_number_groups,2))

#ifdef __MPI__
        write(*,*) "Temp:", 100.0*dble(displs(mpirank)+itemp-1)/dble(final_Itable_size_temp),"%"
#else
        write(*,*) "Temp:", 100.0*dble(itemp-1)/dble(final_Itable_size_temp),"%"
#endif

        do ieta=1,final_Itable_size_eta
           do iinE=final_Itable_size_inE,1,-1

#ifdef __MPI__
              call single_Ipoint_return_all(iinE,Itable_eta(ieta), &
                   Itable_temp_subset(itemp),local_Phi0,local_Phi1,mytable_neutrino_scheme)
#else
              call single_Ipoint_return_all(iinE,Itable_eta(ieta), &
                   Itable_temp(itemp),local_Phi0,local_Phi1,mytable_neutrino_scheme)
#endif
              
              !fill in higher out energies with partner
              !Rout(iinE,E>iinE) is not calculated
              !set equal to Rin(E>iinE,iinE) which was calculated already

              !note, experience has taught us that if you interpolate these
              !kernels you must reapply this symmetry after interpolation, we
              !put them here for completeness.

              do ns=1,number_output_species
                 do ng=iinE+1,final_Itable_size_inE
#ifdef __MPI__
                    local_Phi0(ns,ng) = exp(-(energies(ng)-energies(iinE))/Itable_temp_subset(itemp))* &
                         Itable_Phi0_node(displs(mpirank)+itemp,ieta,ng,ns,iinE)
                    local_Phi1(ns,ng) = exp(-(energies(ng)-energies(iinE))/Itable_temp_subset(itemp))* &
                         Itable_Phi1_node(displs(mpirank)+itemp,ieta,ng,ns,iinE)
#else
                    local_Phi0(ns,ng) = exp(-(energies(ng)-energies(iinE))/Itable_temp(itemp))* &
                         Itable_Phi0(itemp,ieta,ng,ns,iinE)
                    local_Phi1(ns,ng) = exp(-(energies(ng)-energies(iinE))/Itable_temp(itemp))* &
                         Itable_Phi1(itemp,ieta,ng,ns,iinE)
#endif
                 enddo
              enddo

              !calculate and check that the number is not NaN or Inf
              !(.gt.1.0d300)
              do ns=1,number_output_species
                 do ng=1,mytable_number_groups

#ifdef __MPI__
                    if (local_Phi0(ns,ng).ne.local_Phi0(ns,ng)) then
                       write(*,"(a,1P2E18.9,i6,i6,i6)") "We have a NaN in Phi0", &
                            Itable_temp_subset(itemp),Itable_eta(ieta),iinE,ns,ng
                       stop
                    endif
                    if (local_Phi1(ns,ng).ne.local_Phi1(ns,ng)) then
                       write(*,"(a,1P2E18.9,i6,i6,i6)") "We have a NaN in Phi1", &
                            Itable_temp_subset(itemp),Itable_eta(ieta),iinE,ns,ng
                       stop
                    endif

                    if (log10(local_Phi0(ns,ng)).ge.300.0d0) then
                       write(*,"(a,1P3E18.9,i6,i6,i6)") "We have a Inf in Phi0", &
                            local_Phi0(ns,ng),Itable_temp_subset(itemp),Itable_eta(ieta),iinE,ns,ng
                       stop
                    endif
                    if (log10(local_Phi1(ns,ng)).ge.300.0d0) then
                       write(*,"(a,1P3E18.9,i6,i6,i6)") "We have a Inf in Phi1", &
                            local_Phi1(ns,ng),Itable_temp_subset(itemp),Itable_eta(ieta),iinE,ns,ng
                       stop
                    endif
#else
                    if (local_Phi0(ns,ng).ne.local_Phi0(ns,ng)) then
                       write(*,"(a,1P2E18.9,i6,i6,i6)") "We have a NaN in Phi0", &
                            Itable_temp(itemp),Itable_eta(ieta),iinE,ns,ng
                       stop
                    endif
                    if (local_Phi1(ns,ng).ne.local_Phi1(ns,ng)) then
                       write(*,"(a,1P2E18.9,i6,i6,i6)") "We have a NaN in Phi1", &
                            Itable_temp(itemp),Itable_eta(ieta),iinE,ns,ng
                       stop
                    endif

                    if (log10(local_Phi0(ns,ng)).ge.300.0d0) then
                       write(*,"(a,1P3E18.9,i6,i6,i6)") "We have a Inf in Phi0", &
                            local_Phi0(ns,ng),Itable_temp(itemp),Itable_eta(ieta),iinE,ns,ng
                       stop
                    endif
                    if (log10(local_Phi1(ns,ng)).ge.300.0d0) then
                       write(*,"(a,1P3E18.9,i6,i6,i6)") "We have a Inf in Phi1", &
                            local_Phi1(ns,ng),Itable_temp(itemp),Itable_eta(ieta),iinE,ns,ng
                       stop
                    endif
#endif

                 enddo !do ng=1,mytable_number_groups
              enddo !do ns=1,number_output_species

              !set global table
              do ns=1,number_output_species
                 do ng=1,mytable_number_groups
#ifdef __MPI__
                    Itable_Phi0_node(displs(mpirank)+itemp,ieta,iinE,ns,ng) = local_Phi0(ns,ng) !cm^3/s
                    Itable_Phi1_node(displs(mpirank)+itemp,ieta,iinE,ns,ng) = local_Phi1(ns,ng) !cm^3/s
#else
                    Itable_Phi0(itemp,ieta,iinE,ns,ng) = local_Phi0(ns,ng) !cm^3/s
                    Itable_Phi1(itemp,ieta,iinE,ns,ng) = local_Phi1(ns,ng) !cm^3/s
#endif
                 enddo !do ns=1,number_output_species
              enddo !do ng=1,mytable_number_groups
           
              call single_epannihil_kernel_point_return_all(iinE,Itable_eta(ieta), &
                   Itable_temp(itemp),local_Phi0_epannihil,local_Phi1_epannihil,mytable_neutrino_scheme)               

              !calculate and check that the number is not NaN or Inf
              !(.gt.1.0d300)
              do ns=1,number_output_species
                 do ng=1,mytable_number_groups

                    if (local_Phi0_epannihil(ns,ng,1).ne.local_Phi0_epannihil(ns,ng,1)) then
                       write(*,"(a,1P2E18.9,i6,i6,i6)") "We have a NaN in Phi0_epannihil production", &
                            Itable_temp(itemp),Itable_eta(ieta),iinE,ns,ng
                       stop
                    endif
                    if (local_Phi0_epannihil(ns,ng,2).ne.local_Phi0_epannihil(ns,ng,2)) then
                       write(*,"(a,1P2E18.9,i6,i6,i6)") "We have a NaN in Phi0_epannihil annihilation", &
                            Itable_temp(itemp),Itable_eta(ieta),iinE,ns,ng
                       stop
                    endif
                    if (local_Phi1_epannihil(ns,ng,1).ne.local_Phi1_epannihil(ns,ng,1)) then
                       write(*,"(a,1P2E18.9,i6,i6,i6)") "We have a NaN in Phi1_epannihil production", &
                            Itable_temp(itemp),Itable_eta(ieta),iinE,ns,ng
                       stop
                    endif
                    if (local_Phi1_epannihil(ns,ng,2).ne.local_Phi1_epannihil(ns,ng,2)) then
                       write(*,"(a,1P2E18.9,i6,i6,i6)") "We have a NaN in Phi1_epannihil annihilation", &
                            Itable_temp(itemp),Itable_eta(ieta),iinE,ns,ng
                       stop
                    endif
                    
                    if (log10(local_Phi0_epannihil(ns,ng,1)).ge.300.0d0) then
                       write(*,"(a,1P3E18.9,i6,i6,i6)") "We have a Inf in Phi0_epannihil production", &
                            local_Phi0_epannihil(ns,ng,1),Itable_temp(itemp),Itable_eta(ieta),iinE,ns,ng
                       stop
                    endif
                    if (log10(local_Phi0_epannihil(ns,ng,2)).ge.300.0d0) then
                       write(*,"(a,1P3E18.9,i6,i6,i6)") "We have a Inf in Phi0_epannihil annihilation", &
                            local_Phi0_epannihil(ns,ng,2),Itable_temp(itemp),Itable_eta(ieta),iinE,ns,ng
                       stop
                    endif
                    if (log10(local_Phi1_epannihil(ns,ng,1)).ge.300.0d0) then
                       write(*,"(a,1P3E18.9,i6,i6,i6)") "We have a Inf in Phi1_epannihil production", &
                            local_Phi1_epannihil(ns,ng,1),Itable_temp(itemp),Itable_eta(ieta),iinE,ns,ng
                       stop
                    endif
                    if (log10(local_Phi1_epannihil(ns,ng,2)).ge.300.0d0) then
                       write(*,"(a,1P3E18.9,i6,i6,i6)") "We have a Inf in Phi1_epannihil annihilation", &
                            local_Phi1_epannihil(ns,ng,2),Itable_temp(itemp),Itable_eta(ieta),iinE,ns,ng
                       stop
                    endif

                 enddo !do ng=1,mytable_number_groups
              enddo !do ns=1,number_output_species

              !set global table
              do ns=1,number_output_species
                 do ng=1,mytable_number_groups
                    epannihiltable_Phi0(itemp,ieta,iinE,ns,ng,1) = local_Phi0_epannihil(ns,ng,1) !cm^3/s
                    epannihiltable_Phi0(itemp,ieta,iinE,ns,ng,2) = local_Phi0_epannihil(ns,ng,2) !cm^3/s
                    epannihiltable_Phi1(itemp,ieta,iinE,ns,ng,1) = local_Phi1_epannihil(ns,ng,1) !cm^3/s
                    epannihiltable_Phi1(itemp,ieta,iinE,ns,ng,2) = local_Phi1_epannihil(ns,ng,2) !cm^3/s
                enddo !do ns=1,number_output_species
              enddo !do ng=1,mytable_number_groups              

           enddo!do iinE=1,final_Itable_size_inE
        enddo!do ieta=1,final_Itable_size_eta

        deallocate(local_Phi0)
        deallocate(local_Phi1)
        deallocate(local_Phi0_epannihil)
        deallocate(local_Phi1_epannihil)
     enddo!do itemp=1,final_Itable_size_temp
     !$OMP END PARALLEL DO! end do

#ifdef __MPI__
     call mpi_barrier(mpi_comm_world, ierror)
     if(mpirank.eq.0)write(*,*) "Finished Inelastic Table" 
     table_size = final_Itable_size_temp*final_Itable_size_eta* &
          final_Itable_size_inE*number_output_species*mytable_number_groups
     call mpi_reduce(Itable_Phi0_node,Itable_Phi0,&
          table_size,mpi_double,mpi_sum,0,mpi_comm_world,ierror)
     call mpi_reduce(Itable_Phi1_node,Itable_Phi1,&
          table_size,mpi_double,mpi_sum,0,mpi_comm_world,ierror)
     call mpi_barrier(mpi_comm_world, ierror)
     mpitime2 = mpi_wtime()
     if(mpirank.eq.0)write(*,*) "Total inelastic time = ",mpitime2-mpitime1
#endif

  endif

  !write out table in H5 format
#ifdef __MPI__
  !only the first mpi node should write
  if(mpirank.eq.0)then
#endif
     call date_and_time(DATE=date,VALUES=values)
     write(srho,*) final_table_size_rho
     write(stemp,*) final_table_size_temp
     write(sye,*) final_table_size_ye
     write(sng,*) mytable_number_groups
     write(sns,*) number_output_species
     write(sItemp,*) final_Itable_size_temp
     write(sIeta,*) final_Itable_size_eta
     timestamp = dble(values(1))*10000.0d0+dble(values(2))*100.0+dble(values(3)) + &
          (dble(values(5))+dble(values(6))/60.0d0 + dble(values(7))/3600.0d0 )/24.0

     if (doing_inelastic.or.doing_epannihil) then
        finaltable_filename = trim(adjustl(outdir))//trim(adjustl(base))//"_rho"//trim(adjustl(srho))// &
             "_temp"//trim(adjustl(stemp))//"_ye"//trim(adjustl(sye))// &
             "_ng"//trim(adjustl(sng))//"_ns"//trim(adjustl(sns))// &
             "_Itemp"//trim(adjustl(sItemp))//"_Ieta"//trim(adjustl(sIeta))// &
             "_version"//trim(adjustl(vnum))//"_"//trim(adjustl(date))//".h5"
     else
        finaltable_filename = trim(adjustl(outdir))//trim(adjustl(base))//"_rho"//trim(adjustl(srho))// &
             "_temp"//trim(adjustl(stemp))//"_ye"//trim(adjustl(sye))// &
             "_ng"//trim(adjustl(sng))//"_ns"//trim(adjustl(sns))// &
             "_version"//trim(adjustl(vnum))//"_"//trim(adjustl(date))//".h5"
     endif

     call write_h5(finaltable_filename,timestamp)
#ifdef __MPI__
  end if

  call mpi_barrier(mpi_comm_world, ierror)
  call mpi_finalize(ierror)
#endif

contains

  subroutine write_h5(filename,timestamp)
    
    use nulib
    use hdf5
    implicit none

    character(len=512) :: filename
    
    !H5 stuff
    integer :: error,rank,cerror
    integer(HID_T) :: file_id,dset_id,dspace_id
    integer(HSIZE_T) :: dims1(1), dims2(2), dims3(3), dims4(4), dims5(5), dims6(6)!, etc....
    
    real*8 :: timestamp
    character(8) :: date
    integer :: values(8)
    
    cerror = 0
    
    !open HDF5 file, given filename
    !note: H5F_ACC_TRUNC specifies that if the file already exists,
    !  the current contents will be deleted so that the application can
    !  rewrite the file with new data.  H5F_ACC_EXCL specifies that the open
    !  is to fail if the file already exists.  If the file does not already
    !  exist, the file access parameter is ignored.  
    !file_id gets set here for future use in distinguishing files
    !if error != 0 we will know at the end when we ensure cerror (c=cummulative) == 0
    call h5open_f(error)
    cerror = cerror + error 
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,error) 
    cerror = cerror + error
    
    !write scalars (rank=1, dims1(1) = 1)
    rank = 1
    dims1(1) = 1
  
    !first lets write number of species and number of groups and number of rho,temp, and ye points, also timestamp
    call h5screate_simple_f(rank, dims1, dspace_id, error)
    call h5dcreate_f(file_id, "number_species", H5T_NATIVE_INTEGER, &
         dspace_id,dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, number_output_species, dims1, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)  
    cerror = cerror + error
    
    call h5screate_simple_f(rank, dims1, dspace_id, error)
    call h5dcreate_f(file_id, "timestamp", H5T_NATIVE_DOUBLE, &
         dspace_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, timestamp, &
         dims1, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)
    
    call h5screate_simple_f(rank, dims1, dspace_id, error)
    call h5dcreate_f(file_id, "number_groups", H5T_NATIVE_INTEGER, &
         dspace_id,dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, number_groups, dims1, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)  
    cerror = cerror + error
    
    call h5screate_simple_f(rank, dims1, dspace_id, error)
    call h5dcreate_f(file_id, "nrho", H5T_NATIVE_INTEGER, &
         dspace_id,dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, final_table_size_rho, dims1, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)  
    cerror = cerror + error
    
    call h5screate_simple_f(rank, dims1, dspace_id, error)
    call h5dcreate_f(file_id, "ntemp", H5T_NATIVE_INTEGER, &
         dspace_id,dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, final_table_size_temp, dims1, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)  
    cerror = cerror + error
    
    call h5screate_simple_f(rank, dims1, dspace_id, error)
    call h5dcreate_f(file_id, "nye", H5T_NATIVE_INTEGER, &
         dspace_id,dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, final_table_size_ye, dims1, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)  
    cerror = cerror + error
    
    !lets also write neutrino energies, bin bottoms,tops and widths
    !write 1D arrays (rank=1, dims1(1) = oneDarray_size)
    rank = 1
    dims1(1) = number_groups
    call h5screate_simple_f(rank, dims1, dspace_id, error)
    call h5dcreate_f(file_id, "neutrino_energies", H5T_NATIVE_DOUBLE, &
         dspace_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,energies, dims1, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)  
    cerror = cerror + error      
    
    call h5screate_simple_f(rank, dims1, dspace_id, error)
    call h5dcreate_f(file_id, "bin_widths", H5T_NATIVE_DOUBLE, &
         dspace_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,bin_widths, dims1, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)  
    cerror = cerror + error   
    
    call h5screate_simple_f(rank, dims1, dspace_id, error)
    call h5dcreate_f(file_id, "bin_bottom", H5T_NATIVE_DOUBLE, &
         dspace_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,bin_bottom, dims1, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)  
    cerror = cerror + error   
    
    call h5screate_simple_f(rank, dims1, dspace_id, error)
    call h5dcreate_f(file_id, "bin_top", H5T_NATIVE_DOUBLE, &
         dspace_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,bin_top, dims1, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)  
    cerror = cerror + error   
    
    rank = 1
    dims1(1) = final_table_size_rho
    call h5screate_simple_f(rank, dims1, dspace_id, error)
    call h5dcreate_f(file_id, "rho_points", H5T_NATIVE_DOUBLE, &
         dspace_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,table_rho, dims1, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)  
    cerror = cerror + error     
    
    rank = 1
    dims1(1) = final_table_size_temp
    call h5screate_simple_f(rank, dims1, dspace_id, error)
    call h5dcreate_f(file_id, "temp_points", H5T_NATIVE_DOUBLE, &
         dspace_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,table_temp, dims1, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)  
    cerror = cerror + error   
    
    rank = 1
    dims1(1) = final_table_size_ye
    call h5screate_simple_f(rank, dims1, dspace_id, error)
    call h5dcreate_f(file_id, "ye_points", H5T_NATIVE_DOUBLE, &
         dspace_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,table_ye, dims1, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)  
    cerror = cerror + error   
    
    rank = 5
    dims5(1) = final_table_size_rho
    dims5(2) = final_table_size_temp
    dims5(3) = final_table_size_ye
    dims5(4) = number_output_species  
    dims5(5) = number_groups  
    
    call h5screate_simple_f(rank, dims5, dspace_id, error)
    call h5dcreate_f(file_id, "emissivities", H5T_NATIVE_DOUBLE, &
         dspace_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,table_emission, dims5, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)  
    cerror = cerror + error   
    
    call h5screate_simple_f(rank, dims5, dspace_id, error)
    call h5dcreate_f(file_id, "absorption_opacity", H5T_NATIVE_DOUBLE, &
         dspace_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,table_absopacity, dims5, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)  
    cerror = cerror + error   
    
    call h5screate_simple_f(rank, dims5, dspace_id, error)
    call h5dcreate_f(file_id, "scattering_opacity", H5T_NATIVE_DOUBLE, &
         dspace_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,table_scatopacity, dims5, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)  
    cerror = cerror + error   

    if(.not. do_transport_opacities) then
       call h5screate_simple_f(rank, dims5, dspace_id, error)
       call h5dcreate_f(file_id, "scattering_delta", H5T_NATIVE_DOUBLE, &
            dspace_id, dset_id, error)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,table_delta, dims5, error)
       call h5dclose_f(dset_id, error)
       call h5sclose_f(dspace_id, error)
       cerror = cerror + error
    endif
    
    if (doing_inelastic.or.doing_epannihil) then

       rank = 1
       dims1(1) = 1
       call h5screate_simple_f(rank, dims1, dspace_id, error)
       call h5dcreate_f(file_id, "Itemp", H5T_NATIVE_INTEGER, &
            dspace_id,dset_id, error)
       call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, final_Itable_size_temp, dims1, error)
       call h5dclose_f(dset_id, error)
       call h5sclose_f(dspace_id, error)  
       cerror = cerror + error

       rank = 1
       dims1(1) = 1
       call h5screate_simple_f(rank, dims1, dspace_id, error)
       call h5dcreate_f(file_id, "Ieta", H5T_NATIVE_INTEGER, &
            dspace_id,dset_id, error)
       call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, final_Itable_size_eta, dims1, error)
       call h5dclose_f(dset_id, error)
       call h5sclose_f(dspace_id, error)  
       cerror = cerror + error

       rank = 1
       dims1(1) = final_Itable_size_temp
       call h5screate_simple_f(rank, dims1, dspace_id, error)
       call h5dcreate_f(file_id, "temp_Ipoints", H5T_NATIVE_DOUBLE, &
            dspace_id, dset_id, error)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,Itable_temp, dims1, error)
       call h5dclose_f(dset_id, error)
       call h5sclose_f(dspace_id, error)  
       cerror = cerror + error   

       rank = 1
       dims1(1) = final_Itable_size_eta
       call h5screate_simple_f(rank, dims1, dspace_id, error)
       call h5dcreate_f(file_id, "eta_Ipoints", H5T_NATIVE_DOUBLE, &
            dspace_id, dset_id, error)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,Itable_eta, dims1, error)
       call h5dclose_f(dset_id, error)
       call h5sclose_f(dspace_id, error)  
       cerror = cerror + error   
       
    endif

    if(doing_inelastic) then

       rank = 5
       dims5(1) = final_Itable_size_temp
       dims5(2) = final_Itable_size_eta
       dims5(3) = final_Itable_size_inE
       dims5(4) = number_output_species  
       dims5(5) = number_groups  
       
       call h5screate_simple_f(rank, dims5, dspace_id, error)
       call h5dcreate_f(file_id, "inelastic_phi0", H5T_NATIVE_DOUBLE, &
            dspace_id, dset_id, error)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,Itable_Phi0, dims5, error)
       call h5dclose_f(dset_id, error)
       call h5sclose_f(dspace_id, error)  
       cerror = cerror + error   
       
       call h5screate_simple_f(rank, dims5, dspace_id, error)
       call h5dcreate_f(file_id, "inelastic_phi1", H5T_NATIVE_DOUBLE, &
            dspace_id, dset_id, error)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,Itable_Phi1, dims5, error)
       call h5dclose_f(dset_id, error)
       call h5sclose_f(dspace_id, error)  
       cerror = cerror + error
    endif

    if (doing_epannihil) then
       rank = 6
       dims6(1) = final_Itable_size_temp
       dims6(2) = final_Itable_size_eta
       dims6(3) = final_Itable_size_inE
       dims6(4) = number_output_species  
       dims6(5) = number_groups
       dims6(6) = 2
       
       call h5screate_simple_f(rank, dims6, dspace_id, error)
       call h5dcreate_f(file_id, "epannihil_phi0", H5T_NATIVE_DOUBLE, &
            dspace_id, dset_id, error)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,epannihiltable_Phi0, dims6, error)
       call h5dclose_f(dset_id, error)
       call h5sclose_f(dspace_id, error)  
       cerror = cerror + error   
       
       call h5screate_simple_f(rank, dims6, dspace_id, error)
       call h5dcreate_f(file_id, "epannihil_phi1", H5T_NATIVE_DOUBLE, &
            dspace_id, dset_id, error)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,epannihiltable_Phi1, dims6, error)
       call h5dclose_f(dset_id, error)
       call h5sclose_f(dspace_id, error)  
       cerror = cerror + error

    endif


    !must close h5 files, check for error
    if (cerror.ne.0) then
       write(*,*) "We have errors on writing HDF5 file", cerror
       stop
    endif
    
    call h5fclose_f(file_id,error)
    call h5close_f(error)
    
  end subroutine write_h5

        
end program make_table_example
