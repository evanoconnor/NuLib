!-*-f90-*-
program point_example

  use nulib
  use weak_rates
  implicit none

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
  integer :: mypoint_neutrino_scheme = 1

  !number of species for nulib to calculate interactions for, must
  !be six currently, average is done via above parameter
  integer :: mypoint_number_species = 6

  !number of species for nulib to output interactions for, must
  !be commensurate with neutrino_scheme above
  integer :: mypoint_number_output_species = 3

  !number of energy groups
  integer :: mypoint_number_groups = 24

  !EOS table
  character*200 :: eos_filename = "/user/sullivan/gr1dnulib/GitHub/NuLib/src/extra_code_and_tables/Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5"

  !Weak rate data (currently LMP rates only)
  character*200 :: weakrates_filename = "/projects/ceclub/gr1dnulib/GitHub/NuLib/src/extra_code_and_tables/rates-ext.out"

  !local variables
  real*8, allocatable,dimension(:,:) :: local_emissivity
  real*8, allocatable,dimension(:,:) :: local_absopacity
  real*8, allocatable,dimension(:,:) :: local_scatopacity
  real*8, allocatable,dimension(:,:) :: local_Phi0, local_Phi1
  real*8, allocatable,dimension(:,:) :: blackbody_spectra
  real*8, allocatable,dimension(:) :: eos_variables
  real*8 :: matter_prs,matter_ent,matter_cs2,matter_dedt,matter_dpderho,matter_dpdrhoe
  integer :: keytemp,keyerr
  real*8 :: precision = 1.0d-10
  real*8 :: xrho, xtemp, xye, xtime
  integer i,j,k
  real*8 dxfac,mindx

  real*8, allocatable,dimension(:) :: mass_fractions
  real*8, allocatable,dimension(:) :: number_densities
  real*8, allocatable,dimension(:) :: probability_dist
  real*8, allocatable,dimension(:) :: emissivity_ni56
  real*8, allocatable,dimension(:) :: blackbody_spectrum
  real*8, allocatable,dimension(:,:) :: emissivity
  real*8, allocatable,dimension(:) :: dye
  real*8, allocatable,dimension(:) :: normalized_dye

  logical :: cont
  character*200 :: time_string
  integer :: nindex
  real*8 :: current_time
  real*8 :: current_energy
  integer :: IO

  real*8 :: fermidirac_dimensionless
  real*8, dimension(27,1610) :: dye_bin
  real*8, dimension(1610) :: dyesum
  real*8, dimension(1610) :: dyesum_time
  real*8 :: time
  real*8, dimension(1610) :: dye_hydro
  real*8, dimension(26,1610) :: dye_bin_runningsum
  real*8, dimension(1610) :: dyesum_int

  !allocate the arrays for the point values
  allocate(local_emissivity(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_absopacity(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_scatopacity(mypoint_number_output_species,mypoint_number_groups))
  allocate(blackbody_spectra(mypoint_number_output_species,mypoint_number_groups))

  !this sets up many cooefficients and creates the energy grid (one
  !zone + log spacing) see nulib.F90:initialize_nulib
  call initialize_nulib(mypoint_neutrino_scheme,mypoint_number_species,mypoint_number_groups)

  !read in EOS table & set reference mass
  call readtable(eos_filename)
  m_ref = m_amu !for SFHo_EOS (Hempel)
  ! m_ref = m_n !for LS220

  !read in weak rates table and build interpolant functions
  call readrates(weakrates_filename,table_bounds)

  allocate(mass_fractions(nspecies))
  allocate(number_densities(nspecies))
  allocate(probability_dist(number_groups))
  allocate(emissivity(nspecies,number_groups))
  allocate(blackbody_spectrum(number_groups))
  allocate(emissivity_ni56(number_groups))
  allocate(dye(nspecies))
  allocate(normalized_dye(nspecies))

  !set up energies bins
  do_integrated_BB_and_emissivity = .false.
  mindx = 1.0d0
  bin_bottom(1) = 0.0d0 !MeV
  bin_bottom(2) = 1.0d0 !MeV
  bin_bottom(3) = bin_bottom(2)+mindx
  bin_bottom(number_groups) = 250.0d0
  
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

  
  allocate(eos_variables(total_eos_variables))


  open(1,file="src/extra_code_and_tables/dye.dat",status='old')
  open(2,file="src/extra_code_and_tables/integrated_dyedt.dat")
  open(3,file="src/extra_code_and_tables/dyedt_hydro_c_t.dat",status='old')
  open(4,file="src/extra_code_and_tables/dyesum.dat",status='old')
  open(5,file="src/extra_code_and_tables/ye_from_int.dat")
  cont = .true.
  nindex = 1

  do while(cont)
     write(*,*) nindex
     read(1,*)dye_bin(1,nindex),dye_bin(2,nindex),dye_bin(3,nindex),dye_bin(4,nindex),dye_bin(5,nindex),dye_bin(6,nindex),dye_bin(7,nindex),dye_bin(8,nindex),dye_bin(9,nindex),dye_bin(10,nindex),dye_bin(11,nindex),dye_bin(12,nindex),dye_bin(13,nindex),dye_bin(14,nindex),dye_bin(15,nindex),dye_bin(16,nindex),dye_bin(17,nindex),dye_bin(18,nindex),dye_bin(19,nindex),dye_bin(20,nindex),dye_bin(21,nindex),dye_bin(22,nindex),dye_bin(23,nindex),dye_bin(24,nindex),dye_bin(25,nindex),dye_bin(26,nindex),dye_bin(27,nindex)
     read(3,*)time,dye_hydro(nindex)
     read(4,*)dyesum_time(nindex),dyesum(nindex)
     nindex = nindex + 1
     if (nindex.ge.1610) cont = .false.     
  enddo
  
  do i=2,27
     dye_bin(i,:) = dye_bin(i,:) + dye_hydro(:)/26
  end do
  dyesum(:) = dyesum(:) + dye_hydro(:)
  
  !time integration of binned dyedt
  ! dye_bin_runningsum = 0.0d0
  ! do i=2,26
  !    do j=1,nindex-1
  !       do k=1,j
  !          dye_bin_runningsum(i-1,j) = dye_bin_runningsum(i-1,j) + dye_bin(i,k)*(dye_bin(1,k+1)-dye_bin(1,k))
  !       end do
  !    end do
  ! end do

  !time integration of binned dyedt
  dye_bin_runningsum = 0.0d0
  do j=1,1609
     do i=2,27
        if(j.eq.1)then
           dye_bin_runningsum(i-1,j) = dye_bin(i,j)*(dye_bin(1,j+1)-dye_bin(1,j)) 
        else
           dye_bin_runningsum(i-1,j) = dye_bin(i,j)*(dye_bin(1,j+1)-dye_bin(1,j))+dye_bin_runningsum(i-1,j-1)
        end if
     end do
  end do

  !time integration of summed dyedt
  dyesum_int = 0.0d0
  do i=1,1609
     if(i.eq.1)then
        dyesum_int(i) = dyesum(i)*(dyesum_time(i+1)-dyesum_time(i))
     else
        dyesum_int(i) = dyesum(i)*(dyesum_time(i+1)-dyesum_time(i))+dyesum_int(i-1)
     end if
  end do
  
  dyesum_int(:) = dyesum_int(:) + 4.222262676d-1

  do i=1,nindex-1
     write(2,*) dye_bin(1,i),dye_bin_runningsum(1,i),dye_bin_runningsum(2,i),dye_bin_runningsum(3,i),dye_bin_runningsum(4,i),dye_bin_runningsum(5,i),dye_bin_runningsum(6,i),dye_bin_runningsum(7,i),dye_bin_runningsum(8,i),dye_bin_runningsum(9,i),dye_bin_runningsum(10,i),dye_bin_runningsum(11,i),dye_bin_runningsum(12,i),dye_bin_runningsum(13,i),dye_bin_runningsum(14,i),dye_bin_runningsum(15,i),dye_bin_runningsum(16,i),dye_bin_runningsum(17,i),dye_bin_runningsum(18,i),dye_bin_runningsum(19,i),dye_bin_runningsum(20,i),dye_bin_runningsum(21,i),dye_bin_runningsum(22,i),dye_bin_runningsum(23,i),dye_bin_runningsum(24,i),dye_bin_runningsum(25,i),dye_bin_runningsum(26,i)
     write(5,*) dye_bin(1,i),dyesum_int(i)
  end do


  close(1)
  close(2)
  close(3)
  close(4)
  close(5)

    
end program point_example
  
