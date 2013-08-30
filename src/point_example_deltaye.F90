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
  integer i,j
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
  real*8, dimension(25) :: dye_bin
  real*8 :: dye_bin_sum
  real*8 :: dyedt_hydro
  real*8 :: dyedt_ecfreep
  logical :: spectralogical

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


  open(11,file="src/extra_code_and_tables/rho_temp_ye_c.dat",status='old')
  open(22,file="src/extra_code_and_tables/dye.dat")
  open(33,file="src/extra_code_and_tables/M1_nue_enspectra_cenprime.xg",status='old')
  open(44,file="src/extra_code_and_tables/dyesum.dat")
  open(55,file="src/extra_code_and_tables/dyedt_hydro_c_t.dat",status='old')
  open(66,file="src/extra_code_and_tables/nue_enspectra.dat")
  open(77,file="src/extra_code_and_tables/xp_c_t.dat")


  cont = .true.
  read(33,"(A)") time_string
  nindex = 1
  spectralogical = .true.

  do while(cont)
     !read profile vars from gr1d evolution
     read(55,*) xtime,dyedt_hydro
     read(11,*) xtime,xrho,xtemp,xye

     !match with line from nue_enspectra_cen file
     read(time_string(index(time_string(1:),"= ")+5:index(time_string(1:),"= ")+28),*) current_time
     if(current_time.gt.1.0d0) stop

     do i=1,number_groups
        read(33,*) current_energy,probability_dist(i)
     end do
     read(33,*)
     read(33,*)

     !setup eos_variables
     eos_variables = 0.0d0
     eos_variables(rhoindex) = xrho
     eos_variables(tempindex) = xtemp
     eos_variables(yeindex) = xye

     !! EOS stuff
     keytemp = 1
     keyerr = 0
     call nuc_eos_full(eos_variables(rhoindex),eos_variables(tempindex), &
          eos_variables(yeindex),eos_variables(energyindex),matter_prs, &
          matter_ent,matter_cs2,matter_dedt,matter_dpderho,matter_dpdrhoe, &
          eos_variables(xaindex),eos_variables(xhindex),eos_variables(xnindex), &
          eos_variables(xpindex),eos_variables(abarindex),eos_variables(zbarindex), &
          eos_variables(mueindex),eos_variables(munindex),eos_variables(mupindex), &
          eos_variables(muhatindex),keytemp,keyerr,precision)
     if (keyerr.ne.0) then
        write(*,*) "rho: ", eos_variables(rhoindex)
        write(*,*) "temperature: ", eos_variables(tempindex)
        write(*,*) "ye: ", eos_variables(yeindex)
        write(*,*) "eos error", keyerr
        stop "set_eos_variables: us eos error"
     endif
     if(eos_variables(xhindex).lt.1.0d-15) then
        eos_variables(xhindex) = 0.0d0
     endif

     !calculate the blackbody function (dimensionless fermi function in this 
     !case because other terms cancel in the division of E/B (f_avg / f_neq)     
     do i=1,number_groups 
        blackbody_spectrum(i) = fermidirac_dimensionless(energies(i)/eos_variables(tempindex),&
             (eos_variables(mueindex)+eos_variables(mupindex)-eos_variables(munindex))/eos_variables(tempindex))
     end do




     !begin emissivity calculation for each nucleus
     !Hempel EOS and number of species are set up in readrates
     call nuclei_distribution_Hempel(nspecies,nuclei_A,nuclei_Z,mass_fractions,number_densities,eos_variables)
     emissivity = 0.0d0
     emissivity = 0.0d0
     emissivity_ni56 = 0.0d0
     do i=1,nspecies
        if(i.eq.1)then !if LMP data is not provided for a given nucleus, we will use the rates for 56Ni
           emissivity_ni56(:) = emissivity_from_weak_interaction_rates(56,28,1.0d0,eos_variables,1) 
        endif
        if(nucleus_index(nuclei_A(i),nuclei_Z(i)) == 0)then
           emissivity(i,:) = emissivity_ni56(:)*number_densities(i)
        else
           emissivity(i,:) = emissivity_from_weak_interaction_rates(nuclei_A(i),nuclei_Z(i),number_densities(i),&
                eos_variables,1)
        end if
     end do

     !calculate emissivity for electron capture on free protons
     call single_point_return_all(eos_variables, &
          local_emissivity,local_absopacity,local_scatopacity, &
          mypoint_neutrino_scheme)
     

     !calculate dyedt for each nucleus
     dye = 0.0d0
     do i=1,nspecies
!        dye(i)=-Sum((4.0d0*pi/6.02214129d23)*bin_widths(:)*(emissivity(i,:)/mev_to_erg)*(1-probability_dist(:)/blackbody_spectrum(:))/energies(:))
        dye(i)=-Sum((4.0d0*pi/6.02214129d23)*bin_widths(:)*(emissivity(i,:)/mev_to_erg)*probability_dist(:)/energies(:))
!        if (nindex.ge.1605) write(*,*) Sum(emissivity(i,:)), mass_fractions(i)
     end do     
     normalized_dye(:)=dye(:)/eos_variables(rhoindex)
    
     !calculate dyedt for free protons
     dyedt_ecfreep = -Sum((4.0d0*pi/6.02214129d23)*bin_widths(:)*(local_emissivity(1,:)/mev_to_erg)*probability_dist(:)/energies(:))/eos_variables(rhoindex)
     

     !write to file
     dye_bin = 0.0d0
     dye_bin_sum = 0.0d0
     do i=1,nspecies
        if(nuclei_A(i).le.5) then
           dye_bin(1) = dye_bin(1) + normalized_dye(i)
        else if(nuclei_A(i).ge.5.and.nuclei_A(i).lt.25) then
           dye_bin(2) = dye_bin(2) + normalized_dye(i)
        else if(nuclei_A(i).ge.25.and.nuclei_A(i).lt.45) then 
           dye_bin(3) = dye_bin(3) + normalized_dye(i)
        else if(nuclei_A(i).ge.45.and.nuclei_A(i).le.65) then
           dye_bin(4) = dye_bin(4) + normalized_dye(i)
        else if(nuclei_A(i).gt.65.and.nuclei_A(i).lt.85) then
           dye_bin(5) = dye_bin(5) + normalized_dye(i)
        else if(nuclei_A(i).ge.85.and.nuclei_A(i).lt.105) then
           dye_bin(6) = dye_bin(6) + normalized_dye(i)
        else if(nuclei_A(i).ge.105.and.nuclei_A(i).lt.125) then
           dye_bin(7) = dye_bin(7) + normalized_dye(i)
        else if(nuclei_A(i).ge.125.and.nuclei_A(i).lt.145) then
           dye_bin(8) = dye_bin(8) + normalized_dye(i)
        else if(nuclei_A(i).ge.145.and.nuclei_A(i).lt.165) then
           dye_bin(9) = dye_bin(9) + normalized_dye(i)
        else if(nuclei_A(i).ge.165.and.nuclei_A(i).lt.185) then
           dye_bin(10) = dye_bin(10) + normalized_dye(i)
        else if(nuclei_A(i).ge.185.and.nuclei_A(i).lt.205) then
           dye_bin(11) = dye_bin(11) + normalized_dye(i)
        else if(nuclei_A(i).ge.205.and.nuclei_A(i).lt.225) then
           dye_bin(12) = dye_bin(12) + normalized_dye(i)
        else if(nuclei_A(i).ge.225.and.nuclei_A(i).lt.245) then
           dye_bin(13) = dye_bin(13) + normalized_dye(i)
        else if(nuclei_A(i).ge.245.and.nuclei_A(i).lt.265) then
           dye_bin(14) = dye_bin(14) + normalized_dye(i)
        else if(nuclei_A(i).ge.265.and.nuclei_A(i).lt.285) then
           dye_bin(15) = dye_bin(15) + normalized_dye(i)
        else if(nuclei_A(i).ge.285.and.nuclei_A(i).lt.305) then
           dye_bin(16) = dye_bin(16) + normalized_dye(i)
        else if(nuclei_A(i).ge.305.and.nuclei_A(i).lt.325) then
           dye_bin(17) = dye_bin(17) + normalized_dye(i)
        else if(nuclei_A(i).ge.325.and.nuclei_A(i).lt.345) then
           dye_bin(18) = dye_bin(18) + normalized_dye(i)
        else if(nuclei_A(i).ge.345.and.nuclei_A(i).lt.365) then
           dye_bin(19) = dye_bin(19) + normalized_dye(i)
        else if(nuclei_A(i).ge.365.and.nuclei_A(i).lt.385) then
           dye_bin(20) = dye_bin(20) + normalized_dye(i)
        else if(nuclei_A(i).ge.385.and.nuclei_A(i).lt.405) then
           dye_bin(21) = dye_bin(21) + normalized_dye(i)
        else if(nuclei_A(i).ge.425.and.nuclei_A(i).lt.445) then
           dye_bin(22) = dye_bin(22) + normalized_dye(i)
        else if(nuclei_A(i).ge.445.and.nuclei_A(i).lt.465) then
           dye_bin(23) = dye_bin(23) + normalized_dye(i)
        else if(nuclei_A(i).ge.465.and.nuclei_A(i).lt.485) then
           dye_bin(24) = dye_bin(24) + normalized_dye(i)
        else if(nuclei_A(i).ge.485.and.nuclei_A(i).lt.505) then
           dye_bin(25) = dye_bin(25) + normalized_dye(i)
        end if         
     end do
     dye_bin_sum = Sum(dye_bin(:)) + dyedt_ecfreep

     !write to binned file
     write(22,*) xtime,dye_bin(1),dye_bin(2),dye_bin(3),dye_bin(4),dye_bin(5),dye_bin(6),dye_bin(7),dye_bin(8),dye_bin(9),dye_bin(10),dye_bin(11),dye_bin(12),dye_bin(13),dye_bin(14),dye_bin(15),dye_bin(16),dye_bin(17),dye_bin(18),dye_bin(19),dye_bin(20),dye_bin(21),dye_bin(22),dye_bin(23),dye_bin(24),dye_bin(25),dyedt_ecfreep

     !write to sum file
     write(44,*) xtime,dye_bin_sum

     !write to sum file
     write(77,*) xtime,eos_variables(xpindex)
     
     if(xtime.ge.8.5d-2.and.spectralogical) then
        write(66,*) xtime
        do i=1,number_groups
           write(66,*) energies(i) , probability_dist(i), blackbody_spectrum(i), (1-probability_dist(i)/blackbody_spectrum(i))
        end do
        spectralogical = .false.
        close(66)
     endif

     !read next time and stop if at EOF
     nindex = nindex + 1
     write(*,*) nindex
     read(33,"(A)",IOSTAT=IO) time_string
     if(Sum(mass_fractions(:)).le.0.1d0) cont = .false.
     if (IO.lt.0) cont = .false.     
  enddo

  close(11)
  close(22)
  close(33)
  close(44)
  close(55)
  close(77)
    
end program point_example
  
