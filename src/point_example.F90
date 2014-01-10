!-*-f90-*-
program point_example

  use nulib
  use weak_rates
  use sfho_frdm_composition_module, only : sfho_mass
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

  character*200 :: parameters = "/projects/ceclub/gr1dnulib/GitHub/NuLib/parameters"

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

  real*8, allocatable,dimension(:) :: probability_dist
  real*8, allocatable,dimension(:) :: emissivity_ni56
  real*8, allocatable,dimension(:) :: blackbody_spectrum
  real*8, allocatable,dimension(:,:) :: emissivity,emissivity_temp,nuclei_weakrates
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
  real*8 :: qec_eff
  real*8 :: t9
  real*8 :: lrhoye
  real*8 :: analytic_weakrates
  logical :: spectralogical
  logical :: parameterized_rate          


  !allocate the arrays for the point values
  allocate(local_emissivity(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_absopacity(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_scatopacity(mypoint_number_output_species,mypoint_number_groups))
  allocate(blackbody_spectra(mypoint_number_output_species,mypoint_number_groups))

  call input_parser(parameters)

  !this sets up many cooefficients and creates the energy grid (one
  !zone + log spacing) see nulib.F90:initialize_nulib
  call initialize_nulib(mypoint_neutrino_scheme,mypoint_number_species,mypoint_number_groups)
  call set_up_Hempel ! set's up EOS for nuclear abundances
  !read in EOS table & set reference mass
  call readtable(eos_filename)
  m_ref = m_amu !for SFHo_EOS (Hempel)
  ! m_ref = m_n !for LS220

  !read in weak rates table and build interpolant functions
  call readrates(table_bounds)

  allocate(probability_dist(number_groups))
  allocate(emissivity(nspecies,number_groups))
  allocate(nuclei_weakrates(2,nspecies))
  allocate(emissivity_temp(nspecies,number_groups))
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


  open(11,file="/mnt/simulations/ceclub/sullivan/gr1d/dyedtruns/Data_all_vel_inel_12_tvd/rho_temp_ye_c.dat",status='old')
  open(88,file="/mnt/simulations/ceclub/sullivan/gr1d/dyedtruns/Data_all_vel_inel_12_tvd/ecrates_rho_c_temp.dat")

  cont = .true.
  nindex = 1
  spectralogical = .true.
  
  do while(cont)
     !read profile vars from gr1d evolution
     read(11,*,end=101) xtime,xrho,xtemp,xye

     !setup eos_variables
     eos_variables = 0.0d0
     eos_variables(rhoindex) = xrho
     eos_variables(tempindex) = xtemp
     eos_variables(yeindex) = xye
     lrhoYe = log10(eos_variables(rhoindex)*eos_variables(yeindex))
     t9 = (eos_variables(tempindex)/kelvin_to_mev)/(10.0d0**9.0d0)  


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

     !begin emissivity calculation for each nucleus
     !Hempel EOS and number of species are set up in readrates
     call nuclei_distribution_Hempel(nspecies,nuclei_A,nuclei_Z,mass_fractions,number_densities,eos_variables)          

     emissivity = 0.0d0

     do i=1,nspecies
        parameterized_rate = .false.

        !if rate data from a table is not present and 65<A<120 and iapprox is on, use the parameterized rate function, else skip this nucleus
        if(nucleus_index(nuclei_A(i),nuclei_Z(i)).eq.0) then
           if (nuclei_A(i).gt.4) then
              if(file_priority(5).gt.0) then
                 parameterized_rate = .true.
              else
                 cycle
              end if
           else
              cycle
           end if
        end if
        
        !oda table bounds
        if(nucleus_index(nuclei_A(i),nuclei_Z(i)).gt.0.and.nuclei_A(i).lt.40) then !should work for lmp + lmsh, will need to check for other table combinations
           !test to see if outside oda table bounds
           if(lrhoYe.lt.1.0d0.or.lrhoYe.gt.11.0d0.or.t9.lt.1.0d-2.or.t9.gt.30.0d0) then
              if(file_priority(5).gt.0) then
                 parameterized_rate = .true.
              else
                 cycle
              end if
           end if
        else if(nucleus_index(nuclei_A(i),nuclei_Z(i)).gt.0.and.(nuclei_A(i).gt.40.and.nuclei_A(i).le.65)) then !should work for lmp + lmsh, will need to check for other table combinations
           !test to see if outside lmp table bounds
           if(lrhoYe.lt.1.0d0.or.lrhoYe.gt.12.5d0.or.t9.lt.1.0d-2.or.t9.gt.100.0d0) then
              if(file_priority(5).gt.0) then
                 parameterized_rate = .true.
              else
                 cycle
              end if
           end if
        else if(nucleus_index(nuclei_A(i),nuclei_Z(i)).gt.0.and.nuclei_A(i).gt.65) then !should work for lmp + lmsh
           !test to see if outside lmsh table bounds
           if(lrhoYe.lt.9.28493d0.or.lrhoYe.gt.12.42218d0.or.t9.lt.8.12315d0.or.t9.gt.39.04914d0) then
              if(file_priority(5).gt.0) then
                 parameterized_rate = .true.
              else
                 cycle
              end if
           end if
        end if

        if(parameterized_rate)then
           if(sfho_mass(hempel_lookup_table(nuclei_A(i),nuclei_Z(i))).eq.0.0d0.or.sfho_mass(hempel_lookup_table(nuclei_A(i),nuclei_Z(i)-1)).eq.0.0d0) then
              cycle
           end if
        end if


        if(parameterized_rate)then
           qec_eff = return_hempel_qec(nuclei_A(i),nuclei_Z(i),nuclei_Z(i)-1)
           nuclei_weakrates(1,i) = analytic_weakrates(0,eos_variables(tempindex),qec_eff,eos_variables(mueindex))
           nuclei_weakrates(2,i) = analytic_weakrates(1,eos_variables(tempindex),qec_eff,eos_variables(mueindex))
        else
           nuclei_weakrates(1,i) = 10.0d0**weakrates(nuclei_A(i),nuclei_Z(i),t9,lrhoYe,2)
           nuclei_weakrates(2,i) = 10.0d0**weakrates(nuclei_A(i),nuclei_Z(i),t9,lrhoYe,3)
        end if

     end do

     write(*,*)eos_variables(rhoindex),Sum(nuclei_weakrates(1,:))

     
     !read next time and stop if at EOF
     nindex = nindex + 1
     write(*,*) nindex
     if(Sum(mass_fractions(:)).le.0.1d0) cont = .false.
     if (IO.lt.0) cont = .false.     
     
  end do

101  close(11)
  close(88)
  
    
end program point_example
  
