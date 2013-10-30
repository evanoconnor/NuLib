!-*-F90-*- 
program test

  use weak_rates
  use nulib

  implicit none
  
  real*8 :: t9
  real*8 :: lrhoye
  real*8 :: analytic_weakrates,rate,avge(2)
  integer i,j,t,d


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
  real*8 :: xrho, xtemp, xye
  real*8 dxfac,mindx
  real*8, dimension(4) :: lrhoyearr,t9arr
  real*8, dimension(4,4) :: muearr
  real*8 :: qec_eff

  !allocate the arrays for the point values
  allocate(local_emissivity(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_absopacity(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_scatopacity(mypoint_number_output_species,mypoint_number_groups))
  allocate(blackbody_spectra(mypoint_number_output_species,mypoint_number_groups))

  call input_parser(parameters)
  !this sets up many cooefficients and creates the energy grid (one
  !zone + log spacing) see nulib.F90:initialize_nulib
  call initialize_nulib(mypoint_neutrino_scheme,mypoint_number_species,mypoint_number_groups)

  !read in EOS table & set reference mass
  call readtable(eos_filename)
  m_ref = m_amu !for SFHo_EOS (Hempel)
  ! m_ref = m_n !for LS220

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

  !! EOS stuff
  keytemp = 1
  keyerr = 0




!###########################
!######END NULIB SETUP######
!###########################

  open(33,file="analysis/qec_solver_gpq_test.dat")

  call readrates(table_bounds)
  
  lrhoyearr(1)=10.0d0
  lrhoyearr(2)=10.5d0
  lrhoyearr(3)=11.5d0
  lrhoyearr(4)=12.4d0

  t9arr(1)=8.20d0
  t9arr(2)=20.0d0
  t9arr(3)=30.0d0
  t9arr(4)=39.0d0

  rate=0.0d0

  do i=1,78
     write(33,"(F10.6,A)",advance='no') nuclear_species(i,1)-2.50d0, " "
     do t=1,4    
        do d=1,4
           xye = 0.5d0 !dimensionless
           xrho = (10**(lrhoyearr(d)))/xye !g/cm^3
           xtemp = t9arr(t)*1.0d9*kelvin_to_mev !MeV
           eos_variables = 0.0d0
           eos_variables(rhoindex) = xrho
           eos_variables(tempindex) = xtemp
           eos_variables(yeindex) = xye
           call nuc_eos_full(eos_variables(rhoindex),eos_variables(tempindex), &
                eos_variables(yeindex),eos_variables(energyindex),matter_prs, &
                matter_ent,matter_cs2,matter_dedt,matter_dpderho,matter_dpdrhoe, &
                eos_variables(xaindex),eos_variables(xhindex),eos_variables(xnindex), &
                eos_variables(xpindex),eos_variables(abarindex),eos_variables(zbarindex), &
                eos_variables(mueindex),eos_variables(munindex),eos_variables(mupindex), &
                eos_variables(muhatindex),keytemp,keyerr,precision)

           muearr(t,d)=eos_variables(mueindex)
           avge(2)=nuclear_species(i,1) 
           avge(1)=analytic_weakrates(1,(t9arr(t)*10**9)*kelvin_to_mev,avge(2),muearr(t,d))/analytic_weakrates(0,(t9arr(t)*10**9)*kelvin_to_mev,avge(2),muearr(t,d))
           qec_eff = qec_solver(avge,avge(2),eos_variables)
           write(33,"(F15.6,A)",advance='no') qec_eff, " "
        end do
     end do
     write(33,*) " "
  end do


  close(33)
  

end program test
