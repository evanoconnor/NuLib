!-*-f90-*-
module weakrates_interface
  !""" Main interface between weakratelib and nulib """
    
  use class_ratelibrary
  
  implicit none
  private
  public :: initialize_weakratelib, microphysical_electron_capture

  ! members (singleton)
  type(RateLibrary) :: weakratelib
  !$OMP THREADPRIVATE(weakratelib)

  ! methods
contains

!------------------------------------------------------------------------------------!
  
  subroutine initialize_weakratelib(parameters_filename)
    character*(*) parameters_filename
    !$OMP PARALLEL COPYIN(weakratelib)
    weakratelib = new_RateLibrary(parameters_filename)
    !$OMP END PARALLEL
  end subroutine initialize_weakratelib

!------------------------------------------------------------------------------------!

  function emissivity_from_weak_interaction_rates(&
       A,&
       Z,&
       number_density,&
       eos_variables,&
       neutrino_species,&
       idxtable)&
       result(emissivity)

    use class_ratelibrary
    use class_rateapproximation
    use nulib, only : total_eos_variables,energies,number_groups,&
         do_integrated_BB_and_emissivity,mueindex,rhoindex,tempindex,&
         yeindex,GPQ_n32_roots,GPQ_n32_weights
      
    include 'constants.inc'

    integer A,Z
    real*8, intent(in) :: eos_variables(total_eos_variables)
    real*8, intent(in) :: number_density
    integer, intent(in) :: neutrino_species
    integer,intent(in) :: idxtable


    logical :: approx_rate_flag
    real*8 :: emissivity(number_groups) !final answer in erg/cm^/s/srad/MeV
    real*8 :: avgenergy(2)
    real*8 :: GPQ_interval(2)
    real*8 :: GPQ_coef(2)
    real*8 :: qec_eff                   !effective Qec for the approximate neutrino spectra

    !spectrum integration
    real*8 :: normalization_constant    !nu spectra normalization, units of 1/MeV^5/s
    real*8 :: spectra
    real*8 :: t9
    real*8 :: lrhoYe
    real*8 :: nu_spectrum_eval
    integer i,ng,j, omp_get_thread_num

    !local rate variables
    real*8 :: rbeta      !beta decay rate (plus or minus depending on neutrino species)
    real*8 :: rcap       !capture rate (electron or positron for nue or anue)
    real*8 :: rnu        !nue or anue energy loss rate


    approx_rate_flag = .false.
    GPQ_interval = 0.0d0
    GPQ_coef(:) = 0.0d0
    !set local eos_variables for rate interpolation
    lrhoYe = log10(eos_variables(rhoindex)*eos_variables(yeindex))
    t9 = (eos_variables(tempindex)/kelvin_to_mev)/(10.0d0**9.0d0)  ! Conversion from MeV to GK

    ! if no table contains rate, use approximate routine
    if (idxtable.eq.0) then
       approx_rate_flag = .true.
    endif

    if(approx_rate_flag) then
       qec_eff = return_hempel_qec(A,Z,Z-1)
    else
       !interpolating rates for given eos_variables and calculating 
       !average neutrino energy from rates for nue, emissivities are
       !from the betaplus direction; for anue, emissivities in the betaminus direction
       if (neutrino_species.eq.1) then
          rbeta = return_weakrate(weakratelib,A,Z,t9,lrhoYe,idxtable,1)
          rcap = return_weakrate(weakratelib,A,Z,t9,lrhoYe,idxtable,2)
          rnu = return_weakrate(weakratelib,A,Z,t9,lrhoYe,idxtable,3)
          !using Qgs from table as seed for qec_solver
          qec_eff = weakratelib%tables(idxtable)%nuclear_species(weakratelib%tables(idxtable)%nucleus_index(A,Z),1) 
          avgenergy(1) = rnu/(rcap + rbeta) 
          avgenergy(2) = qec_eff !necessary to fulfill the first comparison in qec_solver
       else if (neutrino_species.eq.2) then
          rbeta = return_weakrate(weakratelib,A,Z+1,t9,lrhoYe,idxtable,4)
          rcap = return_weakrate(weakratelib,A,Z+1,t9,lrhoYe,idxtable,5)
          rnu = return_weakrate(weakratelib,A,Z+1,t9,lrhoYe,idxtable,6)
          qec_eff = -weakratelib%tables(idxtable)%nuclear_species(weakratelib%tables(idxtable)%nucleus_index(A,Z+1),1)
          avgenergy(1) = rnu/(rcap + rbeta)   
          avgenergy(2) = qec_eff
       else
          stop "This module only implements electron-neutrino type weak interactions. &
               Please restrict neutrino_species to nue/anue."
       endif

       !solve for the eff Qec that constrains the effective 
       !neutrino spectra to produce the correct avg. energy
       qec_eff = qec_solver(avgenergy,qec_eff,eos_variables)
    end if

    !calculate normalization constant using effective neutrino spectra
    spectra = 0.0d0
    GPQ_interval = GPQ_intervals(qec_eff,eos_variables)   ! dynamic range finder for matching
    GPQ_coef(1) = (GPQ_interval(2)-GPQ_interval(1))/2.0d0 ! integration range to the neutrino
    GPQ_coef(2) = (GPQ_interval(1)+GPQ_interval(2))/2.0d0 ! spectrum width
    do i=1,32
       spectra = spectra + &
            GPQ_coef(1)*GPQ_n32_weights(i)*&
            ec_neutrino_spectra(&
            GPQ_coef(1)*GPQ_n32_roots(i)+GPQ_coef(2),&
            qec_eff,eos_variables(mueindex)-m_e,eos_variables(tempindex))
    end do
    spectra = (eos_variables(tempindex)**5)*spectra          
    if (neutrino_species.eq.1) then
       if (approx_rate_flag) then
          normalization_constant = (return_weakrate(&
               0,eos_variables(tempindex),qec_eff,&
               eos_variables(mueindex)-m_e))/spectra 
       else
          normalization_constant = (rbeta+rcap)/spectra
       end if
    else if (neutrino_species.eq.2) then
       if (approx_rate_flag) then
          !approximation doesn't work for antineutrinos
          normalization_constant = 0.0
       else
          normalization_constant = (rbeta+rcap)/spectra
       end if
    end if

    !integrate over energy bin or take central value of bin and multiply by the bin width
    if (do_integrated_BB_and_emissivity) then
       stop "Integration is not yet supported for electron-capture emissivities" 
    else
       do ng=1,number_groups
          nu_spectrum_eval =&
               (eos_variables(tempindex)**4.0d0)*normalization_constant*&
               ec_neutrino_spectra(energies(ng)/eos_variables(tempindex),&
               qec_eff,eos_variables(mueindex)-m_e,eos_variables(tempindex))
          emissivity(ng) = &  !erg/cm^3/s/MeV/srad
               (energies(ng)*mev_to_erg)*(1.0d39*number_density)*nu_spectrum_eval/(4.0d0*pi) 
          ! spectra can be 0 and normalizaton_constant Inf. which results in NaN emissivity.
          if(isnan(emissivity(ng))) then
               emissivity(ng) = 0
          endif
       end do
    endif
    return

  end function emissivity_from_weak_interaction_rates

!------------------------------------------------------------------------------------!

  function ec_neutrino_spectra(nu_energy_per_T,q,uf,T) result(nu_spectra)
    !definition: n(E) = (T^4)*N*ec_neutrino_spectra, where N is the normalization constant
    ![N] = #/MeV^5/s
    include 'constants.inc'

    !local variables
    real*8 :: T
    real*8 :: nu_energy_per_T
    real*8, intent(in) :: q
    real*8 :: uf
    real*8 :: nu_spectra
    real*8 :: q_T
    real*8 :: uf_T

    !redefining q and the chemical potential to be dimensionless, T must be in MeV
    q_T = q/T
    uf_T = uf/T  

    !ec neutrino spectra (dimensionless)
    if ((nu_energy_per_T-q_T).le.0.0d0) then
       nu_spectra = 0.0d0
    else
       nu_spectra = (nu_energy_per_T**2.0d0)*((nu_energy_per_T-q_T)**2.0d0) &
            /(1.0d0+exp(nu_energy_per_T-q_T-uf_T))
    end if
    !note a factor of T^4 was removed to make the spectrum dimensionless, the returned result
    !should therefore be multiplied by T^4 (T in MeV)

  end function ec_neutrino_spectra

!------------------------------------------------------------------------------------!
  
  function qec_solver(avgenergy,qin,eos_variables) result(qec_eff)

    use nulib, only : GLQ_n32_roots, GLQ_n32_weights, &
         total_eos_variables, mueindex, tempindex, rhoindex, yeindex

    implicit none
    include 'constants.inc'

    !starting point energies using a seed q = Qgs
    real*8, intent(in) :: avgenergy(2)
    real*8, intent(in) :: eos_variables(total_eos_variables)

    !local buffer variables for calculations
    real*8 :: avge_rates
    real*8 :: avge_spectra
    real*8 :: qec_eff
    real*8 :: q_newton_raphson
    real*8 :: q
    real*8, intent(in)  :: qin
    integer i,N

    !integration variables
    real*8 :: energy_density_integral
    real*8 :: number_density_integral
    real*8 :: dq_energy_density_integral
    real*8 :: dq_number_density_integral
    real*8 :: nu_spectra_deriv_coef
    real*8 :: nu_spectra

    !bisection method declerations
    real*8 :: lower_bound
    real*8 :: upper_bound
    real*8 :: avge_spectra_boundary
    real*8 :: tolerance
    real*8 :: GPQ_interval(2)
    real*8 :: GPQ_interval_deriv(2)
    real*8 :: GPQ_coef(2)
    real*8 :: GPQ_coef_deriv(2)
    real*8 :: energy

    integer limiter,nmax_bisections,extrapolation

    avge_rates = avgenergy(1)
    avge_spectra = avgenergy(2)
    tolerance = 0.1d0
    nmax_bisections = 15
    limiter = 0
    q = qin
    qec_eff = q
    GPQ_interval = 0.0d0
    GPQ_interval_deriv = 0.0d0
    GPQ_coef = 0.0d0
    GPQ_coef_deriv = 0.0d0

    do while (abs((avge_rates - avge_spectra)/avge_rates) > 1.0d-8) 
       lower_bound = -100.0d0
       upper_bound = 100.0d0
       avge_spectra_boundary = average_energy(lower_bound,eos_variables)
       N=0
       nmax_bisections = 1000
       tolerance = 1.0d-8

       !low resolution in the shell-model rates can cause the interpolation to produce an
       !average nu energy below the asymptotic limit of <E>_spectra. In this case, the
       !effective q is set to the lower bound -100MeV
       if(avge_spectra_boundary.gt.avge_rates) then          
          qec_eff = lower_bound
          return
       end if
       do !bisection loop
          if(N.ge.900)then                   
             stop "Over 900 bisections, failed to converge."
          endif
          !if the bisection parameter q exceeds double precision
          if(q.eq.(lower_bound + upper_bound)/2.0d0)then
             if(abs(avge_rates - avge_spectra)/avge_rates.lt.1.0d0)then  
                qec_eff = q
                return
             end if
          end if
          q = (lower_bound + upper_bound)/2.0d0
          avge_spectra = average_energy(q,eos_variables)                   
          if (abs(avge_rates - avge_spectra)/avge_rates.le.tolerance) then
             qec_eff = q
             return 
          end if
          if (sign(1.0d0,(avge_rates-avge_spectra)).eq.&
               sign(1.0d0,(avge_rates-avge_spectra_boundary))) then
             lower_bound = q
          else
             upper_bound = q
          end if
          N = N + 1   
       end do
    end do
    qec_eff = q
    return

  end function qec_solver

!------------------------------------------------------------------------------------!
  
  function average_energy(qec_eff,eos_variables) result(avgenergy)

    use nulib, only : GPQ_n32_roots,GPQ_n32_weights,&
         total_eos_variables,m_e,mueindex,tempindex

    implicit none

    real*8, intent(in) :: eos_variables(total_eos_variables)
    real*8, intent(in) :: qec_eff

    !local integration variables
    integer :: i
    real*8 :: avgenergy
    real*8 :: spectra
    real*8 :: number_density_integral
    real*8 :: energy_density_integral
    real*8 :: GPQ_interval(2)
    real*8 :: GPQ_coef(2)

    GPQ_interval = 0.0d0
    GPQ_coef = 0.0d0

    !set weights and roots for quadrature integration, then calculate <E>
    avgenergy=0.0d0
    energy_density_integral=0.0d0
    number_density_integral=0.0d0
    spectra = 0.0d0

    GPQ_interval = GPQ_intervals(qec_eff,eos_variables)
    GPQ_coef(1) = (GPQ_interval(2)-GPQ_interval(1))/2.0d0
    GPQ_coef(2) = (GPQ_interval(1)+GPQ_interval(2))/2.0d0
    do i=1,32
       spectra = GPQ_coef(1)*ec_neutrino_spectra(&
            GPQ_coef(1)*GPQ_n32_roots(i)+GPQ_coef(2),&
            qec_eff,eos_variables(mueindex)-m_e,eos_variables(tempindex))
       energy_density_integral = energy_density_integral + &
            GPQ_n32_weights(i)*(GPQ_coef(1)*GPQ_n32_roots(i)+GPQ_coef(2))*spectra
       number_density_integral = number_density_integral + GPQ_n32_weights(i)*spectra
    end do

    if (number_density_integral.eq.0.0d0.and.energy_density_integral.eq.0.0d0) then             
       avgenergy = (eos_variables(tempindex))*1.0d0
    else
       avgenergy = (eos_variables(tempindex))*&
            (energy_density_integral/number_density_integral)
    end if

  end function average_energy

!------------------------------------------------------------------------------------!

  function GPQ_intervals(q,eos_variables) result (interval)
    use nulib, only : total_eos_variables,tempindex,mueindex,m_e
    implicit none

    real*8, dimension(2) :: interval
    real*8, intent(in) :: eos_variables(total_eos_variables)
    real*8, intent(in) :: q
    real*8 :: spectra
    real*8 :: energy
    real*8 :: centroid
    real*8 :: spectra_centroid
    real*8 :: spectra_shift
    real*8 :: lower_bound,lower_bound_spectra
    real*8 :: upper_bound,upper_bound_spectra
    real*8 :: tolerance
    integer :: N,nmax_bisections,i

    interval = 0.0d0
    spectra = 0.0d0
    !lowerbound
    centroid = (q+eos_variables(mueindex)-m_e)/eos_variables(tempindex)
    spectra_centroid = &
         ec_neutrino_spectra(centroid,q,eos_variables(mueindex)-m_e,eos_variables(tempindex))
    spectra_shift = spectra_centroid/1.0d2          
    spectra_shift = sqrt(spectra_shift)
    !integration range-finder:
    !If (q+mu)/T is less than 3MeV, the spectra will hug the origin,
    !so an interval from 0 to 25 should be sufficient. For greater values
    !of (q+mu)/T the interval should be determined with the quartic part 
    !of the spectra for the lower bound and an auxillary gaussian that is 
    !used to track the high energy bound. The below represents this method.
    if(centroid.gt.3.0d0)then
       interval(1) = (q/eos_variables(tempindex)+sqrt((q/eos_variables(tempindex))**2.0d0 + &
            4.0d0*spectra_shift))/2.0d0         
       interval(2) = (q+eos_variables(mueindex)-m_e)/eos_variables(tempindex) + 27.314d0
    else
       interval(1) = 0.0d0
       interval(2) = 25.0d0
    end if
    ! fail safe in case the above doesn't work 
    if(interval(1).lt.1.0d-5.and.interval(2).lt.0.0d0)then
       interval(2) = max(interval(2),50.0d0)
    end if
    if(interval(2).lt.0.0d0)then
       write(*,*)  "lower bound = ",interval(1), "upper bound = ", interval(2)
       write(*,*) 
       stop
    end if

    return

  end function GPQ_intervals

!------------------------------------------------------------------------------------!
  
  subroutine microphysical_electron_capture(&
       neutrino_species,&
       eos_variables,&
       emissivity)

    use nuclei_hempel
    use nulib, only : total_eos_variables, number_groups, &
         tempindex, mueindex, rhoindex, yeindex, kelvin_to_mev
    use sfho_frdm_composition_module, only : sfho_mass

    integer i
    integer, intent(in) :: neutrino_species
    real*8, intent(in) :: eos_variables(total_eos_variables)
    real*8, dimension(number_groups) :: emissivity
    real*8, dimension(number_groups) :: emissivity_temp
    real*8 :: logrhoYe,t9
    integer :: idxtable, A, Z ! if idxtable=0 here, omp fails

    idxtable = 0

    !Hempel EOS and number of species are set up in readrates
    call nuclei_distribution_Hempel(&
         weakratelib%approx%nspecies,weakratelib%approx%nuclei_A,&
         weakratelib%approx%nuclei_Z,weakratelib%approx%mass_fractions,&
         weakratelib%approx%number_densities,eos_variables)

    emissivity = 0.0d0
    logrhoYe = log10(eos_variables(rhoindex)*eos_variables(yeindex))
    t9 = (eos_variables(tempindex)/kelvin_to_mev)*1.0d-9
    
    do i=1,weakratelib%approx%nspecies 

       if(weakratelib%approx%number_densities(i).eq.0.0d0) cycle

       A = weakratelib%approx%nuclei_A(i)
       Z = weakratelib%approx%nuclei_Z(i)

       !if rate data from a table is not present and A>4 with iapprox nonzero,
       !use the parameterized rate function, else skip this nucleus
       ! check if rate exists in a table
       if (weakratelib%ntables.ne.0)then
          if (neutrino_species .eq. 1) idxtable = in_table(weakratelib,A,Z,logrhoYe,t9)
          if (neutrino_species .eq. 2) idxtable = in_table(weakratelib,A,Z+1,logrhoYe,t9)
       endif

       ! if no table contains the requested rate
       if(idxtable.eq.0) then
          ! and if the approximation is turned on
          if(weakratelib%priority(5).gt.0) then
             ! and the nucleus is above the A=4 isobars and below
             if (A.gt.4) then
             else
                cycle
             end if
             ! check to make sure masses exist for both parent and daughter nucleus
             if(hempel_lookup_table(A,Z).eq.0.or.hempel_lookup_table(A,Z-1).eq.0) then
                cycle
             end if             
          else
             cycle
          end if
       end if
       

       !emissivity calculation
       emissivity_temp = emissivity_from_weak_interaction_rates(A,Z,&
            weakratelib%approx%number_densities(i),eos_variables,&
            neutrino_species,idxtable)

       !add the calculation to the total persistent emissivity array for this eos_variables
       emissivity = emissivity + emissivity_temp

    end do
    weakratelib%approx%number_densities = 0.0d0
    weakratelib%approx%mass_fractions = 0.0d0
    return

  end subroutine microphysical_electron_capture
!------------------------------------------------------------------------------------!
end module weakrates_interface
