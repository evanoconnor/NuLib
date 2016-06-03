!-*-f90-*- 

!If you are using Hempel EOS, these routines use his, go download them
!from http://phys-merger.physik.unibas.ch/~hempel/eos.html.
!Currently, this is setup to use the sfho_frdm_comp.zip composition
!data, modify appropriately

module nuclei_hempel
  
  implicit none
  integer,allocatable :: hempel_lookup_table(:,:)

contains

  subroutine set_up_Hempel

#if NUCLEI_HEMPEL
    use sfho_frdm_composition_module
    implicit none

    !local
    integer :: highestA, highestZ
    integer :: i

    call compdata_readin

    !highest A
    highestA = maxval(az(1:kmax,1))
    !highest Z
    highestZ = maxval(az(1:kmax,2))

    allocate(hempel_lookup_table(highestA,highestZ))

    hempel_lookup_table = 0

    do i=1,kmax
       hempel_lookup_table(az(i,1),az(i,2)) = i
    enddo

#endif

  end subroutine set_up_Hempel

  !so far, all we need is mass fractions of a select number of nuclei,
  !this routine takes those nuclei and returns the mass fractions and
  !number densities of only those nuclei.

  subroutine nuclei_distribution_Hempel(number_nuclei,nuclei_A,nuclei_Z,mass_fractions,number_densities,eos_variables)

#if NUCLEI_HEMPEL
    use sfho_frdm_composition_module
    use nulib
    implicit none

    !inputs and outputs
    integer, intent(in) :: number_nuclei
    integer, intent(in), dimension(number_nuclei) :: nuclei_A
    integer, intent(in), dimension(number_nuclei) :: nuclei_Z
    real*8, intent(out), dimension(number_nuclei) :: mass_fractions
    real*8, intent(out), dimension(number_nuclei) :: number_densities !number/fm^3
    real*8, intent(in), dimension(total_eos_variables) :: eos_variables

    !local
    real*8 :: t,ye,nb
    real*8, dimension(kmax) :: naz, xaz
    real*8 :: xn,xp,nn,np
    integer :: sflag
    integer :: i,inuc

    t=eos_variables(tempindex)
    ye=eos_variables(yeindex)
    nb=eos_variables(rhoindex)*1.0d-39/(m_ref*mev_to_gram)

    sflag = 0

    if (maxval(nuclei_A).gt.maxval(az(1:kmax,1))) stop "At least one value of A is too high for Hempel's table"
    if (maxval(nuclei_Z).gt.maxval(az(1:kmax,2))) stop "At least one value of Z is too high for Hempel's table"
    call sub_dist_interpol(t,ye,nb,xaz,xn,xp,naz,nn,np,sflag)

    if (sflag.eq.1) then
       !worked as expected
    else
       stop "nuclei_distribution_Hempel: interpolation failed"
    endif

    do i=1,number_nuclei
       inuc = hempel_lookup_table(nuclei_A(i),nuclei_Z(i))
       if (inuc.eq.0) stop "you want a nuclei Hempel doesn't have"
       number_densities(i) = naz(inuc) !number/fm^3
       mass_fractions(i) = xaz(inuc)
    enddo

#endif

  end subroutine nuclei_distribution_Hempel


  !use these routines if you want to probe the full Hempel table
  subroutine get_Hempel_number_of_species(species)

#if NUCLEI_HEMPEL
    use sfho_frdm_composition_module  
    implicit none

    integer, intent(out) :: species

    species = kmax

#endif

  end subroutine get_Hempel_number_of_species

  subroutine get_Hempel_As_and_Zs(As,Zs)

#if NUCLEI_HEMPEL
    use sfho_frdm_composition_module  
    implicit none

    integer, intent(out) :: As(kmax)
    integer, intent(out) :: Zs(kmax)

    As = az(1:kmax,1)
    Zs = az(1:kmax,2)

#endif

  end subroutine get_Hempel_As_and_Zs

end module nuclei_hempel
