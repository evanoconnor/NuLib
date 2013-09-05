!-*-f90-*-
program driver
  
  use nulibtable
  implicit none

  character(len=1024) :: filename = "/Users/evanoc/research/git/NuLib/NuLib_LS220_rho10_temp10_ye10_ng24_ns3_version1.0_20111024.h5"

  real*8, allocatable :: eas(:)
  real*8, allocatable :: eas_energy(:,:)
  real*8, allocatable :: eas_species_energy(:,:,:)

  integer i,j,k,ns,ng
  real*8 rho,temp,ye

  call nulibtable_reader(filename)

  allocate(eas(nulibtable_number_easvariables))
  allocate(eas_energy(nulibtable_number_groups,nulibtable_number_easvariables))
  allocate(eas_species_energy(nulibtable_number_species, &
       nulibtable_number_groups,nulibtable_number_easvariables))

  write(*,*) "Total species:", nulibtable_number_species
  write(*,*) "Total groups:", nulibtable_number_groups

  write(*,*) "Energies:", nulibtable_energies
  write(*,*) "Inverse Energies:", nulibtable_inv_energies
  write(*,*) "bin widths:", nulibtable_ewidths
  write(*,*) "bin bottoms:", nulibtable_ebottom
  write(*,*) "bin tops:", nulibtable_etop
  
  write(*,*) "Units are: emissivity, ergs/cm^3/s/srad (I multiplied through the bin width for you..)"
  write(*,*) "Units are: absorption opacity, cm^-1"
  write(*,*) "Units are: scattering opacity, cm^-1"

  if (.true.) then
     !example of single energy single species call, here I loop over rho
     write(*,*) "rho,temp,ye, eas variables"
     do i=1,100
        eas(1:3) = 0.0d0
        rho = 10.0d0**(8.0d0+6.0d0*dble(i)/100.0d0)
        temp = 1.0d0
        ye = 0.35d0
        ns = 1
        ng = 5
        call nulibtable_single_species_single_energy(rho,temp,ye,ns,ng,eas, &
             nulibtable_number_easvariables)
        write(*,"(1P10E18.9)") rho,temp,ye,eas
     enddo
  endif

  write(*,*) 
  write(*,*) 

  if (.true.) then
     !example of all energy, single species call.
     rho = 1.0d10
     temp = 1.0d0
     ye = 0.35d0
     ns = 1
     call nulibtable_single_species_range_energy(rho,temp,ye,ns,eas_energy, &
          nulibtable_number_groups,nulibtable_number_easvariables)
     
     write(*,*) "rho,temp,ye:",rho,temp,ye
     write(*,*) "data below, index, energy, eas variables"

     do ng=1,nulibtable_number_groups
        write(*,"(i4,1P10E18.9)") ng,nulibtable_energies(ng), &
             (eas_energy(ng,j),j=1,nulibtable_number_easvariables)
     enddo
  endif

  write(*,*) 
  write(*,*) 

  if (.true.) then
     !example of all energy, all species call.
     rho = 1.0d10
     temp = 1.0d0
     ye = 0.35d0     
     call nulibtable_range_species_range_energy(rho,temp,ye,eas_species_energy, &
          nulibtable_number_species,nulibtable_number_groups,nulibtable_number_easvariables)
     
     write(*,*) "rho,temp,ye:",rho,temp,ye
     write(*,*) "data below: species index, energy, eas variables"

     do ns=1,nulibtable_number_species
        do ng=1,nulibtable_number_groups
           write(*,"(i4,1P10E18.9)") ns,nulibtable_energies(ng), &
                (eas_species_energy(ns,ng,k),k=1,nulibtable_number_easvariables)
        enddo
     enddo
  endif

end program driver
