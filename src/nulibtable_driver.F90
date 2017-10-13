!-*-f90-*-
program driver
  
  use nulibtable
  use nulibtable_interface
  implicit none

  character(len=1024) :: filename = "testertable.h5"

  real*8, allocatable :: eas(:)
  real*8, allocatable :: eas_energy(:,:)
  real*8, allocatable :: eas_species_energy(:,:,:)
  real*8, allocatable :: ies_energy_energy(:,:,:)
  real*8, allocatable :: ies_species_energy_energy(:,:,:,:)
  real*8, allocatable :: epannihil_energy_energy(:,:,:)
  real*8, allocatable :: epannihil_species_energy_energy(:,:,:,:)

  integer i,j,k,ns,ng,ngprime
  real*8 rho,temp,ye,mue,eta

  call nulibtable_reader(filename,include_Ielectron=.false.,include_epannihil_kernels=.false.,include_scattering_delta=.false.)

  allocate(eas(nulibtable_number_easvariables))
  allocate(eas_energy(nulibtable_number_groups,nulibtable_number_easvariables))
  allocate(eas_species_energy(nulibtable_number_species, &
       nulibtable_number_groups,nulibtable_number_easvariables))
  allocate(ies_energy_energy(nulibtable_number_groups,nulibtable_number_groups,2))
  allocate(ies_species_energy_energy(nulibtable_number_species, &
       nulibtable_number_groups,nulibtable_number_groups,2))
  allocate(epannihil_energy_energy(nulibtable_number_groups,nulibtable_number_groups,4))
  allocate(epannihil_species_energy_energy(nulibtable_number_species, &
       nulibtable_number_groups,nulibtable_number_groups,4))

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
  write(*,*) "Units are: inelastic electron scattering kernels, cm^3/s"
  write(*,*) "Units are: e- + e+ <-> \nu + \bar{\nu} kernels, cm^3/s"
  
  if (.true.) then
     !example of single energy single species call, here I loop over rho
     write(*,*) "rho,temp,ye, eas variables"
     do i=1,100
        eas(:) = 0.0d0
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

  write(*,*) 
  write(*,*) 

  if (.true.) then
     !example of an inelastic kernel call, single species call.
     temp = 1.5d0
     mue = 10.0d0
     eta = mue/temp
     ns = 1

     call nulibtable_inelastic_single_species_range_energy(temp,eta,ns, &
          ies_energy_energy,nulibtable_number_groups,nulibtable_number_groups,2)
     
     write(*,*) "temp,eta:",temp,eta
     write(*,*) "data below: energy_in, energy_out, ies Phi_0, ies Phi_1"
     write(*,*) "below is just a sample of the kernels (every fourth energy)"

     do ng=1,nulibtable_number_groups,4
        do ngprime=1,nulibtable_number_groups,4
           write(*,"(1P10E18.9)") nulibtable_energies(ng),nulibtable_energies(ngprime), &
                ies_energy_energy(ng,ngprime,1),ies_energy_energy(ng,ngprime,2)
        enddo
     enddo

     if (nulibtable_number_groups.ge.10) then     
        write(*,*) "detail balence symmetry check, the following sets of two number should be the same"
        write(*,*) exp(-(nulibtable_energies(7)-nulibtable_energies(3))/temp)* &
             ies_energy_energy(7,3,1),ies_energy_energy(3,7,1)
        write(*,*) exp(-(nulibtable_energies(8)-nulibtable_energies(1))/temp)* &
            ies_energy_energy(8,1,2),ies_energy_energy(1,8,2)
     else
        write(*,*) "detail balence symmetry check, the following sets of two number should be the same"
        write(*,*) exp(-(nulibtable_energies(2)-nulibtable_energies(1))/temp)* &
            ies_energy_energy(2,1,1),ies_energy_energy(1,2,1)
     endif
  endif

  write(*,*) 
  write(*,*) 

  if (.true.) then
     !example of an inelastic kernel call, all species call.
     temp = 1.5d0
     mue = 10.0d0
     eta = mue/temp

     call nulibtable_inelastic_range_species_range_energy(temp,eta, &
          ies_species_energy_energy,nulibtable_number_species, &
          nulibtable_number_groups,nulibtable_number_groups,2)
     
     write(*,*) "temp,eta:",temp,eta
     write(*,*) "data below: species, energy_in, energy_out, ies Phi_0, ies Phi_1"
     write(*,*) "below is just a sample of the kernels (every fourth energy)"

     do ns=1,nulibtable_number_species
        do ng=1,nulibtable_number_groups,4
           do ngprime=1,nulibtable_number_groups,4
              write(*,"(i4,1P10E18.9)") ns,nulibtable_energies(ng),nulibtable_energies(ngprime), &
                   ies_species_energy_energy(ns,ng,ngprime,1),ies_species_energy_energy(ns,ng,ngprime,2)
           enddo
        enddo
     enddo

     if (nulibtable_number_groups.ge.10) then     
        write(*,*) "detail balence symmetry check, the following sets of two number should be the same"
        write(*,*) exp(-(nulibtable_energies(7)-nulibtable_energies(3))/temp)* &
            ies_species_energy_energy(1,7,3,1),ies_species_energy_energy(1,3,7,1)
        write(*,*) exp(-(nulibtable_energies(8)-nulibtable_energies(1))/temp)* &
            ies_species_energy_energy(2,8,1,1),ies_species_energy_energy(2,1,8,1)
        write(*,*) exp(-(nulibtable_energies(10)-nulibtable_energies(9))/temp)* &
            ies_species_energy_energy(3,10,9,1),ies_species_energy_energy(3,9,10,1)
     else
        write(*,*) "detail balence symmetry check, the following sets of two number should be the same"
        write(*,*) exp(-(nulibtable_energies(2)-nulibtable_energies(1))/temp)* &
            ies_species_energy_energy(1,2,1,1),ies_species_energy_energy(1,1,2,1)
     endif
  endif

  write(*,*) 
  write(*,*) 

  if (.true.) then
     !example of an ep-annihilation kernel call, single species call.
     temp = 1.5d0
     mue = 10.0d0
     eta = mue/temp
     ns = 3

     call nulibtable_epannihil_single_species_range_energy(temp,eta,ns, &
          epannihil_energy_energy,nulibtable_number_groups,nulibtable_number_groups,4)
     
     write(*,*) "temp,eta:",temp,eta
     write(*,*) "data below: energy_in, energy_out, epannihil production Phi_0, ", &
          "epannihil annihilation Phi_0, epannihil production Phi_1, epannihil annihilation Phi_1"
     write(*,*) "below is just a sample of the kernels (every fourth energy)"

     do ng=1,nulibtable_number_groups,4
        do ngprime=1,nulibtable_number_groups,4
           write(*,"(1P10E18.9)") nulibtable_energies(ng),nulibtable_energies(ngprime), &
                epannihil_energy_energy(ng,ngprime,1),epannihil_energy_energy(ng,ngprime,2), &
                epannihil_energy_energy(ng,ngprime,3),epannihil_energy_energy(ng,ngprime,4)
        enddo
     enddo

     if (nulibtable_number_groups.ge.7) then
        write(*,*) "in/out epannihil symmetry check, the following two number should be the same (assuming nux = nuxbar)"
        write(*,*) epannihil_energy_energy(3,7,1),epannihil_energy_energy(7,3,1)
        write(*,*) epannihil_energy_energy(3,7,2),epannihil_energy_energy(7,3,2)
        write(*,*) "detailed balance epannihil symmetry check, the following two number should be the same"
        write(*,*) epannihil_energy_energy(3,7,1), &
             exp(-(nulibtable_energies(3)+nulibtable_energies(7))/temp)* &
             epannihil_energy_energy(3,7,2)
        write(*,*) epannihil_energy_energy(1,8,1), &
             exp(-(nulibtable_energies(1)+nulibtable_energies(8))/temp)* &
             epannihil_energy_energy(1,8,2)
     else
        write(*,*) "in/out epannihil symmetry check, the following two number should be the same (assuming nux = nuxbar)"
        write(*,*) epannihil_energy_energy(1,2,1),epannihil_energy_energy(2,1,1)
        write(*,*) "detailed balance epannihil symmetry check, the following two number should be the same"
        write(*,*) epannihil_energy_energy(1,2,1), &
             exp(-(nulibtable_energies(1)+nulibtable_energies(2))/temp)* &
             epannihil_energy_energy(1,2,2)
     endif

  endif

  write(*,*) 
  write(*,*) 

  if (.true.) then
     !example of an ep-annihilation kernel call, all species call.
     temp = 1.5d0
     mue = 10.0d0
     eta = mue/temp

     call nulibtable_epannihil_range_species_range_energy(temp,eta, &
          epannihil_species_energy_energy,nulibtable_number_species, &
          nulibtable_number_groups,nulibtable_number_groups,4)
     
     write(*,*) "temp,eta:",temp,eta
     write(*,*) "data below: species, energy_in, energy_out, epannihil annihilation Phi_0,", &
          "epannihil production Phi_1, epannihil annihilation Phi_1"
     write(*,*) "below is just a sample of the kernels (every fourth energy)"

     do ns=1,nulibtable_number_species
        do ng=1,nulibtable_number_groups,4
           do ngprime=1,nulibtable_number_groups,4
              write(*,"(i4,1P10E18.9)") ns,nulibtable_energies(ng),nulibtable_energies(ngprime), &
                   epannihil_species_energy_energy(ns,ng,ngprime,1),epannihil_species_energy_energy(ns,ng,ngprime,2), &
                   epannihil_species_energy_energy(ns,ng,ngprime,3),epannihil_species_energy_energy(ns,ng,ngprime,4)
           enddo
        enddo
     enddo

     if (nulibtable_number_groups.ge.7) then
        write(*,*) "in/out epannihil symmetry check, the following two number should be the same"
        write(*,*) epannihil_species_energy_energy(1,3,7,1),epannihil_species_energy_energy(2,7,3,1)
        write(*,*) epannihil_species_energy_energy(1,3,7,3),epannihil_species_energy_energy(2,7,3,3)
        write(*,*) epannihil_species_energy_energy(1,1,8,2),epannihil_species_energy_energy(2,8,1,2)

        write(*,*) "detailed balance epannihil symmetry check, the following two number should be the same"
        write(*,*) epannihil_species_energy_energy(1,3,7,1), &
             exp(-(nulibtable_energies(3)+nulibtable_energies(7))/temp)* &
             epannihil_species_energy_energy(1,3,7,2)
        write(*,*) epannihil_species_energy_energy(1,1,8,1), &
             exp(-(nulibtable_energies(1)+nulibtable_energies(8))/temp)* &
             epannihil_species_energy_energy(1,1,8,2)
     else
        write(*,*) "in/out epannihil symmetry check, the following two number should be the same"
        write(*,*) epannihil_species_energy_energy(1,1,2,1),epannihil_species_energy_energy(2,2,1,1)
        write(*,*) "detailed balance epannihil symmetry check, the following two number should be the same"
        write(*,*) epannihil_species_energy_energy(1,1,2,1), &
             exp(-(nulibtable_energies(1)+nulibtable_energies(2))/temp)* &
             epannihil_species_energy_energy(1,1,2,2)
     endif     
   

  endif

end program driver
