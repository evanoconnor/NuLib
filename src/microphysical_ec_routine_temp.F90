!-*-f90-*- 
    subroutine microphysical_electron_capture()!ns,temp_emissivity,eos_variables)      
      use nulib, only : GLQ_n16_roots,GLQ_n16_weights
      implicit none

      real*8 :: avgE(2),leps,lnu,lbetap,uf
      integer i
      
      ! function declaration
      real*8 :: weak_rate_neutrino_spectra

!      do i=1,nuc
         ! call weakrate_interp(eos_variables,rate_of_interest)
         ! These four variables will be set after interpolation, setting manually for now to test
         leps = rates(nuc,1,1,1)    
         lnu = rates(nuc,1,1,2)
         lbetap = rates(nuc,1,1,3)
         uf = rates(nuc,1,1,7)
         avgE(1) = 10.0d0**(lnu)/(10.0d0**(leps) + 10.0d0**(lbetap))
         
         ! Find q_perT parameter in weak_rate_neutrino_spectra that solves avgE-avgE(q)=0
         avgE(2) = 0.0d0
         do i=1,16
            avgE(2) = avgE(2) + &
                 GLQ_n16_weights(i)*exp(GLQ_n16_roots(i))*weak_rate_neutrino_spectra(GLQ_n16_roots(i), &
                 t9dat(5),nucspec(1,1),rates(1,5,1,7))
         end do
         write(*,*) avgE(2)
         

!      end do
         
       contains

         function weak_rate_neutrino_spectra(nu_energy_perT,temp,q_perT,uf_perT) result(nu_spectra)
           include 'constants.inc'
           
           !inputs
           real*8, intent(in) :: nu_energy_perT
           real*8, intent(in) :: temp
           real*8, intent(in) :: q_perT
           real*8, intent(in) :: uf_perT
           
           !output
           real*8 ::  nu_spectra
           nu_spectra = ((kelvin_to_mev*temp)**4.0d0)*(nu_energy_perT**2.0d0)*((nu_energy_perT-q_perT)**2.0d0) &
                /(1+exp(nu_energy_perT-q_perT-uf_perT))
         end function weak_rate_neutrino_spectra
         
    end subroutine microphysical_electron_capture
