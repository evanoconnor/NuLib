!-*-f90-*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! READ THE EOS TABLE INTO MEMORY, SET THE REFERENCE MASS !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_eos_table(local_eos_filename)
  use nulib
  implicit none
  character*200, intent(IN) :: local_eos_filename !ignored if Helmholtz

#if HELMHOLTZ_EOS
  include 'other_eos/helmholtz/vector_eos.dek'
  call read_helm_table
  m_ref = m_amu
#elif NUCLEI_HEMPEL
  call readtable(local_eos_filename)
  m_ref = m_amu !for SFHo_EOS (Hempel)
#else
  call readtable(local_eos_filename)
  m_ref = m_n !for LS220 
#endif
end subroutine read_eos_table


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! CALL THE EOS TO FILL IN MISSING EOS VARIABLES !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_eos_variables(eos_variables)
  use nulib
  implicit none
  real*8, intent(inout) :: eos_variables(total_eos_variables)
  
#if HELMHOLTZ_EOS
  include 'other_eos/helmholtz/vector_eos.dek'

  !for the helmholtz EOS, only three species
  integer, parameter :: HELM_ionmax = 3
  real*8, parameter :: tiny = 1.0d-3
  real*8 HELM_xmass(HELM_ionmax),HELM_aion(HELM_ionmax),HELM_zion(HELM_ionmax),HELM_abar,HELM_zbar

  ! set the mass fractions, z's and a's of the composition, there is a
  ! lot of manual stuff here depending on your application, the
  ! Helmholtz EOS is very general when it comes to compositions

  ! no heavy nuclei
  if(eos_variables(xhindex) .ne. 0.0d0) then
     write(*,*) "ERROR: the helmholtz EOS is not set up to deal with heavy nuclei."
     stop
  end if

  ! get the NSE mass fractions from rho,T,Ye
  if (eos_variables(tempindex)<0.5) then
    eos_variables(xnindex) = 1d0 - 2d0*eos_variables(yeindex)
    eos_variables(xpindex) = 0d0
    eos_variables(xaindex) = 2d0*eos_variables(yeindex)
  else
    call nse_mass_fractions(eos_variables(rhoindex) , &
                            eos_variables(tempindex), &
                            eos_variables(yeindex)  , &
                            eos_variables(xnindex)  , &
                            eos_variables(xpindex)  , &
                            eos_variables(xaindex)  )
  end if

  ! make sure mass fractions add up
  if(abs(eos_variables(xnindex) + eos_variables(xpindex) + eos_variables(xaindex) - 1.0) > tiny) then
     write(*,*) "ERROR: mass fractions don't add up."
     stop
  end if

  ! neutron, protons, alphas - for Helmholtz EOS
  HELM_xmass(1) = eos_variables(xnindex) ; HELM_aion(1)  = 1.0d0 ; HELM_zion(1)  = 0.0d0
  HELM_xmass(2) = eos_variables(xpindex) ; HELM_aion(2)  = 1.0d0 ; HELM_zion(2)  = 1.0d0
  HELM_xmass(3) = eos_variables(xaindex) ; HELM_aion(3)  = 4.0d0 ; HELM_zion(3)  = 2.0d0

  ! average atomic weight and charge, used in Helmholtz EOS, not the same abar and zbar of NuLib 
  HELM_abar   = 1.0d0/sum(HELM_xmass(1:HELM_ionmax)/HELM_aion(1:HELM_ionmax))
  HELM_zbar   = HELM_abar*sum(HELM_xmass(1:HELM_ionmax)*HELM_zion(1:HELM_ionmax)/HELM_aion(1:HELM_ionmax))

  !NuLib does not include neutrons, proton, and alpha in abar and
  !zbar, unlike the Helmholtz EOS
  !setting these to 1 since eos_variables(xhindex) = 0.0d0, if you
  !have heavy nuclei other than neutrons,protons, and alphas, you must
  !set these appropiately
  eos_variables(abarindex) = 1.0d0
  eos_variables(zbarindex) = 1.0d0

  !here zbar and abar are an
  !average for ALL species, so defines ye, not true for
  !stellarcollapse.org EOS
  if(abs(eos_variables(yeindex) - HELM_zbar/HELM_abar) > tiny) then
     write(*,*) "ERROR: ye values don't match"
     stop
  end if

  ! call the eos
  !$OMP CRITICAL
  temp_row(1) = eos_variables(tempindex)/kelvin_to_mev
  den_row(1)  = eos_variables(rhoindex)
  abar_row(1) = HELM_abar
  zbar_row(1) = HELM_zbar
  jlo_eos = 1 ; jhi_eos = 1
  call helmeos
  !$OMP END CRITICAL

  !set eos_variables
  eos_variables(energyindex) = etot_row(1)
  eos_variables(mueindex) = etaele_row(1)*eos_variables(tempindex) + 0.511d0 !add in electron rest mass

  !analytic mu's from EOSmaker on stellarcollapse.org, has correction
  !to mu_p for coulomb, following our convention, we add in the rest
  !mass difference of the neutron and proton into the chemical
  !potentials and into muhat.
  call mu_np(eos_variables(rhoindex),eos_variables(tempindex)/kelvin_to_mev, &
       eos_variables(xnindex),eos_variables(xpindex),eos_variables(yeindex), &
       eos_variables(munindex),eos_variables(mupindex))

  eos_variables(muhatindex) = eos_variables(munindex) - eos_variables(mupindex)

#else
  real*8 :: matter_prs,matter_ent,matter_cs2,matter_dedt,matter_dpderho,matter_dpdrhoe
  integer :: keytemp,keyerr
  real*8 :: precision = 1.0d-10
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
#endif
end subroutine set_eos_variables


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! NEWTON-RHAPSON SOLVE FOR MASS FRACTIONS !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nse_mass_fractions(rho,T,Ye,xn,xp,xa)
  use nulib
  implicit none
  integer, parameter :: max_count = 100
  real*8, parameter :: max_error = 1d-6
  real*8, intent(in)  :: rho !g/ccm
  real*8, intent(in)  :: T   !MeV
  real*8, intent(in)  :: Ye
  real*8, intent(out) :: xn
  real*8, intent(out) :: xp
  real*8, intent(out) :: xa

  integer :: count
  real*8 :: error

  real*8 :: Ja,Jb,Jc,Jd
  real*8 :: Jinv_a,Jinv_b,Jinv_c,Jinv_d
  real*8 :: det
  real*8 :: mc,cc

  ! eta is my notation for (mu-mc^2)/kT, where mu is the chemical potential including rest mass
  real*8 :: eta_n
  real*8 :: eta_p
  real*8 :: eta_alpha

  !set initial conditions. Tricky. set s.t. X_alpha==1 and ( Ye>0.5 ? (X_p==1) : (X_n==1) )
  count = 0
  error = 1d0
  if(Ye<0.5) then
    eta_n = log(rho/(2d0*m_n*mev_to_gram) * (m_n*T/(2.0d0*pi*hbarc_mevcm**2))**(-1.5d0))
    eta_p =  0.5d0 * (log(rho/(m_alpha*mev_to_gram) * (m_alpha*T/(2.0d0*pi*hbarc_mevcm**2))**(-1.5d0)) &
          - (2d0*m_n + 2d0*m_p - m_alpha)/T) - eta_n
  else
    eta_p = log(rho/(2d0*m_p*mev_to_gram) * (m_p*T/(2.0d0*pi*hbarc_mevcm**2))**(-1.5d0))
    eta_n =  0.5d0 * (log(rho/(m_alpha*mev_to_gram) * (m_alpha*T/(2.0d0*pi*hbarc_mevcm**2))**(-1.5d0)) &
          - (2d0*m_n + 2d0*m_p - m_alpha)/T) - eta_p
  end if
  eta_alpha = eta_a(eta_n,eta_p,T)

  !set the mass fractions
  xn = X_n(rho,T,eta_n               )
  xp = X_p(rho,T,eta_p               )
  xa = X_a(rho,T,eta_a(eta_n,eta_p,T))

  do while (count<max_count .and. error>max_error)
     !Jacobian is [[Ja,Jb][Jc,Jd]] (left is d/d(eta_n), right is d/d(eta_p), up is d(mass constraint), down is d(charge constraint) )
     Ja = ddetan_mass_constraint(eta_n,eta_p,rho,T)
     Jb = ddetap_mass_constraint(eta_n,eta_p,rho,T)
     Jc = ddetan_charge_constraint(eta_n,eta_p,rho,T)
     Jd = ddetap_charge_constraint(eta_n,eta_p,rho,T)
     
     !inverse Jacobian
     det = (Ja*Jd - Jb*Jc)
     Jinv_a =  1d0/det * Jd
     Jinv_b = -1d0/det * Jb
     Jinv_c = -1d0/det * Jc
     Jinv_d =  1d0/det * Ja

     !calculate new q's
     mc = mass_constraint(eta_n,eta_p,rho,T)
     cc = charge_constraint(eta_n,eta_p,rho,T,Ye)
     eta_n = eta_n - (Jinv_a*mc + Jinv_b*cc)
     eta_p = eta_p - (Jinv_c*mc + Jinv_d*cc)
     !eta_alpha = eta_a(eta_n,eta_p,T)

     !set the mass fractions
     xn = X_n(rho,T,eta_n               )
     xp = X_p(rho,T,eta_p               )
     xa = X_a(rho,T,eta_a(eta_n,eta_p,T))

     !calculate error
     error = max( abs(mass_constraint(eta_n,eta_p,rho,T)), abs(charge_constraint(eta_n,eta_p,rho,T,Ye)) )
     count = count+1

  end do

  !warn if we didn't converge
  if(count>=max_count) then
     write(*,*) "WARNING: NSE solver could not converge"
  end if

  !set the mass fractions
  xn = X_n(rho,T,eta_n               )
  xp = X_p(rho,T,eta_p               )
  xa = X_a(rho,T,eta_a(eta_n,eta_p,T))

  if((xn.ne.xn) .or. (xp.ne.xp) .or. (xa.ne.xa)) then
    write(*,*) "ERROR: mass fractions are NaN"
    write(*,*) "rho=",rho
    write(*,*) "T=",T
    write(*,*) "Ye=",Ye
    stop
  end if

  contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! CALCULATE NSE MASS FRACTION !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function nse_mass_frac(rho,T,eta,g,m)
      implicit none
      real*8 :: nse_mass_frac
      real*8, intent(in) :: eta !mu/kT
      real*8, intent(in) :: T   !MeV
      real*8, intent(in) :: rho !g/ccm
      real*8, intent(in) :: g   !partition function (2 for nucleons, 1 for alpha)
      real*8, intent(in) :: m   !species mass (MeV)
      nse_mass_frac = g*m*mev_to_gram/rho * exp(eta) * (m*T/(2.0d0*pi*hbarc_mevcm**2))**(1.5d0)
      return
    end function nse_mass_frac
    function X_n(rho,T,eta)
      implicit none
      real*8 :: X_n
      real*8, intent(in) :: eta !mu/kT
      real*8, intent(in) :: T   !MeV
      real*8, intent(in) :: rho !g/ccm
      X_n = nse_mass_frac(rho,T,eta,2d0,m_n)
      return
    end function X_n
    function X_p(rho,T,eta)
      implicit none
      real*8 :: X_p
      real*8, intent(in) :: eta !mu/kT
      real*8, intent(in) :: T   !MeV
      real*8, intent(in) :: rho !g/ccm
      X_p = nse_mass_frac(rho,T,eta,2d0,m_p)
      return
    end function X_p
    function X_a(rho,T,eta)
      implicit none
      real*8 :: X_a
      real*8, intent(in) :: eta !mu/kT
      real*8, intent(in) :: T   !MeV
      real*8, intent(in) :: rho !g/ccm
      X_a = nse_mass_frac(rho,T,eta,1d0,m_alpha)
      return
    end function X_a

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! CALCULATE ETA_ALPHA FROM THE OTHER ETA'S !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function eta_a(eta_n,eta_p,T)
      implicit none
      real*8 :: eta_a
      real*8, intent(in) :: eta_n
      real*8, intent(in) :: eta_p
      real*8, intent(in) :: T     !MeV
      eta_a = 2d0*eta_n + 2d0*eta_p + (2d0*m_n + 2d0*m_p - m_alpha)/T
      return
    end function eta_a
    
    !!!!!!!!!!!!!!!!!!!!!!!
    !!! MASS CONSTRAINT !!!
    !!!!!!!!!!!!!!!!!!!!!!!
    function mass_constraint(eta_n,eta_p,rho,T)
      implicit none
      real*8 :: mass_constraint
      real*8, intent(in) :: eta_n
      real*8, intent(in) :: eta_p
      real*8, intent(in) :: rho !g/ccm
      real*8, intent(in) :: T   !MeV
      mass_constraint = X_n(rho,T,eta_n               ) &
                      + X_p(rho,T,eta_p               ) &
                      + X_a(rho,T,eta_a(eta_n,eta_p,T)) &
                      - 1d0
    end function mass_constraint
    function ddetan_mass_constraint(eta_n,eta_p,rho,T)
      implicit none
      real*8 :: ddetan_mass_constraint
      real*8, intent(in) :: eta_n
      real*8, intent(in) :: eta_p
      real*8, intent(in) :: rho !g/ccm
      real*8, intent(in) :: T   !MeV
      ddetan_mass_constraint = X_n(rho,T,eta_n               ) &
                             + X_a(rho,T,eta_a(eta_n,eta_p,T)) * 2d0
    end function ddetan_mass_constraint
    function ddetap_mass_constraint(eta_n,eta_p,rho,T)
      implicit none
      real*8 :: ddetap_mass_constraint
      real*8, intent(in) :: eta_n
      real*8, intent(in) :: eta_p
      real*8, intent(in) :: rho !g/ccm
      real*8, intent(in) :: T   !MeV
      ddetap_mass_constraint = X_p(rho,T,eta_p               ) &
                             + X_a(rho,T,eta_a(eta_n,eta_p,T)) * 2d0
    end function ddetap_mass_constraint
    
    !!!!!!!!!!!!!!!!!!!!!!!!!
    !!! CHARGE CONSTRAINT !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!
    function charge_constraint(eta_n,eta_p,rho,T,Ye)
      implicit none
      real*8 :: charge_constraint
      real*8, intent(in) :: eta_n
      real*8, intent(in) :: eta_p
      real*8, intent(in) :: rho !g/ccm
      real*8, intent(in) :: T   !MeV
      real*8, intent(in) :: Ye
      charge_constraint = X_p(rho,T,eta_p               ) * 1d0/1d0 &
                        + X_a(rho,T,eta_a(eta_n,eta_p,T)) * 2d0/4d0 &
                        - Ye
    end function charge_constraint
    function ddetan_charge_constraint(eta_n,eta_p,rho,T)
      implicit none
      real*8 :: ddetan_charge_constraint
      real*8, intent(in) :: eta_n
      real*8, intent(in) :: eta_p
      real*8, intent(in) :: rho !g/ccm
      real*8, intent(in) :: T   !MeV
      ddetan_charge_constraint = X_a(rho,T,eta_a(eta_n,eta_p,T)) * 2d0/4d0 * 2d0
    end function ddetan_charge_constraint
    function ddetap_charge_constraint(eta_n,eta_p,rho,T)
      implicit none
      real*8 :: ddetap_charge_constraint
      real*8, intent(in) :: eta_n
      real*8, intent(in) :: eta_p
      real*8, intent(in) :: rho !g/ccm
      real*8, intent(in) :: T   !MeV
      ddetap_charge_constraint = X_p(rho,T,eta_p               ) * 1d0/1d0 &
                               + X_a(rho,T,eta_a(eta_n,eta_p,T)) * 2d0/4d0 * 2d0
    end function ddetap_charge_constraint


end subroutine nse_mass_fractions
