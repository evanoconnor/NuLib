!-*-f90-*-
program driver

  use eosmodule
  implicit none


  real*8 xrho,xye,xtemp,xtemp2
  real*8 xenr,xprs,xent,xcs2,xdedt,xmunu
  real*8 xdpderho,xdpdrhoe
  integer keytemp,keyerr

  ! for full eos call:
  real*8 xabar,xzbar,xmu_e,xmu_n,xmu_p,xmuhat
  real*8 xxa,xxh,xxn,xxp

  integer i

  keytemp = 2
  keyerr  = 0

  xrho = 10.0d0**1.474994d1
  xtemp = 0.7d0
  xye = 0.3d0

  call readtable("/Users/evanoc/research/EOS/HShen_2010.h5")

! keyerr --> error output; should be 0
! rfeps --> root finding relative accuracy, set around 1.0d-10
! keytemp: 0 -> coming in with eps
!          1 -> coming in with temperature
!          2 -> coming in with entropy

  ! short eos call
  call nuc_eos_short(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,&
       xdpderho,xdpdrhoe,xmunu,keytemp,keyerr,precision)

  write(6,*) "######################################"
  write(6,"(1P10E15.6)") xrho,xtemp,xye
  write(6,"(1P10E15.6)") xenr,xprs,xent,sqrt(xcs2)
  write(6,"(1P10E15.6)") xdedt,xdpdrhoe,xdpderho
  write(6,*) "######################################"

  ! full eos call
  
  write(6,*) "Full EOS: ############################"
  write(6,"(1P10E15.6)") xrho,xtemp,xye
  write(6,"(1P10E15.6)") xenr,xprs,xent,sqrt(xcs2)
  write(6,"(1P10E15.6)") xdedt,xdpdrhoe,xdpderho
  write(6,"(1P10E15.6)") xabar,xzbar
  write(6,"(1P10E15.6)") xxa,xxh,xxn,xxp
  write(6,"(1P10E15.6)") xmu_e,xmu_p,xmu_n,xmuhat
  write(6,*) "######################################"

  
  xtemp2 = 2.0d0*xtemp
  keytemp = 0

  call nuc_eos_full(xrho,xtemp2,xye,xenr,xprs,xent,xcs2,xdedt,&
       xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p,&
       xmuhat,keytemp,keyerr,precision)

  write(6,*) "Full EOS: ############################"
  write(6,"(1P10E15.6)") xrho,xtemp2,xtemp,xye
  write(6,"(1P10E15.6)") xenr,xprs,xent,sqrt(xcs2)
  write(6,"(1P10E15.6)") xdedt,xdpdrhoe,xdpderho
  write(6,"(1P10E15.6)") xabar,xzbar
  write(6,"(1P10E15.6)") xxa,xxh,xxn,xxp
  write(6,"(1P10E15.6)") xmu_e,xmu_p,xmu_n,xmuhat
  write(6,*) "######################################"



end program driver
