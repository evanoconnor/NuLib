!-*-f90-*-
#define NUM_TABLES 7
#define gScale 1.0d0
#define scale_all .false.
#define scale_diamond .false.
module class_ratelibrary

  use class_ratetable
  use class_rateapproximation
  
  implicit none
  private
  public :: RateLibrary, new_RateLibrary, in_table, return_weakrate

  ! constructors declaration
  interface new_RateLibrary
     module procedure new_RateLibraryDefault
  end interface new_RateLibrary

  interface return_weakrate
     module procedure return_weakrate_dynamic_search, return_weakrate_from_table, return_weakrate_from_approx
  end interface return_weakrate

  type RateLibrary
     ! number of allocated tables
     integer :: ntables
     ! array of rate tables
     type(RateTable), dimension(:), pointer :: tables
     ! approximate weak rate object
     type(RateApprox) :: approx
     ! path to directory containing weakrate tables
     character*200 :: directory
     ! weak interaction rate data table file names
     character*200, dimension(NUM_TABLES) :: files_to_load
     ! table priorities
     integer, dimension(NUM_TABLES+1) :: priority
     ! file ordering
     integer, dimension(NUM_TABLES) :: ifiles
     ! path to eos file
     character*200 :: eos_path
  end type RateLibrary
  
  ! class instance (singleton)
  type(RateLibrary), save, target :: this

  ! members
  type(RateTable), dimension(NUM_TABLES), save, target :: ratetables
  !integer, dimension(500,500),save :: gindex

  !$OMP THREADPRIVATE(this,ratetables)
  
  ! methods
contains

!------------------------------------------------------------------------------------!
  
  function new_RateLibraryDefault(parameters) result(library)
    !""" Default RateLibrary constructor """
    
    implicit none
    type(RateLibrary), pointer :: library
    character*200 :: parameters
    integer :: idxfiles,nfiles,i,j,k
    character*200 :: filename

    !$OMP CRITICAL
    call weakrate_inputparser(parameters,this)
    !$OMP END CRITICAL
    call print_reference
    
    ! construct rate table objects
    nfiles = 4
    idxfiles = 0
    this%ntables = 0    
    this%ifiles = 0
    do i=1,nfiles
       if(this%priority(i).gt.0)idxfiles = idxfiles + 1
    end do
    nfiles = idxfiles

    do i=1,size(this%ifiles)       ! order to load file
       do j=1,size(this%ifiles)    ! files in order
          if(this%priority(j).eq.i)then
             this%ifiles(i)=j
          end if
       end do
    end do
    if(Sum(this%ifiles).ne.0)then
       do i=1,size(this%ifiles)
          if(this%ifiles(i).eq.0) cycle
          filename=this%files_to_load(this%ifiles(i))          
          ratetables(i) = new_RateTable(this%directory,filename)
          this%ntables = this%ntables + 1
       enddo
    endif

    ! construct rate approximation object
    call print_approx_reference
    this%approx = new_RateApprox()

    
    library => this
    this%tables => ratetables  
    
    !$OMP PARALLEL COPYIN(this,ratetables)
    !$OMP END PARALLEL

    return
  end function new_RateLibraryDefault

  !------

  function scaling_factor(A,Z) result(scale)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!added by rachel!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer j, myIndex, A, Z
    !number of nuclei in the high sensitivity region, 30, 74 or 315
    integer, parameter :: nNuclei = 74
    integer hsA(nNuclei), hsZ(nNuclei)

!!! YOU'RE SCALING ALL RATES AT THE MOMENT
!!! GO CHANGE THE SCALING IF THIS IS NOT CORRECT
    
    real*8 :: scale

    !specific rate scaling happens here
    if (scale_all) then
       ! scale everything
       scale = gScale
       return
    else if (scale_diamond) then
       !nucleus is in the diamond, scale the rate
       hsZ = (/ 26,26,27,27,27,27,28,28,28,28,28,28,29,29,29,29,29,29,29,29,30,30,30,30,30,30,30,30,30,30,31,31,31,31,31,31,31,31,31,31,31,32,32,32,32,32,32,32,32,32,32,33,33,33,33,33,33,33,33,33,33,33,34,34,34,34,34,34,35,35,35,35,36,36 /)

       hsA = (/ 75,76,75,76,77,78,75,76,77,78,79,80,75,76,77,78,79,80,81,82,75,76,77,78,79,80,81,82,83,84,75,76,77,78,79,80,81,82,83,84,85,76,77,78,79,80,81,82,83,84,85,75,76,77,78,79,80,81,82,83,84,85,80,81,82,83,84,85,82,83,84,85,84,85 /)

       myIndex = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do j=1,nNuclei
          if(hsA(j).eq.A.and.hsZ(j).eq.Z) then
             !got a match, rate needs to be scaled
             !print *, "flag worked"
             myIndex = -1
             exit !exit the do loop and continue with the next lines
          else
             !no match, carry on with the j loop searching for a match
          end if
       end do

       if (myIndex.lt.0) then
          scale = gScale
          return 
       else
          scale = 1.0d0
          return
       end if

    else
       ! not scaling anything
       scale = 1.0d0
       return
    endif
    
    
  end function scaling_factor
!------------------------------------------------------------------------------------!

  function return_weakrate_from_table(this,A,Z,query_t9,query_lrhoye,idxtable,idxrate) result(rate) 
    
    implicit none
    type(RateLibrary) :: this
    integer :: A, Z, idxtable, idxrate
    real*8 :: query_t9, query_lrhoye
    real*8 :: rate

    rate = scaling_factor(A,Z) * 10.0d0**(weakrates_table(ratetables(idxtable),ratetables(idxtable)%nucleus_index(A,Z),query_t9,query_lrhoye,idxrate))
    return
    
  end function return_weakrate_from_table
  
!------------------------------------------------------------------------------------!

  function return_weakrate_from_approx(idxrate,xtemp,xq,xmue,A,Z) result(rate) 
    
    implicit none
    type(RateLibrary) :: this
    integer :: idxrate, A, Z
    real*8 :: xtemp, xmue, xq
    real*8 :: rate
    
    rate = scaling_factor(A,Z) * weakrates_approx(idxrate,xtemp,xq,xmue)
    return
    
  end function return_weakrate_from_approx
  
!------------------------------------------------------------------------------------!

  function return_weakrate_dynamic_search(this,A,Z,xrho,xtemp,xye,xmue,idxrate) result(rate) 

    implicit none
    type(RateLibrary) :: this
    real*8 :: xrho, xtemp, xye, xmue, rate
    integer :: A, Z, idxrate, idxtable
    real*8 :: lrhoye, t9, q    
    real*8, parameter :: kelvin_to_mev = 8.6173423d-11 !one K is # MeV

    ! convert to units of table grid
    lrhoye = log10(xrho*xye)
    t9 = xtemp/kelvin_to_mev*1.0d-9
    ! determine which table should be used for a given rate, if any
    idxtable = in_table(this,A,Z,lrhoye,t9)    
    if (idxtable.eq.0) then
       ! use approx if no table contains a rate for (A,Z) at the req. point
       if (idxrate.eq.2.or.idxrate.eq.3)then
          q = return_hempel_qec(A,Z,Z-1)
          rate = scaling_factor(A,Z) * weakrates_approx(idxrate-2,xtemp,q,xmue) ! xmue should be mu_e-m_e
          return
       else
          stop "RateLibrary Error: approximate rates only exist for electron capture and neutrino e-loss"
       endif
    endif
    ! interpolate correct rate table - defined by the priority hierarchy set in parameters
    rate = scaling_factor(A,Z) * 10.0d0**(weakrates_table(ratetables(idxtable),ratetables(idxtable)%nucleus_index(A,Z),t9,lrhoye,idxrate))
    return

  end function return_weakrate_dynamic_search
  
!------------------------------------------------------------------------------------!

  function in_table(this,A,Z,lrhoye,t9) result (idxtable)
    implicit none
    type(RateLibrary) :: this
    integer :: A, Z
    real*8 :: lrhoye, t9
    integer :: idxtable, i
    real*8 :: min, max

!    print *, "hi"
    
    idxtable = 0
    min = 0.0d0
    max = 0.0d0
    ! first check if the nucleus is in a table
    do i=1,this%ntables
       if (A.gt.size(ratetables(i)%nucleus_index,dim=1).or.Z.gt.size(ratetables(i)%nucleus_index,dim=2)) then
!          print *, "if"
          idxtable = 0
       else
!          print *, "else"
!          print *, ratetables(i)%nucleus_index(A,Z)
          idxtable = ratetables(i)%nucleus_index(A,Z)
       endif
       if (idxtable.ne.0)then
          idxtable = i
          exit
       endif       
    end do

!    print *, "hi2"
    
    ! return if the requested nucleus i not found in a table
    if (idxtable.eq.0)then
       return
    endif

    ! check if the requested point is in the grid    
    if (t9.ge.ratetables(idxtable)%range_t9(1).and.t9.le.ratetables(idxtable)%range_t9(2)&
         .and.&
         lrhoye.ge.ratetables(idxtable)%range_lrhoye(1).and.lrhoye.le.ratetables(idxtable)%range_lrhoye(2)) then
    else
       idxtable = 0
    endif

!    print *, "hi3"
    
    return
    
  end function in_table
  
  subroutine weakrate_inputparser(fn,library)

    use inputparser
    
    implicit none
    character*(*) fn
    type(RateLibrary) :: library

    call get_string_parameter(fn,'directory',library%directory)
    call get_string_parameter(fn,'lmp_rates',library%files_to_load(1))
    call get_string_parameter(fn,'lmsh_rates',library%files_to_load(2))
    call get_string_parameter(fn,'oda_rates',library%files_to_load(3))
    call get_string_parameter(fn,'ffn_rates',library%files_to_load(4))
    call get_string_parameter(fn,'pruet_rates1',library%files_to_load(5))
    call get_string_parameter(fn,'pruet_rates2',library%files_to_load(6))
    call get_string_parameter(fn,'pruet_rates3',library%files_to_load(7))
!    call get_string_parameter(fn,'calc_rates',library%files_to_load(8))

    call get_integer_parameter(fn,'ilmp',library%priority(1))
    call get_integer_parameter(fn,'ilmsh',library%priority(2))
    call get_integer_parameter(fn,'ioda',library%priority(3))
    call get_integer_parameter(fn,'iffn',library%priority(4))
    call get_integer_parameter(fn,'ipruet1',library%priority(5))
    call get_integer_parameter(fn,'ipruet2',library%priority(6))
    call get_integer_parameter(fn,'ipruet3',library%priority(7))
!    call get_integer_parameter(fn,'icalc',library%priority(8))
    call get_integer_parameter(fn,'iapprox',library%priority(9))
    call get_string_parameter(fn,'eos_table_name',library%eos_path)
    
  end subroutine weakrate_inputparser

!------------------------------------------------------------------------------------!

  subroutine print_reference
    !$OMP SINGLE
    print *, 
    print *, "Loading weak rate library. Make reference to:"    
    print *, "------------------------------------------------------------------------------------"
    print *, "| Sullivan, C., O'Connor, E., Zegers, R. G. T., Grubb, T., & Austin, S. M. (2015). |"
    print *, "| The Sensitivity of Core-Collapse Supernovae to Nuclear Electron Capture.         |"
    print *, "| http://arxiv.org/abs/1508.07348                                                  |"
    print *, "| Contact: Chris Sullivan <sullivan@nscl.msu.edu>                                  |"
    print *, "------------------------------------------------------------------------------------"
    !$OMP END SINGLE
  end subroutine print_reference
  
!------------------------------------------------------------------------------------!
end module class_ratelibrary

