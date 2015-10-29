!-*-f90-*-
#define NUM_TABLES 4
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

  ! members
  type RateLibrary
     ! array of rate tables
     type(RateTable), dimension(:), pointer :: tables
     ! approximate weak rate object
     type(RateApprox) :: approx
     ! weak interaction rate data table file names
     character*200, dimension(NUM_TABLES) :: files_to_load
     ! table priorities
     integer, dimension(NUM_TABLES+1) :: priority
     ! file ordering
     integer, dimension(NUM_TABLES) :: ifiles   
  end type RateLibrary

  
  ! methods
contains

!------------------------------------------------------------------------------------!
  
  function new_RateLibraryDefault(parameters) result(this)
    !""" Default RateLibrary constructor """
    
    implicit none
    type(RateLibrary) :: this
    character*200 :: parameters
    integer :: idxfiles=0,nfiles=0,i=0,j=0
    character*200 :: filename
    
    call weakrate_inputparser(parameters,this)

    ! construct rate approximation object
    !if(this%priority(5).ne.0)then
    this%approx = new_RateApprox()
    !endif

    ! construct rate table objects
    nfiles = 4
    this%ifiles = 0
    do i=1,nfiles
       if(this%priority(i).gt.0)idxfiles = idxfiles + 1
    end do
    nfiles = idxfiles
    allocate(this%tables(nfiles))

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
          print *, ">> Loading table  = ",filename
          this%tables(i) = new_RateTable(filename)          
       enddo
    endif
    
  end function new_RateLibraryDefault
  
!------------------------------------------------------------------------------------!

  function return_weakrate_from_table(this,A,Z,query_t9,query_lrhoye,idxtable,idxrate) result(rate) 
    
    implicit none
    type(RateLibrary) :: this
    integer :: A, Z, idxtable, idxrate
    real*8 :: query_t9, query_lrhoye
    real*8 :: rate

    if (idxtable.eq.0) print*, "idxtable is 0, how did this happen"
    rate = weakrates_table(this%tables(idxtable),A,Z,query_t9,query_lrhoye,idxrate)
    return
    
  end function return_weakrate_from_table
  
!------------------------------------------------------------------------------------!

  function return_weakrate_from_approx(idxrate,xtemp,xq,xmue) result(rate) 
    
    implicit none
    type(RateLibrary) :: this
    integer :: idxrate
    real*8 :: xtemp, xmue, xq
    real*8 :: rate
    
    rate = weakrates_approx(idxrate,xtemp,xq,xmue)
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
          rate = weakrates_approx(idxrate-2,xtemp,q,xmue) ! xmue should be mu_e-m_e
          return
       else
          stop "RateLibrary Error: approximate rates only exist for electron capture and neutrino e-loss"
       endif
    endif
    ! interpolate correct rate table - defined by the priority hierarchy set in parameters
    rate = weakrates_table(this%tables(idxtable),A,Z,t9,lrhoye,idxrate)
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

    idxtable = 0
    min = 0.0d0
    max = 0.0d0
    ! first check if the nucleus is in a table
    do i=1,size(this%tables)
       print *, size(this%tables(i)%nucleus_index,dim=1), size(this%tables(i)%nucleus_index,dim=1)
       if (A.gt.size(this%tables(i)%nucleus_index,dim=1).or.Z.gt.size(this%tables(i)%nucleus_index,dim=2)) then
          idxtable = 0
       else
          idxtable = this%tables(i)%nucleus_index(A,Z)
       endif
       if (idxtable.ne.0)then
          idxtable = i
          exit
       endif       
    end do
    ! return if the requested nucleus i not found in a table
    if (idxtable.eq.0)then
       return
    endif

    ! check if the requested point is in the grid    
    if (t9.ge.this%tables(idxtable)%range_t9(1).and.t9.le.this%tables(idxtable)%range_t9(2)&
         .and.&
         lrhoye.ge.this%tables(idxtable)%range_lrhoye(1).and.lrhoye.le.this%tables(idxtable)%range_lrhoye(2)) then
    else
       idxtable = 0
    endif
    
    return
    
  end function in_table
  
  subroutine weakrate_inputparser(fn,library)

    use inputparser
    
    implicit none
    character*(*) fn
    type(RateLibrary) :: library

    call get_string_parameter(fn,'lmp_rates',library%files_to_load(1))
    call get_string_parameter(fn,'lmsh_rates',library%files_to_load(2))
    call get_string_parameter(fn,'oda_rates',library%files_to_load(3))
    call get_string_parameter(fn,'ffn_rates',library%files_to_load(4))
    call get_integer_parameter(fn,'ilmp',library%priority(1))
    call get_integer_parameter(fn,'ilmsh',library%priority(2))
    call get_integer_parameter(fn,'ioda',library%priority(3))
    call get_integer_parameter(fn,'iffn',library%priority(4))
    call get_integer_parameter(fn,'iapprox',library%priority(5))
    
  end subroutine weakrate_inputparser
!------------------------------------------------------------------------------------!
end module class_ratelibrary

