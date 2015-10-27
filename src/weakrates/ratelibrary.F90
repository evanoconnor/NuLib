!-*-f90-*- 
module class_ratelibrary

  use class_ratetable
  use class_rateapproximation
  
  implicit none
  private
  public :: RateLibrary, new_RateLibrary

  ! constructors declaration
  interface new_RateLibrary
     module procedure new_RateLibraryDefault
  end interface new_RateLibrary

  ! members
  type RateLibrary
     ! array of rate tables
     type(RateTable), allocatable, dimension(:) :: tables
     ! approximate weak rate object
     type(RateApprox) :: approx
     ! weak interaction rate data table file names
     character*200, dimension(5) :: files_to_load
     ! table priorities
     integer, dimension(5) :: file_priority
     ! file ordering
     integer, dimension(4) :: ifiles   
  end type RateLibrary

  
  ! methods
contains

!------------------------------------------------------------------------------------!
  
  function new_RateLibraryDefault(parameters) result(this)
    !""" Default RateLibrary constructor """
    
    implicit none
    type(RateLibrary) :: this
    integer :: idxfiles=0,nfiles=0,i=0,j=0
    character*200 :: filename
    character*200 :: parameters

    call weakrate_inputparser(parameters,this)
    
    nfiles = 4
    this%ifiles = 0
    do i=1,nfiles
       if(this%file_priority(i).gt.0)idxfiles = idxfiles + 1
    end do
    nfiles = idxfiles
    allocate(this%tables(nfiles))

    do i=1,nfiles       ! order to load file
       do j=1,nfiles    ! files in order
          if(this%file_priority(j).eq.i)then
             this%ifiles(i)=j
          end if
       end do
    end do
    if(Sum(this%ifiles).ne.0)then
       do i=1,nfiles
          if(this%ifiles(i).eq.0) cycle
          filename=this%files_to_load(this%ifiles(i))
          this%tables(i) = new_RateTable(filename)          
       enddo
    endif
    
  end function new_RateLibraryDefault
  
!------------------------------------------------------------------------------------!
  
  subroutine weakrate_inputparser(fn,library)

    use inputparser
    
    implicit none
    character*(*) fn
    type(RateLibrary) :: library

    call get_string_parameter(fn,'lmp_rates',library%files_to_load(1))
    call get_string_parameter(fn,'lmsh_rates',library%files_to_load(2))
    call get_string_parameter(fn,'oda_rates',library%files_to_load(3))
    call get_string_parameter(fn,'ffn_rates',library%files_to_load(4))
    call get_integer_parameter(fn,'ilmp',library%file_priority(1))
    call get_integer_parameter(fn,'ilmsh',library%file_priority(2))
    call get_integer_parameter(fn,'ioda',library%file_priority(3))
    call get_integer_parameter(fn,'iffn',library%file_priority(4))
    call get_integer_parameter(fn,'iapprox',library%file_priority(5))
    
  end subroutine weakrate_inputparser
!------------------------------------------------------------------------------------!
end module class_ratelibrary

