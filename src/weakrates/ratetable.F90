!-*-f90-*-
module class_ratetable

  implicit none
  private
  public :: RateTable, new_RateTable, weakrates_table

  ! constructors declaration
  interface new_RateTable
     module procedure new_RateTableDefault, new_RateTableFromFile
  end interface new_RateTable

  ! members
  type RateTable
     integer nrho,nt9,nnuc
     real*8 :: range_t9(2), range_lrhoye(2)
     real*8, pointer,dimension(:,:,:,:) :: rates  ! rates[species,T,rhoYe,rates+uf(indexed by nrate)]
     real*8, pointer,dimension(:,:,:,:,:,:) :: C ! Matrix of spline coefficients (see desc. below)
     real*8, pointer,dimension(:,:) :: nuclear_species ! nuclear_species[nucleus, (Q, A, Z)]
     real*8, pointer,dimension(:) :: t9dat
     real*8, pointer,dimension(:) :: rhoyedat
     integer, pointer,dimension(:,:) :: nucleus_index ! output array index for a given (A,Z) in tables
  end type RateTable

! methods
contains

!------------------------------------------------------------------------------------!

  function new_RateTableDefault(nnuc,nt9,nrho) result(this)
    !""" Default RateTable constructor """

    implicit none
    ! inputs
    integer, intent(in) :: nnuc, nt9, nrho

    ! return value
    type(RateTable) :: this

    this%nrho = nrho
    this%nnuc = nnuc
    this%nt9 = nt9
    allocate(this%nucleus_index(nnuc,nnuc))
    allocate(this%rates(nnuc,nt9/nrho,nrho,7))
    allocate(this%nuclear_species(nnuc,3))
    allocate(this%t9dat(nt9/nrho))
    allocate(this%rhoYedat(nrho))
    allocate(this%C(nnuc,7,2,nrho,nt9,4))

  end function new_RateTableDefault

!------------------------------------------------------------------------------------!

  function new_RateTableFromFile(directory,filename) result(this)
    ! """ RateTable constructor (from data file) """

    implicit none
    ! return value
    type(RateTable) :: this

    ! inputs
    character*200, intent(in) :: filename
    character*200 :: path
    character*(*) directory
    character lindex
    character*200 :: line,params
    integer i,j,dim,A,Z,nfile,nuc,max_A,max_Z
    real*8 :: t9,lrho,uf,lbetap,leps,lnu,lbetam,lpos,lanu
    real*8 :: lrho_prior,nucA,nucZ,nucQ
    logical continue_reading

    nuc = 0
    this%nrho = 0
    this%nt9 = 0
    dim = 1
    max_A=0
    max_Z=0
    allocate(this%nucleus_index(1000,1000))
    this%nucleus_index = 0
    path = trim(adjustl(directory))//trim(adjustl(filename))

    call print_reference(filename)

    !$OMP CRITICAL
    open(1,file=path,status='old')
    do
       read(1,'(A)',end=10) line
       read(line,*) lindex
       if (lindex.eq.'p') then
          continue_reading = .true.
          read(line(index(line(1:),"Q=")+2:index(line(1:),"Q=")+9),*) nucQ
          read(line(index(line(1:),"a=")+2:index(line(1:),"a=")+5),*) nucA
          read(line(index(line(1:),"z=")+2:index(line(1:),"z=")+4),*) nucZ
          A = int(nucA)
          Z = int(nucZ)+1
          max_A = max(max_A,A)
          max_Z = max(max_Z,Z)
          if(this%nucleus_index(A,Z).ne.0) then
             continue_reading = .false.
             cycle
          end if
          nuc = nuc + 1
          this%nrho = 0
          this%nt9 = 0
       end if
       if (index('0123456789',lindex).ne.0.and.continue_reading) then
          lrho_prior = lrho
          read(line,*) t9,lrho,uf,lbetap,leps,lnu,lbetam,lpos,lanu
          if (lrho.ne.lrho_prior) then
             this%nrho = this%nrho + 1
          end if
          this%nt9 = this%nt9 + 1
       end if
    end do

10  close(1)


    deallocate(this%nucleus_index)
    allocate(this%nucleus_index(max_A,max_A))
    allocate(this%rates(nuc,this%nt9/this%nrho,this%nrho,7))
    allocate(this%nuclear_species(nuc,3))
    allocate(this%t9dat(this%nt9/this%nrho))
    allocate(this%rhoYedat(this%nrho))

    nuc = 0
    lrho = 0.0d0
    this%nrho = 0
    this%nt9 = 0
    this%nucleus_index = 0

    open(1,file=path,status='old')
    do
       read(1,'(A)',end=20) line
       read(line,*) lindex
       if (lindex.eq.'p') then
          continue_reading = .true.
          read(line(index(line(1:),"Q=")+2:index(line(1:),"Q=")+9),*) nucQ
          read(line(index(line(1:),"a=")+2:index(line(1:),"a=")+5),*) nucA
          read(line(index(line(1:),"z=")+2:index(line(1:),"z=")+4),*) nucZ
          A = int(nucA)
          Z = int(nucZ)+1
          if(this%nucleus_index(A,Z).ne.0) then
             continue_reading = .false.
             cycle
          end if
          nuc = nuc + 1
          this%nuclear_species(nuc,1) = nucQ
          this%nuclear_species(nuc,2) = nucA
          this%nuclear_species(nuc,3) = nucZ
          this%nuclear_species(nuc,3) = this%nuclear_species(nuc,3)+1
          !^ currently reading in z of daughter which is
          !1 less than that of parent, hence the addition of one
          this%nucleus_index(A,Z) = nuc
          this%nrho = 0
          this%nt9 = 0
       end if
       if (index('0123456789',lindex).ne.0.and.continue_reading) then
          lrho_prior = lrho
          read(line,*) t9,lrho,uf,lbetap,leps,lnu,lbetam,lpos,lanu
          if (lrho.ne.lrho_prior) then
             this%nrho = this%nrho + 1
             this%nt9 = 0
          end if
          this%nt9 = this%nt9 + 1
          this%rates(nuc,this%nt9,this%nrho,1)=lbetap
          this%rates(nuc,this%nt9,this%nrho,2)=leps
          this%rates(nuc,this%nt9,this%nrho,3)=lnu
          this%rates(nuc,this%nt9,this%nrho,4)=lbetam
          this%rates(nuc,this%nt9,this%nrho,5)=lpos
          this%rates(nuc,this%nt9,this%nrho,6)=lanu
          this%rates(nuc,this%nt9,this%nrho,7)=uf
          this%t9dat(this%nt9)=t9
          this%rhoYedat(this%nrho)=lrho
       end if
    end do

20  close(1)

    !$OMP END CRITICAL

    this%range_t9(1)=minval(this%t9dat)
    this%range_t9(2)=maxval(this%t9dat)
    this%range_lrhoye(1)=minval(this%rhoyedat)
    this%range_lrhoye(2)=maxval(this%rhoyedat)

    !$OMP SINGLE
    write(*,"(A37,I4,A15,I2,A1,I2,A29)") &
         "      - Done reading weak rates for ",nuc,&
         " nuclei across ",this%nrho,"/",this%nt9,&
         " log10(rhoYe)/T9 grid points."
    write(*,"(A18,F5.2,A1,F6.2,A27,F5.2,A1,F5.2,A1)") &
         "      - T9(GK): [", this%range_t9(1),",",this%range_t9(2),&
         "] , log10(rhoYe (g/cm3)): [", this%range_lrhoye(1),",",&
         this%range_lrhoye(2),"]"
    !$OMP END SINGLE

    ! build array of interpolating spline coefficients
    this%nnuc = nuc
    allocate(this%C(this%nnuc,7,2,this%nrho,this%nt9,4))
    call monotonic_interp_2d(this)

    !$OMP PARALLEL FIRSTPRIVATE(this)
    !$OMP END PARALLEL
  end function new_RateTableFromFile

  !------------------------------------------------------------------------------------!

  subroutine monotonic_interp_2d(this)

    implicit none
    type(RateTable), intent(inout) :: this
    integer :: dim
    real*8, allocatable, dimension(:,:) :: Data
    INTEGER i,j,n,r,idxnuc,idxrate
    dim = 1

    allocate(Data(this%nt9,2))

    idxnuc = 1
    do idxnuc=1,this%nnuc
       do i=1,this%nrho
          do idxrate=1,7
             do j=1,this%nt9
                Data(j,1)=this%t9dat(j)
                Data(j,2)=this%rates(idxnuc,j,i,idxrate)
             end do
             call monotonic_interpolator(this,idxnuc,idxrate,dim,i,this%nt9,Data)
          end do
       end do
    end do

    return
  end subroutine monotonic_interp_2d

!------------------------------------------------------------------------------------!

  subroutine monotonic_interpolator(this,idxnuc,idxrate,dim,spl,N,Data)

    type(RateTable), intent(inout) :: this
    integer, intent(in) :: idxnuc, idxrate

    integer N      ! Number of data points (x,y) pairs
    real*8, dimension(N,2) :: Data
    real*8, dimension(N-1) :: S       ! Slope of secant passing through i+1 and i
    real*8, dimension(N-1) :: h       ! Size of interval, x_i+1 - x_i
    real*8, dimension(1:N) :: P     ! Slope of parabola passing through points i-1,i,i+1
    real*8, dimension(N) :: Dy         ! Value of 1st derivative at point i
    real*8 :: unity
    integer i, j, dim, spl
    unity = 1.0d0


    ! Build the array of slopes
    do i=1,N-1
       S(i)=(Data(i+1,2)-Data(i,2))/(Data(i+1,1)-Data(i,1))
    end do
    ! Build the array of intervals
    do i=1,N-1
       h(i)=Data(i+1,1)-Data(i,1)
    end do
    ! Build the array of parabolic tangents
    P(1)=S(1)*(1+h(1)/(h(1)+h(2)))-S(2)*h(1)/(h(1)+h(2))
    do i=2,N-1
       P(i)=(S(i-1)*h(i)+S(i)*h(i-1))/(h(i-1)+h(i))
    end do
    P(N)=S(N-1)*(1+h(N-1)/(h(N-1)+h(N-2)))-S(N-2)*h(N-1)/(h(N-1)+h(N-2))
    ! Build the array of 1st derivatives
    do i=1,N
       if (i==1) then
          if (P(1)*S(1) <= 0) then
             Dy(1)=0.0d0
          else if (abs(P(1)) > 2.0d0*abs(S(1))) then
             Dy(1)=2.0d0*S(1)
          else
             Dy(1)=P(1)
          end if
       else if (i==N) then
          if (P(N)*S(N-1) <= 0.0d0) then
             Dy(N)=0.0d0
          else if (abs(P(N)) > 2.0d0*abs(S(N-1))) then
             Dy(N)=2.0d0*S(N-1)
          else
             Dy(N)=P(N)
          end if
       else
          if (S(i-1)*S(i) <= 0.0d0) then
             Dy(i)=0
          else if (abs(P(i)) > 2.0d0*abs(S(i-1))) then
             Dy(i)=2.0d0*sign(unity,S(i-1))*min(abs(S(i-1)),abs(S(i)))
          else if (abs(P(i)) > 2.0d0*abs(S(i))) then
             Dy(i)=2.0d0*sign(unity,S(i-1))*min(abs(S(i-1)),abs(S(i)))
          else
             Dy(i)=P(i)
          end if
       end if
    end do

    ! Build the matrix of coefficients
    do i=1,N-1
       do j=1,4
          if (j==1) then
             this%C(idxnuc,idxrate,dim,spl,i,j)=(Dy(i)+Dy(i+1)-2.0d0*S(i))/(h(i)**2)
          else if (j==2) then
             this%C(idxnuc,idxrate,dim,spl,i,j)=(3.0d0*S(i)-2.0d0*Dy(i)-Dy(i+1))/h(i)
          else if (j==3) then
             this%C(idxnuc,idxrate,dim,spl,i,j)=Dy(i)
          else if (j==4) then
             this%C(idxnuc,idxrate,dim,spl,i,j)=Data(i,2)
          end if
       end do
    end do

    return

  end subroutine monotonic_interpolator

!------------------------------------------------------------------------------------!

  function weakrates_table(this,idxnuc,query_t9,query_lrhoye,idxrate) result(interp_val)

    type(RateTable), intent(inout) :: this

    real*8, allocatable, dimension(:,:) :: Data2d
    real*8 :: query_t9,query_lrhoye,value,interp_val
    integer i,j,counter,dim,idxnuc,idxrate
    dim = 1
    allocate(Data2d(this%nrho,2))


    do i=1,this%nrho
       do j=1,this%nt9
          if (query_t9 <= this%t9dat(j)) then
             if(j==1) then
                value = this%C(idxnuc,idxrate,dim,i,j,1)*(query_t9-this%t9dat(j))**3 +&
                     this%C(idxnuc,idxrate,dim,i,j,2)*(query_t9-this%t9dat(j))**2 +&
                     this%C(idxnuc,idxrate,dim,i,j,3)*(query_t9-this%t9dat(j)) +&
                     this%C(idxnuc,idxrate,dim,i,j,4)
                exit
             else
                value = this%C(idxnuc,idxrate,dim,i,j-1,1)*(query_t9-this%t9dat(j-1))**3 +&
                     this%C(idxnuc,idxrate,dim,i,j-1,2)*(query_t9-this%t9dat(j-1))**2 +&
                     this%C(idxnuc,idxrate,dim,i,j-1,3)*(query_t9-this%t9dat(j-1)) +&
                     this%C(idxnuc,idxrate,dim,i,j-1,4)
                exit
             end if
          else if (query_t9.gt.this%t9dat(j).and.j.eq.this%nt9) then
             stop "weak_rates.weakrates :: t > tmax"
          end if
       end do
       Data2d(i,1) = this%rhoYedat(i)
       Data2d(i,2) = value
    end do

    call monotonic_interpolator(this,idxnuc,idxrate,2,1,this%nrho,Data2d)
    interp_val = interpolant(this,idxnuc,idxrate,2,1,query_lrhoye)
    deallocate(Data2d)

    return

  end function weakrates_table

!------------------------------------------------------------------------------------!

  function interpolant(this,idxnuc,idxrate,dim,spl,query_lrhoye) result (interp_val)

    type(RateTable), intent(in) :: this
    integer, intent(in) :: idxnuc, idxrate
    real*8 :: interp_val

    real*8 :: query_lrhoye,result
    INTEGER i,counter,dim,spl


    do i=1,this%nrho
       if (query_lrhoye <= this%rhoYedat(i)) then
          if (i==1) then
             result = this%C(idxnuc,idxrate,dim,spl,i,1)*(query_lrhoye-this%rhoYedat(i))**3 +&
                  this%C(idxnuc,idxrate,dim,spl,i,2)*(query_lrhoye-this%rhoYedat(i))**2 +&
                  this%C(idxnuc,idxrate,dim,spl,i,3)*(query_lrhoye-this%rhoYedat(i)) +&
                  this%C(idxnuc,idxrate,dim,spl,i,4)
             exit
          else
             result = this%C(idxnuc,idxrate,dim,spl,i-1,1)*(query_lrhoye-this%rhoYedat(i-1))**3 +&
                  this%C(idxnuc,idxrate,dim,spl,i-1,2)*(query_lrhoye-this%rhoYedat(i-1))**2 +&
                  this%C(idxnuc,idxrate,dim,spl,i-1,3)*(query_lrhoye-this%rhoYedat(i-1)) +&
                  this%C(idxnuc,idxrate,dim,spl,i-1,4)
             exit
          end if
          !extrapolate to lrhoYe at most, consider adding extrapolate flag
       else if (query_lrhoye.gt.this%rhoYedat(i).and.i.eq.this%nrho) then
          stop "weak_rates.interpolant :: rho > rhomax"
       end if

    end do

    interp_val = result
    return

  end function interpolant

!------------------------------------------------------------------------------------!

  subroutine print_reference(filename)
    implicit none
    character*200, intent(in) :: filename
    !$OMP SINGLE
    select case (filename)
    case ("ravlicrates.dat")
       print *, "    Loading Ravlic et al. table. Make reference to: "
       print *, "    ------------------------------------------------------------------------------------------"
       print *, "    | iravlic | Ravlic, A., Giraud, S., Paar, N., Zegers, R. G. T (2024).                    |"
       print *, "    |         | Self-consistent microscopic calculations for electron captures on nuclei in  |"
       print *, "    |         | core-collapse supernovae                                                     |"
       print *, "    |         | arXiv:2412.00650v1                                                           |"
       print *, "    |         | https://arxiv.org/pdf/2412.00650                                             |"
       print *, "    ------------------------------------------------------------------------------------------"
    case ("suzukirates.dat")
       print *, "    Loading Suzuki et al. table, Make reference to: "
       print *, "    ------------------------------------------------------------------------------------------"
       print *, "    | isuzuki | Toshio Suzuki, Hiroshi Toki and Ken'ichi Nomoto (2016).                      |"
       print *, "    |         | ELECTRON-CAPTURE AND beta-DECAY RATES FOR sd-SHELL NUCLEI IN STELLAR         |"
       print *, "    |         | ENVIRONMENTS RELEVANT TO HIGH-DENSITY O–NE–MG CORES.                         |"
       print *, "    |         | Astrophys. J. 817, 163.                                                      |"
       print *, "    |         | https://doi.org/10.3847/0004-637x/817/2/163                                  |"
       print *, "    ------------------------------------------------------------------------------------------"
    case ("lmprates.dat")
       print *, "    Loading LMP table. Make reference to: "
       print *, "    ------------------------------------------------------------------------------------------"
       print *, "    | ilmp    | Langanke, K., & Mart\'{i}nez-Pinedo, G. (2000).                              |"
       print *, "    |         | Shell-model calculations of stellar weak interaction rates:                  |"
       print *, "    |         | II. Weak rates for nuclei in the mass range in supernovae environments.      |"
       print *, "    |         | Nuclear Physics A, 673(1-4), 481-508.                                        |"
       print *, "    |         | http://doi.org/10.1016/S0375-9474(00)00131-7                                 |"
       print *, "    ------------------------------------------------------------------------------------------"
    case ("lmshrates.dat")
       print *, "    Loading LMSH table. Make reference to: "
       print *, "    ------------------------------------------------------------------------------------------"
       print *, "    | ilmsh   | Langanke, K., & Mart\'{i}nez-Pinedo, G. (2003).                              |"
       print *, "    |         | Electron capture rates on nuclei and implications for stellar core collapse. |"
       print *, "    |         | Physical Review Letters 90, 241102.                                          |"
       print *, "    |         | http://prl.aps.org/abstract/PRL/v90/i24/e241102                              |"
       print *, "    ------------------------------------------------------------------------------------------"
    case ("odarates.dat")
       print *, "    Loading Oda et al. table. Make reference to: "
       print *, "    ------------------------------------------------------------------------------------------"
       print *, "    | ioda    | Oda, T., Hino, M., Muto, K., Takahara, M., & Sato, K. (1994).                |"
       print *, "    |         | Rate Tables for the Weak Processes of sd-Shell Nuclei in Stellar Matter.     |"
       print *, "    |         | Atomic Data and Nuclear Data Tables, 56(2), 231-403.                         |"
       print *, "    |         | http://doi.org/10.1006/adnd.1994.1007                                        |"
       print *, "    ------------------------------------------------------------------------------------------"
    case ("ffnrates.dat")
       print *, "    Loading FFN table. Make reference to: "
       print *, "    ------------------------------------------------------------------------------------------"
       print *, "    | iffn    | Fuller, G. M., Fowler, W. A., & Newman, M. J. (1982).                        |"
       print *, "    |         | Stellar weak interaction rates for intermediate-mass nuclei.                 |"
       print *, "    |         | II - A = 21 to A = 60. The Astrophysical Journal, 252, 715.                  |"
       print *, "    |         | http://doi.org/10.1086/159597                                                |"
       print *, "    ------------------------------------------------------------------------------------------"
    case default
       stop "No default"
    end select
    !$OMP END SINGLE

  end subroutine print_reference
!------------------------------------------------------------------------------------!
end module class_ratetable
