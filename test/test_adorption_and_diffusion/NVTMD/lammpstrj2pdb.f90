!ITEM: TIMESTEP
!0
!ITEM: NUMBER OF ATOMS
!816
!ITEM: BOX BOUNDS xy xz yz pp pp pp
!-1.9741981549897282e+01 3.0372448059976723e+01 -1.0016673770959233e+01
!-1.3891208734959575e+01 2.8628128422373553e+01 -1.0032159718961333e+01
!2.3768098780503374e-01 2.4811788012194945e+01 -1.4187558312586022e+01
!ITEM: BOX BOUNDS xy xz yz
!xlo_bound xhi_bound xy
!ylo_bound yhi_bound xz
!zlo_bound zhi_bound yz
!ITEM: ATOMS element xs ys zs 
!Zn 0.250224 0.125095 0.375126 
program lammpstrj2pdb
 implicit none
 integer             :: ierr,i,j,k
 integer             :: timestep = 0
 integer             :: n_atoms = 0
 real                :: vr(3,3),rv(3,3),cell_0(6)
 real :: xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz,lx,ly,lz,cosa,cosb,cosg
 real :: xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound
 real,allocatable    :: xcrystal(:,:),xcarte(:,:)
 character(len=2),allocatable::label(:)
 character(len=80)   :: line
 character(len=40)   :: keyword
 real,parameter :: pi=acos(-1.0)
 real,parameter :: radtodeg = 180.0/pi
 real,parameter :: degtorad = pi/180.0
 open(100,file="out.pdb")
 k=0
 read_lammpstrj:do
  read(5,'(a)',iostat=ierr)line
  if(ierr/=0) exit read_lammpstrj
  if(line(1:5)=="ITEM:")then
   read(line,'(5x,1x,a)') keyword
   select case (keyword)
    case('TIMESTEP')
     read(5,*) timestep
     k=k+1
    case('NUMBER OF ATOMS')
     read(5,*) n_atoms
     allocate(xcrystal(3,n_atoms))
     allocate(xcarte(3,n_atoms))
     allocate(label(n_atoms))
    case('BOX BOUNDS xy xz yz pp pp pp')
     read(5,*)xlo_bound,xhi_bound,xy
     read(5,*)ylo_bound,yhi_bound,xz
     read(5,*)zlo_bound,zhi_bound,yz
     xlo=xlo_bound-MIN(0.0,xy,xz,xy+xz)
     xhi=xhi_bound-MAX(0.0,xy,xz,xy+xz)
     ylo=ylo_bound-MIN(0.0,yz)
     yhi=yhi_bound-MAX(0.0,yz)
     zlo=zlo_bound
     zhi=zhi_bound
     lx=xhi-xlo
     ly=yhi-ylo
     lz=zhi-zlo
     cell_0(1)=lx
     cell_0(2)=sqrt(ly*ly+xy*xy)
     cell_0(3)=sqrt(lz*lz+xz*xz+yz*yz)
     cell_0(4)=radtodeg*acos((xy*xz+ly*yz)/(cell_0(2)*cell_0(3)))
     cell_0(5)=radtodeg*acos(xz/cell_0(3))
     cell_0(6)=radtodeg*acos(xy/cell_0(2))
     call cell(rv,vr,cell_0)
     !write(6,*)timestep,(cell_0(j),j=1,6),volume(rv)
    case('ATOMS element xs ys zs')
     write(100,'(a5,i5)')'MODEL',k
     write(100,'(a6,9(f14.7,1x))')'REMARK',xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz
     write(100,'(a6,1x,a,1x,i6,1x,a,i6)')'REMARK','# atoms',n_atoms,'# timestep',timestep
     write(100,'(a6,1x,a,1x,f14.7)')'REMARK','Volume:',volume(rv)
     write(100,'(a6,3f9.3,3f7.2)')'CRYST1',(cell_0(j),j=1,6)
     call cell(rv,vr,cell_0)
     do i=1,n_atoms
      read(5,*) label(i),(xcrystal(j,i),j=1,3) 
      ! crystallyne -> cartesians
      forall ( j=1:3 )
       xcarte(j,i)=rv(j,1)*xcrystal(1,i) + rv(j,2)*xcrystal(2,i) + rv(j,3)*xcrystal(3,i)
      end forall
      write(100,'(a6,i5,1x,a2,3x,a4,1x,i4,4x,3f8.3,2f6.2,10x,a2)') &
      'ATOM  ',i,label(i),'MOL ',0,(xcarte(j,i),j=1,3),0.0,0.0,label(i) 
     end do
     write(100,'(a6)')'ENDMDL'
     deallocate(xcrystal)
     deallocate(xcarte)
     deallocate(label)
    case default
     STOP 'Keyword ??'
   end select
  end if
 end do read_lammpstrj
contains
SUBROUTINE cell(rv,vr,cell_0)
 implicit none
 integer :: i,j
 real, intent(in)  :: cell_0(6)
 real, intent(out) :: rv(3,3),vr(3,3)
 real :: alp,bet
 real :: cosa,cosb,cosg
 real :: gam,sing
 real :: pi,DEGTORAD
 pi = ACOS(-1.0)
 DEGTORAD=pi/180.0
 IF(cell_0(4) == 90.0) THEN
   cosa = 0.0
 ELSE
   ALP=cell_0(4)*degtorad
   COSA=cos(ALP)
 ENDIF
 IF(cell_0(5) == 90.0) THEN
   cosb = 0.0
 ELSE
   bet = cell_0(5)*degtorad
   cosb = cos(bet)
 ENDIF
 IF(cell_0(6) == 90.0) then
   sing = 1.0
   cosg = 0.0
 ELSE
   gam = cell_0(6)*degtorad
   sing = sin(gam)
   cosg = cos(gam)
 ENDIF
 rv(1,1) = cell_0(1)
 rv(1,2) = cell_0(2)*cosg
 rv(1,3) = cell_0(3)*cosb
 rv(2,1) = 0.0
 rv(2,2) = cell_0(2)*sing
 rv(2,3) = cell_0(3)*(cosa - cosb*cosg)/sing
 rv(3,1) = 0.0
 rv(3,2) = 0.0
 rv(3,3) = sqrt( cell_0(3)*cell_0(3) - rv(1,3)*rv(1,3) - rv(2,3)*rv(2,3)) 
 call inverse(rv,vr,3)
! print*,'Cell:'
! WRITE(*,'(6F14.7)')( cell_0(j), j=1,6 )
! print*,'Box:'
! DO i=1,3
!    WRITE(*,'(F14.7,F14.7,F14.7)')( rv(i,j), j=1,3 )
! ENDDO
! WRITE(*,*)'----------------------------------------'
! WRITE(*,*)'bOX:'
! DO i=1,3
!    WRITE(*,'(F14.7,F14.7,F14.7)')( vr(i,j), j=1,3 )
! ENDDO
 RETURN
END SUBROUTINE cell
!
SUBROUTINE uncell(rv,cell_0)
  implicit none
  real,intent(out) :: cell_0(6)
  real,intent(in)  :: rv(3,3)
  integer  :: i,j
  real     :: temp(6)
  REAL :: radtodeg,PI
  PI=ACOS(-1.0)
  radtodeg=180.0/PI
!
  do i = 1,3
    temp(i) = 0.0
    do j = 1,3
      temp(i) = temp(i) + rv(j,i)**2
    enddo
    temp(i) = sqrt(temp(i))
  enddo
  cell_0(1) = abs(temp(1))
  cell_0(2) = abs(temp(2))
  cell_0(3) = abs(temp(3))
  do i = 1,3
    temp(3+i) = 0.0
  enddo
  do j = 1,3
    temp(4) = temp(4) + rv(j,2)*rv(j,3)
    temp(5) = temp(5) + rv(j,1)*rv(j,3)
    temp(6) = temp(6) + rv(j,1)*rv(j,2)
  enddo
  temp(4) = temp(4)/(temp(2)*temp(3))
  temp(5) = temp(5)/(temp(1)*temp(3))
  temp(6) = temp(6)/(temp(1)*temp(2))
  cell_0(4) = radtodeg*acos(temp(4))
  cell_0(5) = radtodeg*acos(temp(5))
  cell_0(6) = radtodeg*acos(temp(6))
!  Avoid round off errors for 90.0 and 120.0 degrees
  DO i=4,6
     if (abs(cell_0(i) - 90.0 ).lt.0.00001) cell_0(i) = 90.0
     if (abs(cell_0(i) - 120.0).lt.0.00001) cell_0(i) = 120.0
  ENDDO
!
  return
end subroutine uncell
!
subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
 implicit none 
 integer n
 real a(n,n), c(n,n)
 real L(n,n), U(n,n), b(n), d(n), x(n)
 real coeff
 integer i, j, k
! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
 L=0.0
 U=0.0
 b=0.0
! step 1: forward elimination
 do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
 end do
! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
 do i=1,n
  L(i,i) = 1.0
 end do
! U matrix is the upper triangular part of A
 do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
 end do
!
! Step 3: compute columns of the inverse matrix C
 do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
 end do
 return
END SUBROUTINE inverse
!
real function volume(rv)
!
!  Conditions of use:
!
!  GULP is available free of charge to academic institutions
!  and non-commerical establishments only. Copies should be
!  obtained from the author only and should not be distributed
!  in any form by the user to a third party without the express
!  permission of the author. This notice applies to all parts
!  of the program, except any library routines which are
!  distributed with the code for completeness. All rights for
!  such routines remain with the original distributor.
!
!  No claim is made that this program is free from errors and
!  no liability will be accepted for any loss or damage that
!  may result. The user is responsible for checking the validity
!  of their results.
!
!  Copyright Curtin University 2004
!
  implicit none
  real, intent(in)  :: rv(3,3)
  real       :: r1x
  real       :: r1y
  real       :: r1z
  real       :: r2x
  real       :: r2y
  real       :: r2z
  real       :: r3x
  real       :: r3y
  real       :: r3z
  real       :: vol
!
  r1x = rv(1,1)
  r1y = rv(2,1)
  r1z = rv(3,1)
  r2x = rv(1,2)
  r2y = rv(2,2)
  r2z = rv(3,2)
  r3x = rv(1,3)
  r3y = rv(2,3)
  r3z = rv(3,3)
  vol = r1x*(r2y*r3z - r2z*r3y) + r1y*(r3x*r2z - r3z*r2x) + r1z*(r2x*r3y - r2y*r3x)
  volume = abs(vol)
  RETURN
end function 
end program lammpstrj2pdb
