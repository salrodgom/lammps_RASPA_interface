program pdb2cif
 implicit none
 real,    allocatable           :: xcryst(:,:),xinbox(:,:)
 integer, allocatable           :: id(:)
 character (len=4), allocatable :: label(:)
 integer                        :: n_atoms = 0
 integer                        :: n_t = 0
 integer                        :: ierr,i,j,k
 real, parameter                :: pi=acos(-1.0)
 real                           :: atom(3),ouratom(3),cell_0(1:6)
 real                           :: rv(1:3,1:3),vr(1:3,1:3),average(1:3)
 character (len=80)             :: line,string
 character (len=4)              :: mol
 character (len=6)              :: atomc
 character (len=5)              :: model
 character (len=6)              :: charct(6)
 character (len=4)              :: typ(3)
 character (len=10)             :: spacegroup = "p1"
 open(100,file="input.pdb",iostat=ierr)
 if( ierr /= 0 ) stop "[error] no se puede abrir input.pdb"
 fileopen_pdb: if( ierr == 0) then
   do
    read (100,'(a)',iostat=ierr) line
    if( ierr /= 0 ) exit
    if (line(1:5)=='MODEL') then
      n_t = n_t + 1
      n_atoms = 0
    endif
    if (line(1:4)=='ATOM') n_atoms = n_atoms + 1
   end do
 endif fileopen_pdb
 allocate(xcryst(0:3,n_atoms) ,stat=ierr)
 allocate(xinbox(3,n_atoms)   ,stat=ierr)
 allocate(id(n_atoms)         ,stat=ierr)
 allocate(label(n_atoms)      ,stat=ierr)
 if( ierr /= 0 ) stop "problema al alicatar variables"
 rewind(100)
 do
  read (100,'(a)') line
  if(line(1:4)=="CRYS") exit
 end do
 read (line,'(6x,3f9.3,3f7.2,1x,a10)') &
 cell_0(1),cell_0(2),cell_0(3),cell_0(4),cell_0(5),cell_0(6),spacegroup
 spacegroup = "p1"
 call cell(rv,vr,cell_0)
 read_coor_pdb: do i=1,n_atoms
    read(100, '(a)', iostat = ierr ) line
    if (ierr/=0) exit read_coor_pdb          ! formatos de lectura para pdb 
       atomc = line(1:6)                     ! leemos
       mol   = line(18:20)
       id(i) = i
       read(line(31:38),*) xinbox(1,id(i))
       read(line(39:46),*) xinbox(2,id(i))
       read(line(47:54),*) xinbox(3,id(i))
       label(id(i))= line(13:14)
       forall ( j = 1:3 )
        average(j)=average(j)+xinbox(j,id(i))/real(n_atoms)
       end forall
 enddo read_coor_pdb
 read (100,'(a)') line
 coor_pdb: do i=1,n_atoms ! definimos las coordenadas cartesianas
    forall ( j=1:3 )
     xinbox(j,id(i)) = xinbox(j,id(i)) - average(j)
    end forall
    forall ( j=1:3 )
     xcryst(j,id(i)) =mod(vr(j,1)*xinbox(1,id(i))+vr(j,2)*xinbox(2,id(i))+vr(j,3)*xinbox(3,id(i))+100.0,1.0)
    end forall
 enddo coor_pdb
! escribo el cif
 call escritura_cif(xcryst,n_atoms,label,cell_0,rv)
 close(100)
 stop "bob esponja!"
end program pdb2cif
subroutine escritura_cif(xcryst,n_atoms,label,cell_0,rv)
 implicit none
 integer           :: n_atoms
 real              :: xcryst(0:3,n_atoms)
 real              :: volume,cell_0(6),rv(1:3,1:3)
 integer           :: i,u
 character (len=4) :: label(1:n_atoms)
 u=1000
 open(u,file="p1.cif")
 write(u,'(a)')'data_subtitutions'
 write(u,'(a)')'_audit_creation_method    igor'
 write(u,'(a)')"_audit_author_name 'sponge bob'"
 write(u,'(a,f14.7)')'_cell_length_a',cell_0(1)
 write(u,'(a,f14.7)')'_cell_length_b',cell_0(2)
 write(u,'(a,f14.7)')'_cell_length_c',cell_0(3)
 write(u,'(a,f14.7)')'_cell_angle_alpha',cell_0(4)
 write(u,'(a,f14.7)')'_cell_angle_beta',cell_0(5)
 write(u,'(a,f14.7)')'_cell_angle_gamma',cell_0(6)
 write(u,'(a,f14.7)')'_cell_volume',volume(rv)
 write(u,'(a)')"_symmetry_space_group_name_hall 'p 1'"
 write(u,'(a)')"_symmetry_space_group_name_h-m 'p 1'"
 write(u,'(a)')'_symmetry_int_tables_number 1'
 write(u,'(a)')"_symmetry_equiv_pos_as_xyz 'x,y,z'"
 write(u,'(a)')'loop_'
 write(u,'(a)')'_atom_site_label'
 write(u,'(a)')'_atom_site_fract_x'
 write(u,'(a)')'_atom_site_fract_y'
 write(u,'(a)')'_atom_site_fract_z'
 !write(u,'(a)')'_atom_site_charge' 
 atoms_: do i=1,n_atoms
   write(u,'(a4,1x,4(f10.8,1x))')label(i),xcryst(1,i),xcryst(2,i),xcryst(3,i)
 end do atoms_
 close(u)
 return
end subroutine escritura_cif
!
subroutine cell(rv,vr,cell_0)
 implicit none
 integer :: i,j
 real, intent(in)  :: cell_0(6)
 real, intent(out) :: rv(3,3),vr(3,3)
 real :: alp,bet
 real :: cosa,cosb,cosg
 real :: gam,sing
 real :: pi,degtorad
 pi = acos(-1.0)
 degtorad=pi/180.0
 if(cell_0(4) == 90.0) then
   cosa = 0.0
 else
   alp=cell_0(4)*degtorad
   cosa=cos(alp)
 endif
 if(cell_0(5) == 90.0) then
   cosb = 0.0
 else
   bet = cell_0(5)*degtorad
   cosb = cos(bet)
 endif
 if(cell_0(6) == 90.0) then
   sing = 1.0
   cosg = 0.0
 else
   gam = cell_0(6)*degtorad
   sing = sin(gam)
   cosg = cos(gam)
 endif
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
 return
end subroutine cell
!
subroutine uncell(rv,cell_0)
  implicit none
  real,intent(out) :: cell_0(6)
  real,intent(in)  :: rv(3,3)
  integer  :: i,j
  real     :: temp(6)
  real :: radtodeg,pi
  pi=acos(-1.0)
  radtodeg=180.0/pi
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
  do i=4,6
     if (abs(cell_0(i) - 90.0 ).lt.0.00001) cell_0(i) = 90.0
     if (abs(cell_0(i) - 120.0).lt.0.00001) cell_0(i) = 120.0
  end do
  return
end subroutine uncell
!
subroutine inverse(a,c,n)
 implicit none 
 integer n
 real a(n,n), c(n,n)
 real l(n,n), u(n,n), b(n), d(n), x(n)
 real coeff
 integer i, j, k
 l=0.0
 u=0.0
 b=0.0
 do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      l(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
 end do
 do i=1,n
  l(i,i) = 1.0
 end do
 do j=1,n
  do i=1,j
    u(i,j) = a(i,j)
  end do
 end do
 do k=1,n
  b(k)=1.0
  d(1) = b(1)
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - l(i,j)*d(j)
    end do
  end do
  x(n)=d(n)/u(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-u(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
 end do
 return
end subroutine inverse
!
real function volume(rv)
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
  return
end function
