program zif_cif2gin
 use iso_fortran_env
 implicit none
! locals
 integer             :: i,j,k,l,h,m,ierr
 integer             :: ii,jj,kk,ll
 real                :: r = 1.0e12
!parameters
 real,parameter      :: k_B = 8.617332478e-5
 real,parameter      :: r_min_criteria_connectivity=0.2
 INTEGER, PARAMETER  :: muchisimo=100000
! variables
 integer             :: num_args
 integer             :: n_atoms = 0
 real                :: r0=r_min_criteria_connectivity
 real                :: ratom(3),rouratom(3)
 real :: rv(3,3),vr(3,3),cell_0(1:6) = 0.0
 real :: xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound
 real :: xy,xz,yz
 character(len=100)  :: line
 character(len=20)   :: spam,simulation_type="single"
 character(len=80)   :: string_stop_head= "_atom_site_charge" !"_atom_site_charge"
 character(len=100)  :: CIFFilename=" "
 character(len=100)  :: filename=" "
 CHARACTER (LEN=18)  :: string(muchisimo)
 integer,parameter   :: n_atom_types_max=100
 integer             :: n_atom_types = 0
 character(len=4)    :: atom_types(n_atom_types_max)
 ! type
 integer,parameter   :: max_n_componets=2
 type  :: node
  integer          :: element
  integer          :: type_
  character(len=2) :: label_element
  character(len=4) :: label
  character(len=4) :: new_label="Xxxx"
  integer          :: degree
  real             :: charge 
  real             :: radius
  real             :: mass
  integer          :: n_components
  character(5)     :: clabel(max_n_componets)
  real             :: xyzs(1:3,max_n_componets)
  real             :: xyzc(1:3,max_n_componets)
 end type 
 ! allocatable
 type(node),allocatable,dimension(:) :: atom
 type(node),allocatable,dimension(:) :: guest_atom
 real,allocatable    :: DistanceMatrix(:,:)
 logical,allocatable :: ConnectedAtoms(:,:)
 integer :: bond_types_max=100, bend_types_max=100, tors_types_max=100, impr_types_max=100
! bonds
 character(len=8),dimension(:), allocatable   :: bond_type_string
 integer,allocatable                          :: bond_type_histogram(:)
 character(len=8)                             :: bond_string
 integer                                      :: n_bond_types = 0
 character(len=80),dimension(:), allocatable  :: bond
 !character(len=8),dimension(:), allocatable   :: bond_type
 integer                                      :: n_bonds = 0
 integer                                      :: n_bonds_max = 100000
! bends
 character(len=12),dimension(:), allocatable  :: bend_type_string
 integer,allocatable                          :: bend_type_histogram(:)
 character(len=12)                            :: bend_string
 integer                                      :: n_bend_types = 0
 character(len=80),dimension(:), allocatable  :: bend
 !character(len=20),dimension(:), allocatable  :: bend_type
 integer                                      :: n_bends = 0
 integer                                      :: n_bends_max = 100000
! torsion (dihedrals)
 character(len=16),dimension(:), allocatable  :: tors_type_string
 integer,allocatable                          :: tors_type_histogram(:)
 character(len=16)                            :: tors_string
 integer                                      :: n_tors_types = 0
 character(len=80),dimension(:), allocatable  :: tors
 !character(len=25),dimension(:), allocatable  :: tors_type
 integer                                      :: n_torss = 0
 integer                                      :: n_torss_max = 100000
! torsion (impropers)
 character(len=16),dimension(:), allocatable  :: impr_type_string
 integer,allocatable                          :: impr_type_histogram(:)
 character(len=16)                            :: impr_string
 integer                                      :: n_impr_types = 0
 character(len=80),dimension(:), allocatable  :: impr
 !character(len=25),dimension(:), allocatable  :: impr_type
 integer                                      :: n_imprs = 0
 integer                                      :: n_imprs_max = 100000
! arguments in line
 character(len=100),dimension(:), allocatable :: args
 logical                                      :: charges_flag         = .true.
 logical                                      :: modify_topology_flag = .false.
 logical                                      :: fill_with_Ar_flag    = .false.
 !
 num_args = command_argument_count()
 allocate(args(num_args))
 do i = 1, num_args
  call get_command_argument(i,args(i))
 end do
 write(6,'(a,1x,i2)')'Arguments:',command_argument_count()
 write(6,'(a)')( args(i),i=1,num_args)
 if(num_args==0) then
  call print_help()
  stop
 end if
 do i=1,num_args
  select case(args(i))
   case ('-h','--help')
    call print_help()
    stop
   case ('-c','--cif')
    CIFFilename=args(i+1)
    filename=CIFFilename(1:Clen_trim(CIFFilename)-4)
    write(6,'(a)') filename
   case ('-wq','--no-charges')
    charges_flag=.false.
    string_stop_head= "_atom_site_fract_z"
   case ('-S','--search-and-modify-topology')
    modify_topology_flag = .true.
   case ('-F','--fill-with-Ar')
    fill_with_Ar_flag = .true.
  end select
 end do
 open(100,file=CIFfilename,status='old',iostat=ierr)
 if(ierr/=0) stop 'Error opening CIF file'
 read_cif: do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0) exit read_cif
  if(line(1:14)=="_cell_length_a")then
   read(line,*)spam,cell_0(1)
   cycle read_cif
  end if
  if(line(1:14)=="_cell_length_b")then
   read(line,*)spam,cell_0(2)
   cycle read_cif
  end if
  if(line(1:14)=="_cell_length_c")then
   read(line,*)spam,cell_0(3)
   cycle read_cif
  end if
  if(line(1:17)=="_cell_angle_alpha")then
   read(line,*)spam,cell_0(4)
   cycle read_cif
  end if
  if(line(1:16)=="_cell_angle_beta")then
   read(line,*)spam,cell_0(5)
   cycle read_cif
  end if
  if(line(1:17)=="_cell_angle_gamma")then
   read(line,*)spam,cell_0(6)
   cycle read_cif
  end if
  if(line(1:)==string_stop_head) exit read_cif
 end do read_cif
 call cell(rv,vr,cell_0)
 n_atoms=0
 read_natoms: do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0) exit read_natoms
  n_atoms=n_atoms+1
 end do read_natoms
!
 allocate(atom(n_atoms))
 atom(1:n_atoms)%new_label="Xxxx"
 allocate(ConnectedAtoms(n_atoms,n_atoms))
 allocate(DistanceMatrix(n_atoms,n_atoms))
!
 rewind(100)
 write(6,*)'Atoms:',n_atoms
 do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0) exit
  if(line(1:)==string_stop_head) exit
 end do
 forall (i=1:n_atom_types_max) 
  forall (j=1:4)
   atom_types(i)(j:j)=" "
  end forall
 end forall
 i=0
 read_atoms: do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0) exit
  i=i+1
  atom(i)%n_components=1
  if(charges_flag) then
   read(line,*)atom(i)%label,(atom(i)%xyzs(j,1),j=1,3),atom(i)%charge
  else
   read(line,*)atom(i)%label,(atom(i)%xyzs(j,1),j=1,3)
   atom(i)%charge=0.0
  end if
  call CheckAtom(atom(i)%label,atom(i)%mass,atom(i)%radius,atom(i)%element,atom(i)%label_element)
  do j=1,3
   atom(i)%xyzs(j,1)=mod(atom(i)%xyzs(j,1)+100.0,1.0)
   atom(i)%xyzc(j,1)=rv(j,1)*atom(i)%xyzs(1,1)+rv(j,2)*atom(i)%xyzs(2,1)+rv(j,3)*atom(i)%xyzs(3,1)
  end do
 end do read_atoms
 close(100)
 if ( fill_with_Ar_flag ) then
  call FillStructureWithAr()
 end if
! topology:
 DistanceMatrix=0.0
 ConnectedAtoms=.false.
 do i=1,n_atoms
  k=0
  do j=1,n_atoms 
   forall (h=1:3)
    rouratom(h)=atom(i)%xyzs(h,1)
    ratom(h)=atom(j)%xyzs(h,1)
   end forall
   call make_distances(cell_0,ratom,rouratom,rv,r)
   DistanceMatrix(i,j)=r
   if(r>0.1.and.r<=atom(i)%radius+atom(j)%radius+r_min_criteria_connectivity)then
    !.and.&
    !  (atom(i)%element/=1.and.atom(j)%element/=1) )then
    k=k+1
    ConnectedAtoms(i,j)=.true.
   end if
  end do
  atom(i)%degree=k
 end do
! degree of each node:
 write(6,'(a)') 'Connectivity array for each atom:'
 write(6,'(1000(i2))')(atom(i)%degree,i=1,n_atoms)
! rename atoms with for certain forcefield
 if ( modify_topology_flag ) then
! charge for topologies:
!charge +0.6918  # Zn
!charge +0.2045  # C4
!charge -0.3879  # N1
!charge -0.0839  # C2
!charge +0.1310  # H2
!charge +0.0     # C5
!charge -0.10745 # C6
!charge +0.1310  # H1
!
  ! {{
  scan_for_rename_first_phase: do i=1,n_atoms
   !first N and Zn:
   if( atom(i)%element==7 )then
    atom(i)%new_label='N1  '
    atom(i)%charge = -0.3879
   else if( atom(i)%element==30) then
    atom(i)%new_label='Zn  '
    atom(i)%charge = +0.6918
   end if
  end do scan_for_rename_first_phase

  scan_for_rename_carbon_atoms: do i=1,n_atoms
   !second, C-atoms:
   if( atom(i)%element==6 ) then
    h=0 ! N-counter
    l=0 ! C-counter
    scan_for_N_C_atoms: do j=1,n_atoms
     if( i/=j.and.ConnectedAtoms(i,j) )then
      if( atom(j)%element==7 ) h=h+1 
      if( atom(j)%element==6 ) l=l+1
     end if
    end do scan_for_N_C_atoms
    if( h==1 .and. l==1 ) then
     atom(i)%new_label = "C2  "
     atom(i)%charge = -0.0839
    else if ( h==2 .and. l==0 ) then
     atom(i)%new_label = "C4  "
     atom(i)%charge =  +0.259300001
    else if ( h==2 .and. l==1 ) then
     atom(i)%new_label = "C1  "
     atom(i)%charge = +0.4291
    else if ( h==0 .and. l==1 ) then
     atom(i)%new_label = "C3  "
     atom(i)%charge = -0.4526
    else if ( h==1 .and. l==2 ) then
     atom(i)%new_label = "C5  "
     atom(i)%charge = +0.0
    else if ( h==0 .and. l==2 ) then
     atom(i)%new_label = "C6  "
     atom(i)%charge = -0.09835
    else
     stop 'unknow C-atom'
    end if
   end if
  end do scan_for_rename_carbon_atoms
  
  scan_for_rename_H_atoms: do i = 1,n_atoms
   if ( atom(i)%element == 1 ) then
    do j=1,n_atoms
     if(j/=i.and.ConnectedAtoms(i,j))then
      if( atom(j)%new_label == "C2  " )      then
       atom(i)%new_label = "H2  "
       atom(i)%charge = +0.1128
      else if( atom(j)%new_label == "C4  " ) then
       atom(i)%new_label = "H2  "
       atom(i)%charge = +0.1128
      else if( atom(j)%new_label == "C3  " ) then
       atom(i)%new_label = "H3  "
       atom(i)%charge = +0.131866667
      else if( atom(j)%new_label == "C6  " ) then
       atom(i)%new_label = "H1  "
       atom(i)%charge = +0.1128
      else
       stop 'unknow H-atom'
      end if
      cycle scan_for_rename_H_atoms
     end if
    end do
   end if
  end do scan_for_rename_H_atoms 
  do i=1,n_atoms
   if(atom(i)%new_label/="Xxxx")then
    atom(i)%label=atom(i)%new_label
   end if
  end do 
  ! }}
 end if
! atom types:
 do i = 1,n_atoms
  check_atom_type: if(i==1)then
   n_atom_types=1
   atom_types(n_atom_types)=atom(i)%label
   atom(i)%type_=n_atom_types
  else
   n_atom_types=n_atom_types+1
   do h=1,n_atom_types-1
    if(atom(i)%label==atom_types(h))then
     n_atom_types=n_atom_types-1
     atom(i)%type_=h
     exit check_atom_type
    end if
   end do
   atom_types(n_atom_types)=atom(i)%label
   atom(i)%type_=n_atom_types
  end if check_atom_type
  write(6,'(a4,1x,3(f10.8,1x),a2,1x,3(f10.6,1x))') &
   atom(i)%label,(atom(i)%xyzs(j,1),j=1,3),&
   atom(i)%label_element,atom(i)%charge,atom(i)%radius,atom(i)%mass
 end do
 write(6,*)'=========='
 write(6,*)'atom_types:',n_atom_types
 do i=1,n_atom_types
  write(6,*)atom_types(i)
 end do
 open(987,file="atom_types_for_dump.txt")
 write(987,'(100(a4,1x))') ( atom_types(i), i=1,n_atom_types)
 close(987)
 write(6,*)'=========='
! bonds: 
 allocate(bond_type_string(bond_types_max))
 allocate(bond_type_histogram(bond_types_max))
 allocate(bond(n_bonds_max))
 bond_type_histogram=0
 bond_string(1:8)="        "
 do i=1,bond_types_max
  do j=1,8 
   write(bond_type_string(j:j),'(a1)')" "
  end do 
 end do
 forall (i=1:n_bonds_max)
  forall (j=1:80)
   bond(i)(j:j)=" "
  end forall
 end forall
 do_i: do i=1,n_atoms
  do_j: do j=i+1,n_atoms
   if(ConnectedAtoms(i,j))then
    write(bond_string(1:8),'(2a4)') atom(i)%label(1:4),atom(j)%label(1:4)
    forall (l=1:100)
     line(l:l)=" "
    end forall
    if(n_bond_types==0)then
     n_bonds=n_bonds+1
     n_bond_types=n_bond_types+1
     bond_type_string(1)=bond_string
     bond_type_histogram(1)=1
     write(bond(n_bonds)(1:15),'(3i5)')n_bond_types,i,j
     write(bond(n_bonds)(16:18),'(a3)')' # '
     write(bond(n_bonds)(19:),'(a)') bond_string
     cycle do_j
    else
     n_bonds=n_bonds+1
     n_bond_types=n_bond_types+1 ! try
     check_bond: do h=n_bond_types-1,1,-1
      if(bond_string(1:8)==bond_type_string(h)(1:8).or.&
         bond_string(1:8)==bond_type_string(h)(5:8)//bond_type_string(h)(1:4))then
       n_bond_types=n_bond_types-1
       bond_type_histogram(h)=bond_type_histogram(h)+1
       write(bond(n_bonds)(1:15),'(3i5)')h,i,j
       write(bond(n_bonds)(16:18),'(a3)')' # '
       write(bond(n_bonds)(19:),'(a)') bond_string
       cycle do_j
      end if
     end do check_bond
     bond_type_histogram(n_bond_types)=1
     bond_type_string(n_bond_types)=bond_string
     write(bond(n_bonds)(1:15),'(3i5)')n_bond_types,i,j
     write(bond(n_bonds)(16:18),'(a3)')' # '
     write(bond(n_bonds)(19:),'(a)') bond_string
    end if
   end if
  end do do_j
 end do do_i
! bends:
 allocate(bend_type_string(bend_types_max))
 allocate(bend_type_histogram(bend_types_max))
 allocate(bend(n_bends_max))
 bend_type_histogram=0
 bend_string(1:12)="            "
 do i=1,bend_types_max
  do j=1,12
   write(bend_type_string(j:j),'(a1)')" "
  end do
 end do
 forall (i=1:n_bends_max)
  forall (j=1:80)
   bend(i)(j:j)=" "
  end forall
 end forall
 do_i_bend: do i=1,n_atoms ! central atom
  do_j_bend: do j=1,n_atoms
   if(ConnectedAtoms(i,j).and.j/=i)then
    do_k_bend: do k=1,n_atoms
     if(ConnectedAtoms(i,k).and.k/=j.and.k/=i)then
      write(bend_string(1:12),'(3a4)') atom(j)%label(1:4),atom(i)%label(1:4),atom(k)%label(1:4)
      if(n_bend_types==0)then
       n_bends=n_bends+1
       n_bend_types=n_bend_types+1
       bend_type_string(1)=bend_string
       bend_type_histogram(1)=1
       write(bend(n_bends)(1:20),'(4i5)')n_bend_types,j,i,k
       write(bend(n_bends)(21:23),'(a3)')' # '
       write(bend(n_bends)(24:),'(a)') bend_string
       cycle do_k_bend
      else
       n_bends=n_bends+1
       n_bend_types=n_bend_types+1 ! try
       check_bend_type: do h=n_bend_types-1,1,-1
        if(bend_string(1:12)==bend_type_string(h)(1:12).or.&
           bend_string(1:12)==bend_type_string(h)(9:12)//bend_type_string(h)(5:8)//&
           bend_type_string(h)(1:4))then
         n_bend_types=n_bend_types-1
         write(bend(n_bends)(1:20),'(4i5)')h,j,i,k
         write(bend(n_bends)(21:23),'(a3)')' # '
         write(bend(n_bends)(24:),'(a)') bend_string
         check_bend: do l=1,n_bends-1
          read(bend(l)(1:15),'(3i5)')jj,ii,kk
          if((jj==j.and.ii==i.and.kk==k).or.&
             (kk==j.and.ii==i.and.jj==k))then
           n_bends=n_bends-1
           cycle do_k_bend
          end if
         end do check_bend
         bend_type_histogram(h)=bend_type_histogram(h)+1
         cycle do_k_bend
        end if
       end do check_bend_type
       bend_type_histogram(n_bend_types)=1
       bend_type_string(n_bend_types)=bend_string
       write(bend(n_bends)(1:20),'(4i5)')n_bend_types,j,i,k
       write(bend(n_bends)(21:23),'(a3)')' # '
       write(bend(n_bends)(24:),'(a)') bend_string
      end if
     end if
    end do do_k_bend
   end if
  end do do_j_bend
 end do do_i_bend
! dihedrals
 allocate(tors_type_string(tors_types_max))
 allocate(tors_type_histogram(tors_types_max))
 allocate(tors(n_torss_max))
 tors_type_histogram=0
 tors_string(1:16)="                "
 do i=1,tors_types_max
  do j=1,16
   write(tors_type_string(j:j),'(a1)')" "
  end do
 end do
 forall (i=1:n_torss_max)
  forall (j=1:80)
   tors(i)(j:j)=" "
  end forall
 end forall
 do_i_tors: do i=1,n_atoms
  do_j_tors: do j=1,n_atoms
   if(ConnectedAtoms(i,j).and.j/=i)then
    do_k_tors: do k=1,n_atoms
     if(ConnectedAtoms(j,k).and.k/=j.and.k/=i)then
      do_l_tors: do l=1,n_atoms
       if(ConnectedAtoms(k,l).and.l/=k.and.l/=j.and.l/=i)then
        write(tors_string(1:16),'(4a4)')atom(i)%label(1:4),atom(j)%label(1:4),&
                                        atom(k)%label(1:4),atom(l)%label(1:4)
        if(n_tors_types==0)then
         n_torss=n_torss+1
         n_tors_types=n_tors_types+1
         tors_type_string(1)=tors_string
         tors_type_histogram(1)=1
         write(tors(n_torss)(1:25),'(5i5)')n_tors_types,i,j,k,l
         write(tors(n_torss)(26:28),'(a3)')' # '
         write(tors(n_torss)(29:),'(a)') tors_string
         cycle do_l_tors
        else
         n_torss=n_torss+1
         n_tors_types=n_tors_types+1 ! try
         check_tors_type: do h=n_tors_types-1,1,-1
          if(tors_string(1:16)==tors_type_string(h)(1:16).or.&
             tors_string(1:16)==tors_type_string(h)(13:16)//tors_type_string(h)(9:12)//&
                                tors_type_string(h)(5:8)//tors_type_string(h)(1:4))then
           n_tors_types=n_tors_types-1
           write(tors(n_torss)(1:25),'(5i5)')h,i,j,k,l
           write(tors(n_torss)(26:28),'(a3)')' # '
           write(tors(n_torss)(29:),'(a)') tors_string
           check_tors: do m=1,n_torss-1
            read(tors(m)(1:20),'(4i5)')ii,jj,kk,ll
            if((ii==i.and.jj==j.and.kk==k.and.ll==l).or.&
               (ll==i.and.kk==j.and.jj==k.and.ii==l))then
             n_torss=n_torss-1
             cycle do_l_tors
            end if
           end do check_tors
           tors_type_histogram(h)=tors_type_histogram(h)+1
           cycle do_l_tors
          end if
         end do check_tors_type
         tors_type_histogram(n_tors_types)=1
         tors_type_string(n_tors_types)=tors_string
         write(tors(n_torss)(1:25),'(5i5)') n_tors_types,i,j,k,l
         write(tors(n_torss)(26:28),'(a3)')' # '
         write(tors(n_torss)(29:),'(a)') tors_string
        end if
       end if
      end do do_l_tors
     end if
    end do do_k_tors
   end if
  end do do_j_tors
 end do do_i_tors
! impropers
 allocate(impr_type_string(impr_types_max))
 allocate(impr_type_histogram(impr_types_max))
 allocate(impr(n_imprs_max))
 impr_type_histogram=0
 impr_string(1:16)="                "
 do i=1,impr_types_max
  do j=1,16
   write(impr_type_string(j:j),'(a1)')" "
  end do
 end do
 forall (i=1:n_imprs_max)
  forall (j=1:80)
   impr(i)(j:j)=" "
  end forall
 end forall
 do_i_impr: do i=1,n_atoms
  do_j_impr: do j=1,n_atoms
   if(ConnectedAtoms(i,j).and.j/=i)then
    do_k_impr: do k=1,n_atoms
     if(ConnectedAtoms(i,k).and.k/=j.and.k/=i)then
      do_l_impr: do l=1,n_atoms
       if(ConnectedAtoms(i,l).and.l/=k.and.l/=j.and.l/=i)then
        write(impr_string(1:16),'(4a4)')atom(i)%label(1:4),atom(j)%label(1:4),&
                                        atom(k)%label(1:4),atom(l)%label(1:4)
        if(n_impr_types==0)then
         n_imprs=n_imprs+1
         n_impr_types=n_impr_types+1
         impr_type_string(1)=impr_string
         impr_type_histogram(1)=1
         write(impr(n_imprs)(1:25),'(5i5)')n_impr_types,i,j,k,l
         write(impr(n_imprs)(26:28),'(a3)')' # '
         write(impr(n_imprs)(29:),'(a)') impr_string
         cycle do_l_impr
        else
         n_imprs=n_imprs+1
         n_impr_types=n_impr_types+1 ! try
         check_impr_type: do h=n_impr_types-1,1,-1
          if(impr_string(1:16)==impr_type_string(h)(1:16).or.&
             impr_string(1:16)==impr_type_string(h)(1:4)//impr_type_string(h)(9:12)//&
                                impr_type_string(h)(5:8)//impr_type_string(h)(13:16).or.&
             impr_string(1:16)==impr_type_string(h)(1:4)//impr_type_string(h)(5:8)//&
                                impr_type_string(h)(13:16)//impr_type_string(h)(9:12).or.&
             impr_string(1:16)==impr_type_string(h)(1:4)//impr_type_string(h)(13:16)//&
                                impr_type_string(h)(9:12)//impr_type_string(h)(5:8) )then
           n_impr_types=n_impr_types-1
           write(impr(n_imprs)(1:25),'(5i5)')h,i,j,k,l
           write(impr(n_imprs)(26:28),'(a3)')' # '
           write(impr(n_imprs)(29:),'(a)') impr_string
           check_impr: do m=1,n_imprs-1
            read(impr(m)(1:20),'(4i5)')ii,jj,kk,ll
            if((ii==i.and.jj==j.and.kk==k.and.ll==l).or.&
               (ii==i.and.kk==j.and.jj==k.and.ll==l).or.&
               (ii==i.and.jj==j.and.ll==k.and.kk==l).or.&
               (ii==i.and.ll==j.and.jj==k.and.kk==l).or.&
               (ii==i.and.kk==j.and.ll==k.and.jj==l).or.&
               (ii==i.and.ll==j.and.kk==k.and.jj==l))then
             n_imprs=n_imprs-1
             cycle do_l_impr
            end if
           end do check_impr
           impr_type_histogram(h)=impr_type_histogram(h)+1
           cycle do_l_impr
          end if
         end do check_impr_type
         impr_type_histogram(n_impr_types)=1
         impr_type_string(n_impr_types)=impr_string
         write(impr(n_imprs)(1:25),'(5i5)') n_impr_types,i,j,k,l
         write(impr(n_imprs)(26:28),'(a3)')' # '
         write(impr(n_imprs)(29:),'(a)') impr_string
        end if
       end if
      end do do_l_impr
     end if
    end do do_k_impr
   end if
  end do do_j_impr
 end do do_i_impr
 write(6,*)'Bond types:',n_bond_types
 write(6,*)'bonds:',n_bonds
 do i=1,n_bond_types 
  write(6,*) bond_type_string(i),bond_type_histogram(i)
 end do
 write(6,*)'Bend types:',n_bend_types
 write(6,*)'bends:',n_bends
 do i=1,n_bend_types 
  write(6,*) bend_type_string(i),bend_type_histogram(i)
 end do
 write(6,*)'Dihedral types:',n_tors_types
 write(6,*)'dihedrals:',n_torss
 do i=1,n_tors_types
  write(6,*) tors_type_string(i),tors_type_histogram(i)
 end do
 write(6,*)'Improper types:',n_impr_types
 write(6,*)'impropers:',n_imprs
 do i=1,n_impr_types
  write(6,*) impr_type_string(i),impr_type_histogram(i)
 end do
! output:
 call cellnormal2lammps(cell_0,xlo_bound,ylo_bound,zlo_bound,&
      xhi_bound,yhi_bound,zhi_bound,xy,xz,yz)
 call output_lammps()
 call output_gulp()
 call output_pdb()
 deallocate(bend_type_string)
 deallocate(bend_type_histogram)
 deallocate(bend)
 deallocate(bond_type_string)
 deallocate(bond_type_histogram)
 deallocate(DistanceMatrix)
 deallocate(ConnectedAtoms)
 deallocate(atom)
 stop 'Done'
 contains
 !
 subroutine FillStructureWithAr()
  implicit none
  write(6,*)'Filling structure with Ar!'
  call FillCellWithAr()
 end subroutine FillStructureWithAr
 !
 subroutine FillCellWithAr()
  implicit none
  integer :: i,kx,ky,kz
  logical :: Ar_allocated =.false.
  integer :: n_filling_atoms = 0
  real    :: x,y,z
  real    :: Ar_coor(3,1000) = 0.0
  real    :: cell_1(6)
  cell_1(1)=5.235
  cell_1(2)=5.235
  cell_1(3)=5.235
  cell_1(4)=90.0 
  cell_1(5)=90.0
  cell_1(6)=90.0
  x=0.0
  y=0.0
  z=0.0
  kx=1
  ky=1
  kz=1
  filling_cycle: do
   if(kx==1.and.ky==1.and.kz==1)then
    n_filling_atoms = 4
!
    Ar_coor(1,1)=0.0
    Ar_coor(2,1)=0.0
    Ar_coor(3,1)=0.0
    !
    Ar_coor(1,2)=0.0
    Ar_coor(2,2)=0.5
    Ar_coor(3,2)=0.5
    !
    Ar_coor(1,3)=0.5
    Ar_coor(2,3)=0.0
    Ar_coor(3,3)=0.5
!
    Ar_coor(1,4)=0.5
    Ar_coor(2,4)=0.5
    Ar_coor(3,4)=0.0
   end if
  end do filling_cycle
  allocate(guest_atom(n_filling_atoms))
  return
 end subroutine FillCellWithAr
!
 subroutine search_forcefield(string,interaction,label)
  implicit none
  character(len=80),intent(out):: string
  real                         :: p1,p2,p3
  integer                      :: p0,p4
  character(len=10)            :: type_fff
  character(len=4),intent(in)  :: interaction
  character(len=4),intent(in)  :: label(1:4)
  character(len=4)             :: ourlabel(1:4)
  character(len=80)            :: units
  character(len=80)            :: passing
  character(len=80)            :: variables
  integer                      :: u=123,ierr=0
  forall (ii=1:80)
   string(ii:ii)=" "
   units(ii:ii)=" "
   passing(ii:ii)=" "
   variables(ii:ii)=" "
  end forall
  open(u,file='forcefield.lib')
  fff: select case (interaction)
   case('bond')
    initbond: do
     read(u,'(a)') line
     if(line(1:11)=='Bond Coeffs') then
      read(line(12:),'(a)') units
      units=adjustl(trim(units))
      exit initbond
     end if
    end do initbond
    do
     read(u,'(a)') line
     if(line(1:11)=='Angle Coeffs'.or.line(1:3)=="End") then
      string(1:4)="none"
      exit fff
     else if(line(1:1)/="#".or.line(1:1)/=" ")then
      read(line,'(2(a4,1x),a)')ourlabel(1),ourlabel(2),passing
      if((adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(2)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(1)))))then
       read(passing,'(a,a)')type_fff,variables
       if(adjustl(trim(type_fff))/='none')then 
        read(variables,*) p1,p2
        select case(units)
         case('kcal/mol','kcal/mol/A/A')
          p1=kcalmol2ev(p1) ! kcal/mol/A/A
          p2=p2             ! A
        end select
        write(string,*) p1,p2,' # ',type_fff
       else
        write(string,*) 0.0,0.0,' # ',type_fff
       end if
       exit fff
      end if
     end if
    end do
   case('bend')
    initbend: do
     read(u,'(a)') line
     if(line(1:12)=='Angle Coeffs') then
      read(line(13:),'(a)') units
      units=adjustl(trim(units))
      !write(6,*) units
      exit initbend
     end if
    end do initbend
    do
     read(u,'(a)') line
     if(line(1:15)=='Dihedral Coeffs'.or.line(1:3)=="End") then
      string(1:4)="none"
      exit fff
     end if
     if(line(1:1)/="#".or.line(1:1)/=" ")then
      read(line,'(3(a4,1x),a)')ourlabel(1),ourlabel(2),ourlabel(3),passing
      if((adjustl(trim(ourlabel(2)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(3)))).or.&
         (adjustl(trim(ourlabel(2)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(1)))==adjustl(trim(label(3))).and.&
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(1)))))then
       read(passing,'(a,a)')type_fff,variables
       if(adjustl(trim(type_fff))/='none')then
        read(variables,*) p1,p2
        select case(units)
         case('kcal/mol','kcal/mol/rad/rad')
          p1=kcalmol2ev(p1) ! kcal/mol/rad/rad
          p2=p2             ! deg
        end select
        write(string,*) p1,p2, ' # ',type_fff
       else
        write(string,*)0.0,0.0,' # ',type_fff
       end if
       exit fff
      end if
     end if
    end do
   case('tors')
    inittors: do
     read(u,'(a)') line
     if(line(1:15)=='Dihedral Coeffs') then
      read(line(16:),'(a)') units
      units=adjustl(trim(units))
      exit inittors
     end if
    end do inittors
    do
     read(u,'(a)') line
     if(line(1:15)=='Improper Coeffs'.or.line(1:3)=="End") then
      string(1:4)="none"
      exit fff
     else if(line(1:1)/="#".or.line(1:1)/=" ")then
      read(line,'(4(a4,1x),a)')(ourlabel(ii),ii=1,4),passing
      if((adjustl(trim(ourlabel(2)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(3))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(4)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(4))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(3))).and.&
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(1)))))then
       read(passing,'(a,a)')type_fff,variables
       if(adjustl(trim(type_fff))/='none')then
        read(variables,*) p0,p1,p4,p2
        select case(units)
         case('kcal/mol')
          p0=p0
          p1=kcalmol2ev(p1) ! kcal/mol/rad/rad
          p4=p4
          p2=p2             ! deg
        end select
        write(string,*)p0,p1,p4,p2,' # ',type_fff
       else
        write(string,*)1,0.0,1,0.0,' # ',type_fff
       end if
       exit fff
      end if
     end if
    end do
   case('impr')
    initimpr: do
     read(u,'(a)') line
     if(line(1:15)=='Improper Coeffs') then
      read(line(16:),'(a)') units
      units=adjustl(trim(units))
      exit initimpr
     end if
    end do initimpr
    do
     read(u,'(a)',iostat=ierr) line
     if (ierr/=0) exit fff
     if(line(1:11)=='Pair Coeffs'.or.line(1:3)=="End") then
      string(1:4)="none"
      exit fff
     else if(line(1:1)/="#".or.line(1:1)/=" ")then
      read(line,'(4(a4,1x),a)')(ourlabel(ii),ii=1,4),passing
      if((adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(3))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(4)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(4))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(3)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(4))).and.&
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(3))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(2)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(3))).and.&
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(4)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(4))).and.&
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(3)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(3))).and.&
          adjustl(trim(ourlabel(3)))==adjustl(trim(label(4))).and.&
          adjustl(trim(ourlabel(4)))==adjustl(trim(label(2)))) )then
       read(passing,'(a,a)')type_fff,variables
       if(adjustl(trim(type_fff))/='none')then
        read(variables,*) p1,p2,p3
        select case(units)
         case('kcal/mol')
          p1=kcalmol2ev(p1) ! kcal/mol/rad/rad
          p2=p2             ! sign (+/-)
          p3=p3             ! -
        end select
        write(string,*) p1,p2,p3,' # ',type_fff
       else
        write(string,*)0.0,-1,0,' # ',type_fff
       end if
       exit fff
      end if
     end if
    end do
   case('pair')
    initpair: do
     read(u,'(a)') line
     if(line(1:11)=='Pair Coeffs') then
      read(line(12:),'(a)') units
      units=adjustl(trim(units))
      exit initpair
     end if
    end do initpair
    do
     read(u,'(a)',iostat=ierr) line
     if(ierr/=0) exit fff
     if(line(1:3)=="End") then
      string(1:4)="none"
      exit fff
     else if(line(1:1)/="#".or.line(1:1)/=" ")then
      read(line,'(2(a4,1x),a)')ourlabel(1),ourlabel(2),passing
      if((adjustl(trim(ourlabel(1)))==adjustl(trim(label(1))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(2)))).or.&
         (adjustl(trim(ourlabel(1)))==adjustl(trim(label(2))).and.&
          adjustl(trim(ourlabel(2)))==adjustl(trim(label(1)))))then
       read(passing,'(a,a)')type_fff,variables
       if(adjustl(trim(type_fff))/='none')then
        read(variables,*) p1,p2
        !write(6,*)p1,p2
        select case(units)
         case('kcal/mol')
          p1=kcalmol2ev(p1) ! kcal/mol
          p2=p2             ! A
        end select
        write(string,*) p1,p2,' # ',type_fff
       else
        write(string,*) 0.0,0.0,' # ',type_fff
       end if
       exit fff
      end if
     end if
    end do
  end select fff
  close(u)
  return
 end subroutine search_forcefield
 subroutine output_pdb()
  implicit none
  integer           :: u=333,i
  character(len=4)  :: extension=".pdb"
  character(len=80) :: PDBFilename
  PDBFilename=filename(1:Clen_trim(filename))//extension
  PDBFilename=adjustl(PDBfilename)
  open(u,file=PDBFilename)
  write(u,'(a5,i5)')'MODEL',k
  write(u,'(a6,9(f14.7,1x))')'REMARK',&
   xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound,xy,xz,yz
  write(u,'(a6,1x,a,1x,i6,1x,a,i6)')'REMARK','# atoms',n_atoms,'# timestep',int(0.0)
  write(u,'(a6,1x,a,1x,f14.7)')'REMARK','Volume:',volume(rv)
  write(u,'(a6,3f9.3,3f7.2)')'CRYST1',(cell_0(j),j=1,6)
  do i=1,n_atoms
   write(u,'(a6,i5,1x,a2,3x,a4,1x,i4,4x,3f8.3,2f6.2,10x,a2)') &
   'ATOM  ',i,atom(i)%label,'MOL ',0,(atom(i)%xyzc(j,1),j=1,3),0.0,0.0,atom(i)%label
  end do
  write(u,'(a6)')'ENDMDL'
  close(u)
 end subroutine output_pdb
 subroutine output_gulp()
  implicit none
  character(len=80) :: GULPFilename
  integer           :: u=444
  integer           :: i
  real              :: mmm,rrr
  integer           :: zzz
  character(len=2)  :: zlz
  character(len=4)  :: extension=".gin"
  GULPFilename=filename(1:Clen_trim(filename))//extension
  GULPFilename=adjustl(GULPfilename)
  !adjustl(
  open(u,file=GULPFilename)
  write(u,'(a)')'single conv molq qok noenergy'
  write(u,'(A)')'cell'
  write(u,'(6(f9.5,1x))') (cell_0(j) , j=1,6)
  write(u,'(A)')'fractional'
  do i=1,n_atoms
   write(u,'(a4,1x,3(f14.7,1x),1x,f14.7)')atom(i)%label,(atom(i)%xyzs(j,1),j=1,3),atom(i)%charge
  end do
  write(u,'(A,1x,i5)')'species',n_atom_types
  do i=1,n_atom_types
   write(u,'(a)')atom_types(i)
  end do
  write(u,'(a,1x,a)')'output lammps ','test'
  write(u,'(a,1x,a)')'output cif ','test'
  close(u)
 end subroutine output_gulp
 subroutine output_lammps()
  implicit none
  integer :: day,hour,i4_huge,milli,minute,month,second,year
  integer :: ii,jj
  character(len=10) :: time
  character(len=8)  :: date
  integer           :: u=222
  real              :: mmm,rrr
  integer           :: zzz
  character(len=2)  :: zlz
  character(len=80) :: DataFilename
  character(len=5)  :: extension=".data"
  character(len=4)  :: label(4)
  character(len=80) :: str
  DataFilename=filename(1:Clen_trim(filename))//extension
  DataFilename=adjustl(trim(DataFilename))
  forall (ii=1:4)
   forall (jj=1:4)
    label(ii)(jj:jj)=" "
   end forall
  end forall
  open(u,file=DataFilename)
  call date_and_time (date,time)
  read (date,'(i4,i2,i2)')year,month,day
  read (time,'(i2,i2,i2,1x,i3)')hour,minute,second,milli
  write(u,'(a,1x,i2,1x,i2,1x,i4,1x,i2,a1,i2,a1,i2)')'Created on',day,month,year,hour,':',minute,':',second
  write(u,*)' '
  write(u,*)n_atoms,' atoms'
  write(u,*)n_bonds,' bonds'
  write(u,*)n_bends,' angles' 
  write(u,*)n_torss,' dihedrals'
  write(u,*)n_imprs,' impropers'
  write(u,*)' '
  write(u,*)n_atom_types,' atom types'
  write(u,*)n_bond_types,' bond types'
  write(u,*)n_bend_types,' angle types'
  write(u,*)n_tors_types,' dihedral types'
  write(u,*)n_impr_types,' improper types'
  write(u,*)' '
  write(u,*)xlo_bound,xhi_bound,' xlo xhi'
  write(u,*)ylo_bound,yhi_bound,' ylo yhi'
  write(u,*)zlo_bound,zhi_bound,' zlo zhi'
  write(u,*)xy,xz,yz,' xy xz yz'
  write(u,*)' '
  write(u,'(a)')'Masses'
  write(u,*)' '
  do i=1,n_atom_types
   call CheckAtom(atom_types(i),mmm,rrr,zzz,zlz)
   write(u,'(i4,1x,f14.7,1x,a,a)')i,mmm,' # ',atom_types(i)
  end do
  write(u,*)' '
  write(u,'(a)')'Bond Coeffs'
  write(u,*)' '
  do i=1,n_bond_types
   read(bond_type_string(i),'(2a4)')(label(ii),ii=1,2)
   forall (ii=1:80)
    str(ii:ii)=" "
   end forall
   call search_forcefield(str,'bond',label)
   write(u,'(i4,3x,a,a,a)')i,adjustl(str(1:Clen_trim(str))),' # ',bond_type_string(i)
  end do
  write(u,*)' '
  write(u,'(a)')'Angle Coeffs'
  write(u,*)' '
  do i=1,n_bend_types
   read(bend_type_string(i),'(3a4)')(label(ii),ii=1,3)
   forall (ii=1:80)
    str(ii:ii)=" "
   end forall
   call search_forcefield(str,'bend',label)
   write(u,'(i4,3x,a,a,a)')i,adjustl(str(1:Clen_trim(str))),' # ',bend_type_string(i)
  end do
  write(u,*)' '
  write(u,'(a)')'Dihedral Coeffs'
  write(u,*)' '
  do i=1,n_tors_types
   read(tors_type_string(i),'(4a4)')(label(ii),ii=1,4)
   forall (ii=1:80)
    str(ii:ii)=" "
   end forall
   call search_forcefield(str,'tors',label)
   write(u,'(i4,3x,a,a,a)')i,adjustl(str(1:Clen_trim(str))),' # ',tors_type_string(i)
  end do
  write(u,*)' '
  !write(u,'(a)')'Improper Coeffs'
  !write(u,*)' '
  !do i=1,n_impr_types
  ! read(tors_type_string(i),'(4a4)')(label(ii),ii=1,4)
  ! forall (ii=1:80)
  !  str(ii:ii)=" "
  ! end forall
  ! call search_forcefield(str,'impr',label)
  ! write(u,'(i4,3x,a,a,a)')i,adjustl(str(1:Clen_trim(str))),' # ',impr_type_string(i)
  !end do
  !write(u,*)' '
  write(u,'(a)')'Pair Coeffs'
  write(u,*)' '
  do i=1,n_atom_types
   label(1)=atom_types(i)
   label(2)=atom_types(i)
   forall (ii=1:80)
    str(ii:ii)=" "
   end forall
   call search_forcefield(str,'pair',label)
   write(u,*)i,str,' # ', (label(ii),ii=1,2)
   !write(u,'(i4,3x,a,a,a)')i,str,' # ',(label(ii),ii=1,2)
  end do
  write(u,*)' '
  write(u,'(a)')'Atoms'
  write(u,*)' '
  do i=1,n_atoms
   write(u,'(i5,1x,i3,1x,i3,1x,f10.7,1x,3(f14.7,1x))')i,1,atom(i)%type_,atom(i)%charge,(atom(i)%xyzc(j,1),j=1,3)
  end do 
  write(u,*)' '
  write(u,'(a)')'Bonds'
  write(u,*)' '
  do i=1,n_bonds
   write(u,'(i5,a)')i,bond(i)
  end do
  write(u,*)' '
  write(u,'(a)')'Angles'
  write(u,*)' '
  do i=1,n_bends
   write(u,'(i5,a)')i,bend(i)
  end do
  write(u,*)' '
  write(u,'(a)')'Dihedrals'
  write(u,*)' '
  do i=1,n_torss
   write(u,'(i5,a)')i,tors(i)
  end do
  write(u,*)' '
  write(u,'(a)')'Impropers'
  write(u,*)' '
  do i=1,n_imprs
   write(u,'(i5,a)')i,impr(i)
  end do
  write(u,*)' '
  close(u)
 end subroutine output_lammps
 subroutine cellnormal2lammps(cell_0,xlo_bound,ylo_bound,zlo_bound,&
                              xhi_bound,yhi_bound,zhi_bound,xy,xz,yz)
  implicit none
  real,intent(in)  :: cell_0(6)
  real,intent(out) :: xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound
  real,intent(out) :: xy,xz,yz
  real,parameter :: pi=acos(-1.0)
  real,parameter :: radtodeg = 180.0/pi
  real,parameter :: degtorad = pi/180.0
  real :: xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz,cosa,cosb,cosg
  real :: alp,bet,gam,sing
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
  lx=cell_0(1)
  xy=cell_0(2)*cosg
  xz=cell_0(3)*cosb
  ly=sqrt(cell_0(2)*cell_0(2)-xy*xy)
  yz=(cell_0(2)*cell_0(3)*cosa-xy*xz)/ly
  lz=sqrt(cell_0(3)*cell_0(3)-xz*xz-yz*yz)
  xlo_bound=0.0
  ylo_bound=0.0
  zlo_bound=0.0
  xhi_bound=lx
  yhi_bound=ly
  zhi_bound=lz
  !xlo_bound = xlo + MIN(0.0,xy,xz,xy+xz)
  !xhi_bound = xhi + MAX(0.0,xy,xz,xy+xz)
  !ylo_bound = ylo + MIN(0.0,yz)
  !yhi_bound = yhi + MAX(0.0,yz)
  !zlo_bound = zlo
  !zhi_bound = zhi 
  return
 end subroutine cellnormal2lammps
 subroutine checkatom(Label,m,s,Z,Zlabel)
  implicit none
  character(len=4),intent(in)  :: Label
  real,intent(out)             :: m,s
  integer,intent(out)          :: Z
  character(len=2),intent(out) :: ZLabel
  select case(Label)
   case('C   ','C0  ':'C999')
    Z=6
    m=12.0107
    s=0.720
    ZLabel=' C'
   case('H   ','H0  ':'H999')
    Z=1
    m=1.00794
    s=0.320
    Zlabel=' H'
   case('N   ','N0  ':'N999')
    Z=7
    m=14.00674
    s=0.7
    Zlabel=' N'
   case('O   ','O0  ':'O999','OH  ')
    Z=8
    m=15.9994
    s=0.7
    Zlabel=' O'
   case('Si  ','Si0 ':'Si99')
    Z=14
    m=28.0855
    s=1.14
    zlabel='Si'
   case('Al  ')
    Z=13
    m=26.9800
    s=1.14
    zlabel='Al'
   case('Ar  ')
    Z=18
    m=39.948
    zlabel="Ar"
    s=0.0
   case('Xe  ')
    Z=54
    m=131.293
    s=0.0001
    Zlabel='Xe'
   case('Zn  ','Zn0 ':'Zn99')
    Z=30
    m=65.37
    s=1.6
    Zlabel='Zn'
   case('Cl  ',' Cl ','Cl0 ':'Cl99')
    Z=17
    m=35.453
    s=1.0
    Zlabel='Cl'
   case default
    write(6,'(a1,a4,a1)')"'",label,"'"
    STOP 'Atom unknowed'
  end select
 end subroutine checkatom

 PURE INTEGER FUNCTION Clen(s)      ! returns same result as LEN unless:
 CHARACTER(*),INTENT(IN) :: s       ! last non-blank char is null
 INTEGER :: i
 Clen = LEN(s)
 i = LEN_TRIM(s)
 IF (s(i:i) == CHAR(0)) Clen = i-1  ! len of C string
 END FUNCTION Clen
 PURE INTEGER FUNCTION Clen_trim(s) ! returns same result as LEN_TRIM unless:
 CHARACTER(*),INTENT(IN) :: s       ! last char non-blank is null, if true:
 INTEGER :: i                       ! then len of C string is returned, note:
                                    ! Ctrim is only user of this function
 i = LEN_TRIM(s) ; Clen_trim = i
 IF (s(i:i) == CHAR(0)) Clen_trim = Clen(s)   ! len of C string
 END FUNCTION Clen_trim
!
 SUBROUTINE cell(rv,vr,cell_0)
 implicit none
 integer :: i,j
 real, intent(in)  :: cell_0(6)
 real, intent(out) :: rv(3,3),vr(3,3)
 real, parameter   :: pi = ACOS(-1.0)
 real :: alp,bet
 real :: cosa,cosb,cosg
 real :: gam,sing
 real :: DEGTORAD
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
 WRITE(6,'(a)') 'Cell:'
 WRITE(6,'(6F14.7)')( cell_0(j), j=1,6 )
 WRITE(6,'(a)')'Linear Transformation Operator:'
 DO i=1,3
  WRITE(6,'(F14.7,F14.7,F14.7)')( rv(i,j), j=1,3 )
 ENDDO
 WRITE(6,'(a)')'----------------------------------------'
 WRITE(6,'(a)')'Inverse Linear Transformation Operator:'
 DO i=1,3
  WRITE(6,'(F14.7,F14.7,F14.7)')( vr(i,j), j=1,3 )
 ENDDO
 WRITE(6,'(a)')'----------------------------------------'
 RETURN
 END SUBROUTINE cell
!
 SUBROUTINE uncell(rv,cell_0)
  implicit none
  real,intent(out)   :: cell_0(6)
  real,intent(in)    :: rv(3,3)
  integer            :: i,j
  real               :: temp(6)
  REAL               :: radtodeg
  REAL, PARAMETER    :: pi=ACOS(-1.0) 
  radtodeg=180.0/PI
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
  DO i=4,6
     if (abs(cell_0(i) - 90.0 ).lt.0.00001) cell_0(i) = 90.0
     if (abs(cell_0(i) - 120.0).lt.0.00001) cell_0(i) = 120.0
  ENDDO
  RETURN
 END SUBROUTINE uncell
!
 SUBROUTINE inverse(a,c,n)
 implicit none
 integer n
 real a(n,n), c(n,n)
 real L(n,n), U(n,n), b(n), d(n), x(n)
 real coeff
 integer i, j, k
 L=0.0
 U=0.0
 b=0.0
 do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
 end do
 do i=1,n
  L(i,i) = 1.0
 end do
 do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
 end do
 do k=1,n
  b(k)=1.0
  d(1) = b(1)
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
 end do
 RETURN
 END SUBROUTINE inverse
!
subroutine make_dist_matrix(n,cell_0,rv,vr,x,dist_matrix)
 implicit none
 integer,intent(in) :: n
 real,intent(in)    :: cell_0(6),rv(3,3),vr(3,3),x(3,n)
 real,intent(out)   :: dist_matrix(n,n)
 integer            :: i,j,k
 real               :: r1(3),r2(3),s
 DO i=1,n
    dist_matrix(i,i)=0.0
    DO j=i+1,n
       forall ( k=1:3 )
        r1(k)=x(k,i)
        r2(k)=x(k,j)
       end forall
       call make_distances(cell_0,r1,r2,rv,s)
       dist_matrix(i,j)=s
       dist_matrix(j,i)=dist_matrix(i,j)
    END DO
 END DO
 return
end subroutine make_dist_matrix
!
 SUBROUTINE make_distances(cell_0,r2,r1,rv,dist)
 IMPLICIT NONE
 REAL,    intent(in)  :: r1(3),r2(3),rv(3,3),cell_0(6)    ! coordenadas y matriz de cambio
 REAL,    intent(out) :: dist
 REAL                 :: d_image(1:27),image(3,27)        ! array de distancias
 INTEGER              :: k,l,m,n,o,i,j                    ! variables mudas
 REAL                 :: atom(3),ouratom(3)               ! coordenadas preparadas
  k=0
  do l=-1,1
   do m=-1,1
      do n=-1,1
         k = k + 1
         ouratom(1) = r1(1)
         ouratom(2) = r1(2)
         ouratom(3) = r1(3)
         atom(1) = r2(1) + l
         atom(2) = r2(2) + m
         atom(3) = r2(3) + n
         d_image(k) = distance(atom,ouratom,rv)
         forall ( i=1:3)
           image(i,k) = atom(i)
         end forall
     enddo
   enddo
  enddo
  dist=MINVAL(d_image)
  RETURN
 END SUBROUTINE
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
  RETURN
end function
!
 REAL FUNCTION DISTANCE(atom,ouratom,rv)
  IMPLICIT NONE
  INTEGER :: j
  REAL :: atom(3),ouratom(3),per_atom(3),dist(3),o_atom(3),o_ouratom(3)
  REAL :: rv(3,3)
  FORALL ( j=1:3 )
   o_ouratom(j) = rv(j,1)*ouratom(1)  + rv(j,2)*ouratom(2)  + rv(j,3)*ouratom(3)
   o_atom(j)    = rv(j,1)*atom(1) + rv(j,2)*atom(2) + rv(j,3)*atom(3)
   dist(j) = o_ouratom(j) - o_atom(j)
  END FORALL
  DISTANCE = sqrt(dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3))
 END FUNCTION
!
 real function kcalmol2ev(x)
  implicit none
  real            :: x
  kcalmol2ev = 0.04336*x
  return
 end function  kcalmol2ev
 real function kjmol2ev(x)
  implicit none
  real,intent(in) :: x
  kjmol2ev=x*(1036.427/100000.0)
  return
 end  function kjmol2ev
 real function kjmolnmnm2evAA(x)
  real, intent(in) :: x
  kjmolnmnm2evaa = kjmol2ev(x)/100.0
  return
 end function  kjmolnmnm2evaa
 subroutine print_help()
    print '(a)', '  -h, --help   print usage information and exit'
    print '(a)', '  -c, --cif    CIF File input'
    print '(a)', '  -l, --linker [imi] imi, mimi'
    print '(a)', '  -GCMC, --adsorption [1e5] pressure Pa'
    print '(a)', '  -T, --temperatue  -T, [298] temperature in Kelvin'
    print '(a)', '  -fff, --forcefield [WHCJ] WHCJ, HZJ'
 end subroutine print_help
end program zif_cif2gin

