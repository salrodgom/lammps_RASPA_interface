install:
	gfortran cif2lammps.f90 -o cif2lammps
	gfortran lammpstrj2pdb.f90  -o lammpstrj2pdb
	gfortran pdb2cif.f90  -o pdb2cif
clean:
	rm cif2lammps lammpstrj2pdb pdb2cif
