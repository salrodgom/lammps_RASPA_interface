#/bin/bash -x
# Inputs:
CIFFile=$1
structure=$(echo $CIFFile | sed 's/\.cif//g')
CIFTemporallyFile=${structure}_tmp.cif
temperature=$2
pressure=$3
guest=$4
MCCycles=1
MCEMMDCycles=$5
nCPU=3
# Files:
raspa_files_folder=$(pwd)/lib/fff_raspa
lammps_files_folder=$(pwd)/lib/fff_lammps
src_files_folder=$(pwd)/src
# Functions:
function mc_raspa {
 # mu V T
 # RASPA
 for i in $(seq 1 ${MCCycles}) ; do
  iname=$(echo $i | awk '{ printf("%02d\n", $1) }')
  NewNameFile=${structure}_${guest}_${pressure}_${temperature}
  CyclesNameFile=${cycle_name}_${iname}_${NewNameFile}
  folder=${CyclesNameFile}_adsorption
  if [ ! -d $folder ] ; then
   mkdir $folder
   cp ${raspa_files_folder}/*.def $folder
   cp ${raspa_files_folder}/INPUT ${folder}
   cp $CIFTemporallyFile $folder
   cd $folder
    sed "s/STRUCTURE/${structure}_tmp/g" INPUT > simulation.input
    sed -i "s/TEMPERATURE/${temperature}/g" simulation.input
    sed -i "s/PRESSURE/${pressure}/g" simulation.input
    sed -i "s/GUEST/${guest}/g" simulation.input
    sed -i "s/SUPERCELL/2 2 2/g" simulation.input
    simulate > /dev/null
    sed '/MODEL    2/,$d' Movies/System_0/Movie*allcomponents.pdb > c
    sed 's/MODEL    2/MODEL    1/g' c > ../${CyclesNameFile}.pdb
    rm c
    rm -rf VTK Movies/System_0/Framework* Movies/System_0/Movie_*_component_*.pdb Movies/System_0/Movie_*_frameworks.pdb 
   cd ..
  fi
 done
}
function make_binaries {
 cp ${src_files_folder}/*.f90 .
 cp ${src_files_folder}/Makefile .
 make install
}
function clean_binaries {
 make clean
}
function em_md_lammps {
 # Optimization
 # NVE
 # NPTPR
 folder=${CyclesNameFile}_emmd
 mkdir $folder
 cp ${CyclesNameFile}.lmp $folder/.
 cp ${lammps_files_folder}/in.lmp $folder/in.lmp
 cp atom_types_for_dump.txt $folder/.
 cd $folder
  sed -i "s/TEMPERATURE/${temperature}/g" in.lmp
  sed -i "s/PRESSURE/${pressure}/g" in.lmp
  sed -i "s/FILENAME/${CyclesNameFile}/g" in.lmp
  elements=$(cat atom_types_for_dump.txt | sed 's/[0-9]//g' | sed 's/  / /g')
  sed -i "s/ELEMENTS/${elements}/g" in.lmp
  mpirun --np $nCPU lmp_mpi -in in.lmp -sf opt
 cd ..
}
function raspa_lammps {
 cp ${CyclesNameFile}.pdb input.pdb
 ./pdb2cif
 mv p1.cif ${CyclesNameFile}.cif
 cp ${CyclesNameFile}.cif ${CIFTemporallyFile}
 cp lib/forcefield.lib .
 ./cif2lammps -c ${CyclesNameFile}.cif -wq -S
 mv ${CyclesNameFile}.data ${CyclesNameFile}.lmp
}
function lammps_raspa {
 ./lammpstrj2pdb < $folder/movs/all.lammpstrj
 n_lines=$(wc -l out.pdb |awk '{print $1}')
 line=$(sed -n '/MODEL/{=;p}' out.pdb | sed '{N;s/\n/ /}' | tail -n1 | awk '{print $1}')
 end=$((n_lines - line + 1))
 tail -n$end out.pdb > input.pdb
 ./pdb2cif 
 mv p1.cif ${CyclesNameFile}.cif
 cp ${CyclesNameFile}.cif ${CIFTemporallyFile}
}
##############################################################
# main program:
##############################################################
cp ${CIFFile} ${CIFTemporallyFile}
for cycle in $(seq 1 ${MCEMMDCycles}) ; do
  cycle_name=$(echo $cycle | awk '{ printf("%02d\n", $1) }')
  make_binaries
  mc_raspa
  raspa_lammps
  em_md_lammps
  lammps_raspa
  clean_binaries
done
exit 0
