#/bin/bash
# ==================================================================
# Author:  Salvador Rodriguez-Gomez Balestra
# License: MIT License
# ==================================================================

# Inputs:  CIF, Temperature, Pressure, Guest_Molecule, #MCEMMDCycles
# ==================================================================
CIFFile=$1
structure=$(echo $CIFFile | sed 's/\.cif//g')
CIFTemporallyFile=${structure}_tmp.cif
temperature=$2
pressure=$3
guest=$4
MCCycles=1
MCEMMDCycles=$5
nCPU=16
# Files:
raspa_files_folder=$(pwd)/lib/fff_raspa
lammps_files_folder=$(pwd)/lib/fff_lammps
src_files_folder=$(pwd)/src


# Functions:
# ==========
function raspa {
 nohup simulate > /dev/null &
}

function lammps {
 nohup mpirun --np $((${nCPU}-2)) lmp_mpi -in in.lmp -sf opt &
}

function check_supercell {
# Check cell size for correct calculation of energies considering cutoff
 cutoff=12.0
 ua=1
 ub=1
 uc=1
 a_cell=$(grep "_cell_length_a" $CIFTemporallyFile | awk '{print $2}')
 b_cell=$(grep "_cell_length_b" $CIFTemporallyFile | awk '{print $2}')
 c_cell=$(grep "_cell_length_c" $CIFTemporallyFile | awk '{print $2}')
 a_cell=$(echo "$ua*$a_cell" | bc -l)
 while [ $(echo "$a_cell < 2*$cutoff" | bc -l) == 1 ] ; do
  let "ua++"
  a_cell=$(echo "$ua*$a_cell" | bc -l)
 done
 while [ $(echo "$b_cell < 2*$cutoff" | bc -l) == 1 ] ; do
  let "ub++"
  b_cell=$(echo "$ub*$b_cell" | bc -l)
 done
 while [ $(echo "$c_cell < 2*$cutoff" | bc -l) == 1 ] ; do
  let "uc++"
  c_cell=$(echo "$uc*$c_cell" | bc -l)
 done
 echo "Make Supercell: ${ua}x${ub}x${uc} > $a_cell $b_cell $c_cell"
}

function mc_muVT_raspa {
 # mu V T, ensemble
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
    check_supercell
    sed -i "s/SUPERCELL/$ua $ub $uc/g" simulation.input
    go_raspa
    wait_for_raspa
    sed '/MODEL    2/,$d' Movies/System_0/Movie*allcomponents.pdb > c
    sed 's/MODEL    2/MODEL    1/g' c > ../${CyclesNameFile}.pdb
    rm c
    rm -rf VTK Movies/System_0/Framework* Movies/System_0/Movie_*_component_*.pdb Movies/System_0/Movie_*_frameworks.pdb 
   cd ..
  fi
 done
}

function make_binaries {
# Make binaries
 cp ${src_files_folder}/*.f90 .
 cp ${src_files_folder}/Makefile .
 make install
}

function clean_binaries {
 make clean
}

function count_used_CPUs {
# 
 n_used=0
 for process in "simulate" "lmp_mpi" "gulp" ; do
  n=$(ps aux | grep ${process} | sed '/grep/d' | wc -l | awk '{print $1}')
  n_used=$((${n}+${n_used}))
 done
}

function go_raspa {
 count_used_CPUs
 while [ $(echo "${n_used} >= ${nCPU}" | bc -l) == 1 ] ; do
  sleep 30
  count_used_CPUs
 done
 raspa
 sleep 1
}

function wait_for_raspa {
 raspa_end=$(grep "Simulation finished" Output/System_0/output_*.data | wc -l | awk '{print $1}')
 while [ $(echo $raspa_end < 1 | bc -l) == 1 ] ; do
  sleep 30
  raspa_end=$(grep "Simulation finished" Output/System_0/output_*.data | wc -l | awk '{print $1}')
 done
}

function go_lammps {
 count_used_CPUs
 while [ $(echo "${n_used} >= ${nCPU}" | bc -l) == 1 ] ; do
  sleep 30
  count_used_CPUs
 done
 lammps
 sleep 1
}

function wait_for_lammps {
 lammps_end=$( grep "Total wall time:" logs/log.main | wc -l | awk '{print $1}')
 while [ $(echo ${lammps_end} < 1 | bc -l) == 1 ] ; do
  sleep 30
  lammps_end=$( grep "Total wall time:" logs/log.main | wc -l | awk '{print $1}')
 done
}

function em_md_lammps {
 # Optimization
 # NVE
 # NPTPR
 folder=${CyclesNameFile}_emmd
 if [ ! -d $folder ] ; then 
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
   go_lammps
   wait_for_lammps
  cd ..
 fi
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
 # Remove guest!!!
 sed -i '/ Ar /d' input.pdb
 #
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
  mc_muVT_raspa
  raspa_lammps
  em_md_lammps
  lammps_raspa
  clean_binaries
done
exit 0
