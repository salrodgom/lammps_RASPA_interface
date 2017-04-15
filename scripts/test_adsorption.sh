#!/bin/bash
temperature=77
structure="sod_mim.cif"
for pressure in 1 10 100 1000 2000 2500 3000 3500 4000 5000 7500 10000 100000 ; do
 bash lammps_raspa_interface.sh $structure $temperature $pressure argon 5
done
exit 0 
