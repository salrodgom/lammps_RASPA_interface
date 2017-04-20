#!/bin/bash
temperature=85
structure="ZIF-8.cif"
for pressure in 10000 ; do
 # 10000 17500 55000 60000 10
 bash lammps_raspa_interface.sh $structure $temperature $pressure argon 2
done
exit 0 
