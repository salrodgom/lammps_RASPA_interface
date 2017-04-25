#!/bin/bash
temperature=85
structure="ZIF-8.cif"
for pressure in 55000 17500 1400 60000 10 ; do
 # 55000 17500 1400 60000 10
 bash lammps_raspa_interface.sh $structure $temperature $pressure argon 5
done
exit 0 
