# assume mol.cml file in Angstroem exists in the same directory
# this bash command converts your mol.cml file to bohr
#

awk '{ {FS="\""}; {OFS="\""}; if ($1 ~ "<atom id") {print $1, $2, $3, $4, $5, $6/0.526, $7, $8/0.526, $9, $10/0.526, $11 } else print $0; }' <mol.cml  > mol-bohr.cml
