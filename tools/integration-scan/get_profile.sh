#!/bin/sh

grep -i 'Induced current (nA/T)' gimic_h0* | awk '{printf "%s %20.9f\n", $1,$6}' | sed "s/gimic_h0/ /g" | sed "s/:/ /g" | sort -n > c_profile.out
grep -i 'Positive contribution' gimic_h0* | awk '{printf "%s %20.9f\n", $1,$6}' | sed "s/gimic_h0/ /g" | sed "s/:/ /g" | sed -n 'g;n;p' | sort -n > positive.out
grep -i 'Negative contribution' gimic_h0* | awk '{printf "%s %20.9f\n", $1,$6}' | sed "s/gimic_h0/ /g" | sed "s/:/ /g" | sed -n 'g;n;p' | sort -n > negative.out
# that would work for modulus:
# grep -i 'Positive contribution' gimic_h0* | awk '{printf "%s %20.9f\n", $1,$6}' | sed "s/gimic_h0/ /g" | sed "s/:/ /g" | sed -n 'p;N' | sort -n > positive.out
# grep -i 'Negative contribution' gimic_h0* | awk '{printf "%s %20.9f\n", $1,$6}' | sed "s/gimic_h0/ /g" | sed "s/:/ /g" | sed -n 'p;N' | sort -n > negative.out

# assuming a step size of 0.1 bohr gnuslice.cmd gives dF_j(x)/dx
gnuplot /home/heike/bin/my_bin/gnuslice.cmd
