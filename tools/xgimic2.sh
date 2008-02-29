#!/bin/sh
#

ACES2=/opt/aces2

PATH=$ACES2/bin:$PATH
GENBAS=$ACES2/basis/GENBAS

set -e
#ulimit -s unlimited
#ulimit -d unlimited
#ulimit -m unlimited
#ulimit -c 0

cleaner () {
	for i in FILES GRD JOBARC NEWMOS MOABCD GAMLAM \
		IIII IIJJ IJIJ IJKL\
		MOINTS OPTARC JAINDX MOLDEN VMLSYM DERGAM  \
		PPPHAA PPHHAA PHPHAA HHHHAA \
		DERINT VMSCR1 HF2 DTRIPLES LAMBDA3 MOLECULE.INP DLAMBDA3 TRIPLES; do
		[ -f $i ] && (rm -f $i; echo $i)
	done
	rm -f fort.* 
}

case "$1" in
	--xjoda)
		cleaner
		if [ ! -f ./GENBAS ]; then
			cp $GENBAS .
		fi	
		echo "*** xjoda ***"; xjoda 
	;;
	--clean)
		cleaner
	;;
	--scf)
		if [ ! -f ./GENBAS ]; then
			cp $GENBAS .
		fi	
		cleaner
		echo "*** xjoda ***"; xjoda 
		xvmol
		xvmol2ja
		xvscf
		xvtran
		xintprc
		xvdint
		xcphf
		xjoda
		xcpdens
		cleaner
	;;
	--mbpt)
		if [ ! -f ./GENBAS ]; then
			cp $GENBAS .
		fi	
		cleaner
		echo "*** xjoda ***"; xjoda 
		echo "*** xvmol"; xvmol
		echo "*** xvmol2ja"; xvmol2ja
		echo "*** xvscf"; xvscf
		echo "*** xvtran"; xvtran
		echo "*** xintprc"; xintprc
		echo "*** xvcc"; xvcc
		echo "*** xlambda"; xlambda
		echo "*** xdens"; xdens
		echo "*** xvdint"; xvdint
		echo "*** xcphf"; xcphf
		echo "*** xsdcc"; xsdcc
		echo "*** xvdint"; xvdint
		echo "*** xsdcc"; xsdcc
		echo "*** xvdint"; xvdint
		echo "*** xsdcc"; xsdcc
		echo "*** xjoda"; xjoda
		xcpdens
		cleaner
	;;
	--cc)
		if [ ! -f ./GENBAS ]; then
			cp $GENBAS .
		fi	
		cleaner
		echo "*** xjoda ***"; xjoda 
		echo "*** xvmol"; xvmol
		echo "*** xvmol2ja"; xvmol2ja
		echo "*** xvscf"; xvscf
		echo "*** xvtran"; xvtran
		echo "*** xintprc"; xintprc
		echo "*** xecc"; xecc
		echo "*** xlcc"; xlcc
		echo "*** xdens"; xdens
		echo "*** xvdint"; xvdint
		echo "*** xcphf"; xcphf
		echo "*** xsdcc"; xsdcc
		echo "*** xvdint"; xvdint
		echo "*** xsdcc"; xsdcc
		echo "*** xvdint"; xvdint
		echo "*** xsdcc"; xsdcc
		echo "*** xjoda"; xjoda
		xcpdens
		cleaner
	;;
	*)
		echo "usage: jaces2 {--xjoda|--scf|--mbpt|--cc|--clean}"
	;;
esac

