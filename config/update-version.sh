#!/bin/sh

test -f version || { echo "Not in top level source dir!"; exit 1; }

if [ "x$1" != "x" ]; then
	echo "$1" | grep '^[0-9]\+\.[0-9].*'
	[ $? != 0 ] && { echo "Invalid version format: $1"; exit 1;}
	echo $1 >version
fi

version=`cat version`
tmpfile=/tmp/update-version.$$ 

sed "s/AC_INIT(\[\(.*\)\],\[.*\],\[\(.*\)\])/\
AC_INIT([\\1],[$version],[\\2])/" configure.ac > $tmpfile
mv $tmpfile configure.ac
autoconf

sed "s/GIMIC_VERSION='.*'/GIMIC_VERSION='$version'/" globals.f90 >$tmpfile
mv $tmpfile globals.f90

sed "s/package_version='.*'/package_version='$version'/" gimic.in >$tmpfile
mv $tmpfile gimic.in

