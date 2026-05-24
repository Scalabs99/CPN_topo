#!/bin/bash

if [ "$#" -gt 2 ] || [ "$#" -eq 0 ]; then
	echo "Wrong number of arguments! Correct usage: $0 CPN_NCOLORS (TARGET_EXECUTABLE) If no target is given, all targets are compiled"
	exit 1
fi

if [ "$#" -eq 2 ]; then
	CPN_NCOLORS=$1
	TARGET_EXEC=$2
	echo "Compiling $2 with CPN_NCOLORS=$1"
fi

if [ "$#" -eq 1 ]; then
	CPN_NCOLORS=$1
	TARGET_EXEC=''
	echo "Compiling all src files with CPN_NCOLORS=$1"
fi

# change value of N in include/macro.h
aux=$( grep '#define N' include/macro.h )
sed "s/${aux}/#define N ${CPN_NCOLORS}/g" include/macro.h > temp
mv temp include/macro.h

if [[ "$OSTYPE" == "darwin"* ]]; then
    OPENSSL_DIR=$(brew --prefix openssl 2>/dev/null)
    if [ -d "$OPENSSL_DIR" ]; then
        MY_CFLAGS="-g -O0 -I${OPENSSL_DIR}/include -Wno-deprecated-declarations"
        MY_LDFLAGS="-L${OPENSSL_DIR}/lib"
    else
        echo "Warning: OpenSSL non trovato tramite Homebrew. Procedo coi parametri di base."
        MY_CFLAGS="-g -O0 -Wno-deprecated-declarations"
        MY_LDFLAGS=""
    fi
else
    MY_CFLAGS="-g -O0 -Wno-deprecated-declarations"
    MY_LDFLAGS=""
fi

# compile code
./configure CFLAGS="$MY_CFLAGS" LDFLAGS="$MY_LDFLAGS"
make ${TARGET_EXEC}
