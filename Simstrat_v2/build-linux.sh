#!/bin/bash
#  Simple build script.
#
#  Requires: FoBiS and Ford
#

#build using FoBiS:

if hash FoBiS.py 2>/dev/null; then

    echo "Building simstrat and dependencies..."

    FoBiS.py build


    #echo "Building test programs..."

    #FoBiS.py build -compiler ${FCOMPILER} -cflags "${FCOMPILERFLAGS}" -dbld ${BINDIR} -s ${TESTSRCDIR} -dmod ./ -dobj ./ -colors -libs ${LIBDIR}${LIBOUT} --include ${LIBDIR}

else
    echo "FoBiS.py not found! Cannot build library. Install using: sudo pip install FoBiS.py"
fi

# build the documentation using FORD:

#if hash ford 2>/dev/null; then
#
#    echo "Building documentation..."
#
#    ford ${FORDMD}
#
#else
#    echo "Ford not found! Cannot build documentation. Install using: sudo pip install ford"
#fi
