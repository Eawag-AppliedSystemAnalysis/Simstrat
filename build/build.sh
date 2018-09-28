#!/bin/bash

# by Davide Vanzo (davide.vanzo@eawag.ch), 2018


cd $(dirname $0)

usage() {
    str=$'\nUsage: '
    str+=$0
    str+=$' [-m <string>] [-c] [-d]\n\n'
    str+=$'  -m <string>  compile mode: release (default), debug\n'
    str+=$'  -c           clean previous compiling files before building\n'
    str+=$'  -d           generate code documentation with FORD after building\n'
    str+=$' in case of no arguments, "./build.sh -m release" is called\n'
    echo "$str" 1>&2
    exit 1
}

printHeadFoot() {
    N=60
    text=$1
    size=${#text}
    if [ $size -gt $N ]; then echo "string too long!"; exit; fi
    diff1=$((N-size-2))
    diff2=$((diff1/2))
    diff3=$((diff2+diff1%2))
    printf_new "*" $N "\n*"
    printf_new " " $diff2
    printf "$text"
    printf_new " " $diff3 "*\n"
    printf_new "*" $N "\n"
}

printf_new() {
    str=$1
    num=$2
    v=$(printf "%-${num}s" "$str")
    printf "${v// /$str}"$3
}

printErrorAndExit() {
    printHeadFoot "$1"
    exit 0
}

cleanup() {
    printHeadFoot "cleaning build/bin directories"

    find ../bin -type f \( -iname "*" ! -iname "README.md" \) -delete
    find . -type f \( -iname "*" ! -iname "README.md" ! -iname "build.sh" ! -iname "fobos" \) -delete
    find . -type d -exec rm -r "{}" \;
    find ../lib/csv_fortran -name "build" -type d -exec rm -r "{}" \;
    find ../lib/json_fortran -name "build" -type d -exec rm -r "{}" \;

    printHeadFoot "successfully cleaned build/bin dirs"
}


docuFORD() {
    printHeadFoot "Building documentation with FORD..."

    cd ../doc/developer/ford

    if hash ford 2>/dev/null; then
        ford ford_projectfile.md
    else
        printErrorAndExit "Ford not found! Cannot build documentation"
    fi

    # go back
    cd -
    printHeadFoot "successful documentation generation"
}


# default values
modeType=release
docu=0
cleaning=0

# parse command line options
optCounter=0
while getopts ":m:dc" o; do
    optCounter=$((optCounter+1))
    case "${o}" in
        m)
            modeType=${OPTARG}
            ;;
        d)
            docu=1
            ;;
        c)
            cleaning=1
            ;;
        *)
            usage
            ;;
    esac
done

# eventually clean up
if [ $cleaning -eq 1 ]; then
    cleanup
fi


#build using FoBiS:
printHeadFoot "Building simstrat and dependencies with Fobis..."
if hash FoBiS.py 2>/dev/null; then
    FoBiS.py build -mode "$modeType"-gnu
else
    printErrorAndExit "FoBiS.py not found! Cannot build library"
fi

# move the binaries
{
    mv simstrat ../bin/. >/dev/null 2>&1
} || {
    printErrorAndExit "unsuccessful executable generation"
}

# eventually generate the documentation with FORD
if [ $docu -eq 1 ]; then
    docuFORD
fi

printHeadFoot "simstrat is ready... what a wonderful day!"

