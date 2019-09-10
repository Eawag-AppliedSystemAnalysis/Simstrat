#!/bin/bash
rm TestCases_Results/*.dat

../build/simstrat TestCase_LakeZurich.par

diff -x .gitignore TestCases_Results_expected/ TestCases_Results/ > diff_results.txt
if [[ $? == 0 ]]; then
    echo "====== Successful test ======"
else
    echo "====== Test failed ====== see diff_results.txt for more details"
fi
