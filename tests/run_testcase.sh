#!/bin/bash
echo "====== Delete previous results ======"

rm TestCases_Results/*.dat

echo "====== First run from day 35001 until day 35702 ======"
awk '{sub(/35732/,"35702")}; 1' TestCase_LakeZurich.par > TestCase_LakeZurich_1.par
diff TestCase_LakeZurich.par TestCase_LakeZurich_1.par
time ../build/simstrat TestCase_LakeZurich_1.par
rm TestCase_LakeZurich_1.par

echo "====== Second run from day 35702 until day 35732 ======"
time ../build/simstrat TestCase_LakeZurich.par

diff -x .gitignore -x simulation-snapshot.dat TestCases_Results_expected/ TestCases_Results/ > diff_results.txt
if [[ $? == 0 ]]; then
    echo "====== Successful test ======"
else
    echo "====== Test failed ====== see diff_results.txt for more details"
fi
