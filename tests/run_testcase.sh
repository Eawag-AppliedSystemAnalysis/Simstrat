#!/bin/bash
run()
{
    local end_time=$1
    echo "====== Run until day $end_time ======"
    
    awk  -v end_time=$end_time '$1 <= end_time || $1 == "t"' TestCase_LakeZurich/forcing_lake_zurich_1981_2013_waedenswil_homogenized.dat3 > TestCase_LakeZurich/forcing.dat
    awk  -v end_time=$end_time '{sub(/35732/,end_time)}; 1' TestCase_LakeZurich.par \
        | awk '{sub(/forcing_lake_zurich_1981_2013_waedenswil_homogenized.dat3/,"forcing.dat")}; 1' \
        > TestCase_LakeZurich_1.par
    time ../build/simstrat TestCase_LakeZurich_1.par
    tail -n1 TestCase_LakeZurich/forcing.dat
    tail -n1 TestCases_Results/HK_out.dat
    rm TestCase_LakeZurich_1.par
    rm TestCase_LakeZurich/forcing.dat
}

echo "====== Delete previous results ======"

rm TestCases_Results/*.dat

run 35150
run 35250.125
run 35350.375

echo "====== Run until day 35732 ======"
time ../build/simstrat TestCase_LakeZurich.par

diff -x .gitignore -x simulation-snapshot.dat TestCases_Results_expected/ TestCases_Results/ > diff_results.txt
if [[ $? == 0 ]]; then
    echo "====== Successful test ======"
else
    echo "====== Test failed ====== see diff_results.txt for more details"
fi
