REM Start all slaves
START /B cmd /c echo model.bat^> "%temp%\answer.tmp" ^& (pslave ^< "%temp%\answer.tmp") ^& del "%temp%\answer.tmp"
cd cpu2
START /B cmd /c echo model.bat^> "%temp%\answer.tmp" ^& (pslave ^< "%temp%\answer.tmp") ^& del "%temp%\answer.tmp"
cd ..
cd cpu3
START /B cmd /c echo model.bat^> "%temp%\answer.tmp" ^& (pslave ^< "%temp%\answer.tmp") ^& del "%temp%\answer.tmp"
cd ..
cd cpu4
START /B cmd /c echo model.bat^> "%temp%\answer.tmp" ^& (pslave ^< "%temp%\answer.tmp") ^& del "%temp%\answer.tmp"
cd ..
REM Launch PEST in parallel mode
ppest keps_calib
REM Clear PEST result files
del "results\*.dat" 
del "cpu2\results\*.dat" 
del "cpu3\results\*.dat" 
del "cpu4\results\*.dat" 
