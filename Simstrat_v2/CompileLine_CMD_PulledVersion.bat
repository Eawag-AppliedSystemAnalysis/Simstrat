REM clear comand window
cls
cd D:\InstaledPrograms\Simstrat\Simstrat\Simstrat_v2
REM Deleats the bin directory
rmdir /S /Q bin
REM clean upp previus builds
python D:\InstaledPrograms\Phython\Scripts\FoBis.py clean
REM compile simstrat
python D:\InstaledPrograms\Phython\Scripts\FoBis.py build
REM Delet old exe file
DEL "D:\InstaledPrograms\Simstrat\TestIce\simstrat.exe"
REM move Simstrat execution file
xcopy /y "D:\InstaledPrograms\Simstrat\Simstrat\Simstrat_v2\bin\simstrat.exe" "D:\InstaledPrograms\Simstrat\TestIce\"