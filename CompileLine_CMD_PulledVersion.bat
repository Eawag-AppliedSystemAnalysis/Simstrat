cls
cd D:\InstaledPrograms\Simstrat\Simstrat\Simstrat_v2

rmdir /S /Q bin

python D:\InstaledPrograms\Phython\Scripts\FoBis.py clean

python D:\InstaledPrograms\Phython\Scripts\FoBis.py build

DEL "D:\InstaledPrograms\Simstrat\TestIce\simstrat.exe"

xcopy /y "D:\InstaledPrograms\Simstrat\Simstrat\Simstrat_v2\bin\simstrat.exe" "D:\InstaledPrograms\Simstrat\TestIce\"