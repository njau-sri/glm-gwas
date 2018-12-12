
rmdir /s /q Release

set REDIST=C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Redist\MSVC\14.16.27012

set VCOMPDLL=%REDIST%\x86\Microsoft.VC141.OPENMP\vcomp140.dll

call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvars32.bat"

MSBuild.exe /p:Configuration=Release /p:Platform=Win32

cd Release
mkdir glm-gwas-win32
copy glm-gwas.exe glm-gwas-win32\
copy "%VCOMPDLL%" glm-gwas-win32\

pause
