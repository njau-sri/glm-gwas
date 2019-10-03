@echo off

set /p my_target=Specify Target Platform [win32/win64]:

rem Distributable Code for Visual Studio 2019
rem https://docs.microsoft.com/en-us/visualstudio/releases/2019/redistribution#visual-c-runtime-files
set "VisualStudioFolder=C:\Program Files (x86)\Microsoft Visual Studio\2019\Community"

set my_vcvar_bat=

if /i "%my_target%" == "win32" (
    rmdir /s /q Release glm-gwas-win32
    set "my_vcvar_bat=%VisualStudioFolder%\VC\Auxiliary\Build\vcvars32.bat"
    set "my_vcomp_dll=%VisualStudioFolder%\VC\Redist\MSVC\14.23.27820\x86\Microsoft.VC142.OPENMP\vcomp140.dll"
)

if /i "%my_target%" == "win64" (
    rmdir /s /q x64 glm-gwas-win64
    set "my_vcvar_bat=%VisualStudioFolder%\VC\Auxiliary\Build\vcvars64.bat"
    set "my_vcomp_dll=%VisualStudioFolder%\VC\Redist\MSVC\14.23.27820\x64\Microsoft.VC142.OPENMP\vcomp140.dll"
)

if "%my_vcvar_bat%" == "" (
    echo ERROR: invalid target platform: %my_target%
    pause
    exit /b 1
)

call "%my_vcvar_bat%"

MSBuild.exe -m -p:Configuration=Release

if errorlevel 1 (
    pause
    exit /b 1
)

if "%my_target%" == "win32" (
    mkdir glm-gwas-win32
    copy "%my_vcomp_dll%" glm-gwas-win32\
    copy Release\glm-gwas.exe glm-gwas-win32\
)

if "%my_target%" == "win64" (
    mkdir glm-gwas-win64
    copy "%my_vcomp_dll%" glm-gwas-win64\
    copy x64\Release\glm-gwas.exe glm-gwas-win64\
)

pause

set my_target=
set my_vcvar_bat=
set my_vcomp_dll=
set VisualStudioFolder=
