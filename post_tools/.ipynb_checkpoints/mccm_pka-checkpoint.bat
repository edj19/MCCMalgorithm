set "ifor1=C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2018.5.274"

set "ifor2=C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2018.3.210"

if exist "%ifor1%" (set "ifor=%ifor1%")

if exist "%ifor2%" (set "ifor=%ifor2%")

set "PATH=%ifor%\windows\bin;%PATH%"
ifortvars.bat intel64 vs2017 && f2py -m --overwrite-signature pka_number -h pka_number.pyf  pka_number.f90 && f2py -c --fcompiler=intelvem  pka_number.pyf  pka_number.f90

REM X "ifortvars.bat" intel64 vs2017 
REM X set "PATH=%ifor%\windows\bin;PATH"
set "PATH=C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2018.5.274\windows\bin;%PATH%"