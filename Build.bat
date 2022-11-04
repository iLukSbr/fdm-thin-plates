@CHCP 65001 >NUL

@echo off

x86_64-w64-mingw32-gfortran.exe PLACAS.f90 -o  PLACAS.exe -std=f2018 -fdefault-real-8 -fcheck=all liblapack.dll

pause