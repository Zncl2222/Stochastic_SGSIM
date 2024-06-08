@echo off

mkdir cbuild
cd cbuild

IF /I "%1"=="-s" (
  cmake -DIS_EXECUTE=OFF -G "MinGW Makefiles" ..
  make
  %~dp0cbuild/c_example.exe
) ELSE (
  cmake -DIS_EXECUTE=ON -G "MinGW Makefiles" ..
)

cd ..

exit /b 0
