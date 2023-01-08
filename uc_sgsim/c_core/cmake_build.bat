mkdir build
cd build
cmake ..
cd ..
cmake --build build/ -- /p:Configuration=Release
move %~dp0build\Release\uc_sgsim.dll %~dp0
