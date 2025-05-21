# Proyecto-TTI---24-25

##  Compilación y ejecución de la aplicación principal
g++ tests/EFK_GEOS3.cpp src/*cpp -lm -std=c++23 -o bin/main.exe
cd bin
main.exe
pause

##  Compilación y ejecución de los tests unitarios
g++ tests/tests.cpp src/*.cpp -lm -std=c++23 -o bin/tests.exe
cd bin
tests.exe
pause