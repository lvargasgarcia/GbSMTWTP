#!/bin/bash
set -e

# Crear el directorio build si no existe
mkdir -p build
cd build

# Configurar el proyecto con CMake
cmake .. -DCMAKE_BUILD_TYPE=Release

# Compilar e instalar (requiere permisos de superusuario)
sudo make install