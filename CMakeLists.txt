cmake_minimum_required(VERSION 3.10)
project(GbSMTWTP VERSION 1.0.0 LANGUAGES CXX)

# Configuración de compilación
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Flags de compilación para optimización y rendimiento
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -funroll-loops -ffast-math -flto")

# Crear la biblioteca header-only
add_library(GbSMTWTP INTERFACE)

# Configurar rutas de inclusión para la instalación
target_include_directories(GbSMTWTP INTERFACE
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Opciones de compilación para quien use la biblioteca
target_compile_options(GbSMTWTP INTERFACE -O3 -funroll-loops -ffast-math -flto)

# Exportar la configuración de la biblioteca
install(TARGETS GbSMTWTP EXPORT GbSMTWTPConfig)
install(EXPORT GbSMTWTPConfig
    FILE GbSMTWTPConfig.cmake
    DESTINATION lib/cmake/GbSMTWTP
)

# Instalar los archivos de cabecera en /usr/local/include/gb
install(DIRECTORY include/ DESTINATION include/gb)