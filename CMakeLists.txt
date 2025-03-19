cmake_minimum_required(VERSION 3.5)

project(diffusion LANGUAGES Fortran)

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-wrap-margin -Ofast")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback")
  endif ()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -gdwarf -Og -fbacktrace -fcheck=all")
  endif ()
endif ()

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}")

get_directory_property(parent_dir PARENT_DIRECTORY)

set(sources
  src/module_tool.F90
  src/module_hdiff.F90
  src/module_vdiff.F90
  src/diffusion.F90
)

add_library(diffusion ${sources})

if (NOT parent_dir)
  add_executable(diffusion.exe app/main.f90)
  target_link_libraries(diffusion.exe diffusion)
endif ()
