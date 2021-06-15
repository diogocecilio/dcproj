# CMakeList.txt : CMake project for dcproj, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("dcproj")

# Add source to this project's executable.
add_executable (dcproj
"dcproj.cpp" "dcproj.h" 
src/mesh.cpp src/mesh.h
src/shapequad.cpp src/shapequad.h
src/slopeproject.cpp src/slopeproject.h 
src/vonmises.cpp src/vonmises.h
src/druckerprager.cpp src/druckerprager.h
src/postprocess.cpp src/postprocess.h 
src/material.cpp  src/material.h 
src/KLGalerkinRF.cpp src/KLGalerkinRF.h
src/gridmesh.cpp src/gridmesh.h
src/elastoplasticbase.cpp src/elastoplasticbase.h
src/elastoplastic2D.cpp src/elastoplastic2D.h
src/elastmat2D.cpp src/elastmat2D.h
src/eigen_unsym.cpp src/eigen_unsym.h

src/eigen_sym.h src/cholesky.h src/nr3.h  src/roots.h src/mins_ndim.h  src/ludcmp.h
)

include_directories('src')

target_include_directories(dcproj PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(dcproj PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_include_directories(dcproj PUBLIC ${PROJECT_SOURCE_DIR}/eigen)

set(HEADER_FILES ${PROJECT_SOURCE_DIR}/src)

#add_executable(dcproj execfile.cpp ${HEADER_FILES})

#add_library(mylib libsrc.cpp ${HEADER_FILES})
#target_include_directories(mylib PRIVATE ${YOUR_DIRECTORY})
#add_executable(myexec execfile.cpp ${HEADER_FILES})
#target_include_directories(myexec PRIVATE ${YOUR_DIRECTORY})



#set(INCLUDE_DIRECTORIES src/include)

# TODO: Add tests and install targets if needed.
