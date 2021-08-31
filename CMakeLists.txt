# CMakeList.txt : CMake project for dcproj, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("dcproj")
#target_link_libraries (dcproj Eien3::Eigen)
# Add source to this project's executable.
add_executable (dcproj
src/beam3dtools.cpp src/beam3dtools.h
src/dcproj.cpp src/dcproj.h
src/mesh.cpp src/mesh.h
src/shapequad.cpp src/shapequad.h
src/shapehexahedron.cpp src/shapehexahedron.h
src/slopeproject.cpp src/slopeproject.h 
src/vonmises.cpp src/vonmises.h
src/druckerprager.cpp src/druckerprager.h
src/mohrcoulomb.cpp src/mohrcoulomb.h
src/postprocess.cpp src/postprocess.h 
src/material.cpp  src/material.h 
src/KLGalerkinRF.cpp src/KLGalerkinRF.h
src/gridmesh.cpp src/gridmesh.h
src/elastoplasticbase.cpp src/elastoplasticbase.h
src/elastoplastic2D.cpp src/elastoplastic2D.h
src/elastoplastic3D.cpp src/elastoplastic3D.h
src/elastmat2D.cpp src/elastmat2D.h
src/eigen_unsym.cpp src/eigen_unsym.h
src/error.cpp src/error.h
src/pressurizedhole.cpp src/pressurizedhole.h
src/vtkmesh.cpp src/vtkmesh.h
src/eigen_sym.h src/cholesky.h src/nr3.h  src/roots.h src/mins_ndim.h  src/ludcmp.h
)
#find_package (Eigen3 3.3.9 REQUIRED NO_MODULE)
 
include_directories('src')
include_directories("/usr/include/eigen3")
#include_directories("/opt/intel/mkl")
#include_directories("/opt/intel")
#include_directories("/opt/intel/oneapi/mkl/2021.3.0/include")
#include_directories(${MKLROOT})
#include_directories("/opt/intel/oneapi/mkl/2021.3.0/include" "/opt/intel/oneapi/mkl/2021.3.0/lib/intel64/" "/usr/include/mkl" "/usr/lib/x86_64-linux-gnu/mkl")
#include_directories( ${Boost_INCLUDE_DIR} ${MKLROOT} "/opt/intel/mkl")

target_include_directories(dcproj PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(dcproj PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_include_directories(dcproj PUBLIC ${PROJECT_SOURCE_DIR}/eigen)
target_include_directories(dcproj PUBLIC ${PROJECT_SOURCE_DIR}/eigen-3.3.9)


set(THREADS_PREFER_PTHREAD_FLAG ON)
set(HEADER_FILES ${PROJECT_SOURCE_DIR}/src)
#set(MKLROOT "/opt/intel/oneapi/mkl/2021.3.0")

find_package(Threads REQUIRED)
#find_package(MKL REQUIRED)

#set(USING_MKL ON)

#TARGET_LINK_LIBRARIES( dcproj LINK_PUBLIC     ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_cdft_core.a ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl)

TARGET_LINK_LIBRARIES( dcproj LINK_PUBLIC ${Boost_LIBRARIES} Threads::Threads ${M_LIBRARIES} )

find_path(MKL_INCLUDE NAMES mkl.h PATHS ${SEARCH_DIRS} /opt/intel/mkl/include /softwares/intel/mkl/include)
		find_library(MKL_LIB_INTEL NAMES libmkl_intel_lp64.so libmkl_intel_lp64.dylib PATHS
				${SEARCH_DIRS} /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64)
		find_library(MKL_LIB_CORE NAMES libmkl_core.so libmkl_core.dylib PATHS ${SEARCH_DIRS}
				/opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64/)
		find_library(MKL_LIB_THREAD NAMES libmkl_intel_thread.so libmkl_intel_thread.dylib PATHS
				${SEARCH_DIRS} /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64/)
		find_library(COMPOSER_OMP NAMES libiomp5.so libiomp5.dylib PATHS ${SEARCH_DIRS}
				/opt/intel/composer_xe/compiler/lib
				/opt/intel/composerxe/lib/intel64
				/opt/intel/lib /softwares/intel/lib/intel64
                                /opt/intel/compilers_and_libraries/linux/lib/intel64_lin )

                link_libraries(${MKL_LIB_INTEL})
                link_libraries(${MKL_LIB_CORE})
                link_libraries(${MKL_LIB_THREAD})
                link_libraries(${COMPOSER_OMP})


add_compile_options( -DMKL_ILP64  -m64  -I"/opt/intel/oneapi/mkl/2021.3.0/include")








