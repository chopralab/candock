set(gsl_source  "${CMAKE_CURRENT_BINARY_DIR}/gsl-src")
set(gsl_install "${CMAKE_BINARY_DIR}/stage")
set(gsl_version "2.5.0")
set(gsl_url "https://github.com/ampl/gsl/archive/v${gsl_version}.tar.gz")
set(gsl_md5 "78419573df6db3186c651ac9ca3d5d38")

ExternalProject_Add(gsl
  DOWNLOAD_DIR ${download_dir}
  SOURCE_DIR "${gsl_source}"
  INSTALL_DIR "${gsl_install}"
  URL ${gsl_url}
  URL_MD5 ${gsl_md5}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    -DCMAKE_BUILD_TYPE:STRING=Release
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON
  CMAKE_CACHE_ARGS
    -DGSL_DISABLE_TESTS:BOOL=ON
    -DGSL_DISABLE_WARNINGS:BOOL=ON
    -DBUILD_SHARED_LIBS:BOOL=OFF
)

set(GSL_LIBRARIES ${CMAKE_BINARY_DIR}/stage/lib/libgsl.a ${CMAKE_BINARY_DIR}/stage/lib/libgslcblas.a)

