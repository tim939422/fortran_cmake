cmake_minimum_required(VERSION 3.21)
project(FortranFFTW LANGUAGES Fortran)

# Find thirdparty packages
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

message(STATUS, "SET UP FFT")
find_package(FFTW)
message(STATUS "FFTW_FOUND     : ${FFTW_FOUND}")
message(STATUS "FFTW_LIBRARIES : ${FFTW_LIBRARIES}")
message(STATUS "FFTW_INCLUDE   : ${FFTW_INCLUDE_DIRS}")
# Make sure we use the double precision only
if(FFTW_FOUND)
    if(FFTW_DOUBLE_LIB_FOUND)
        message(STATUS "FFTW_DOUBLE     : ${FFTW_DOUBLE_LIB}")
    endif(FFTW_DOUBLE_LIB_FOUND)    
endif(FFTW_FOUND)

add_executable(main
    main.f90
    m_tests.f90
    m_fftw3.f90
)
target_link_libraries(main PRIVATE ${FFTW_DOUBLE_LIB})
