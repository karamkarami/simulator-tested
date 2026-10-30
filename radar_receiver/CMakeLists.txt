# Radar receiver: processes the DDC output that the MATLAB simulator streams
# over UDP (DBF, matched filter, MTI, Doppler, 2-D CFAR, plot extraction),
# with an optional Qt PPI and an optional serial output of the detections.
# CMake 3.21: C_STANDARD 11 gives MSVC /std:c11 (needed for 'restrict').
cmake_minimum_required(VERSION 3.21)
project(radar_receiver LANGUAGES C)

option(RADAR_WITH_QT     "Build the Qt PPI display (Qt 6 or Qt 5 Widgets)" ON)
option(RADAR_OPENMP      "Use OpenMP threads"                              ON)
option(RADAR_NATIVE      "Optimise for the CPU of the build machine"       ON)
option(RADAR_BUILD_TESTS "Build the unit tests (ctest)"                    ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Default Qt location on Windows; appended, so -DCMAKE_PREFIX_PATH=<Qt> still wins
list(APPEND CMAKE_PREFIX_PATH "C:/Qt/6.12.0/mingw_64")

# Compiler flags. No -ffast-math: the elevation of a plot is NaN without
# beamforming and the processing relies on IEEE behaviour.
include(CheckCCompilerFlag)
if(MSVC)
    add_compile_options(/W3 /fp:precise)
    string(APPEND CMAKE_C_FLAGS_RELEASE " /O2")
    if(RADAR_NATIVE)
        # MSVC has no "native": AVX2 (every Intel Core CPU since 2013)
        add_compile_options(/arch:AVX2)
    endif()
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
else()
    add_compile_options(-Wall -Wextra -fno-math-errno)
    set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
    if(RADAR_NATIVE)
        check_c_compiler_flag(-march=native RADAR_HAS_MARCH_NATIVE)
        if(RADAR_HAS_MARCH_NATIVE)
            add_compile_options(-march=native)
        else()
            message(STATUS "The compiler does not accept -march=native: generic optimisation")
        endif()
    endif()
endif()

# ---------------------------------------------------------------------------
# Processing library
# ---------------------------------------------------------------------------
add_library(radar_core STATIC
    src/config.c
    src/fft.c
    src/extractor.c
    src/tracker.c
    src/processor.c
    src/udp_input.c
    src/serial_out.c
    src/platform.c
    src/input_thread.c
)
target_include_directories(radar_core PUBLIC src)

if(RADAR_OPENMP)
    find_package(OpenMP COMPONENTS C)
    if(OpenMP_C_FOUND)
        target_link_libraries(radar_core PUBLIC OpenMP::OpenMP_C)
        # OpenMP::OpenMP_C only adds -fopenmp at link time when the C linker
        # is used ($<LINK_LANGUAGE:C>). Programs that also link the Qt PPI use
        # the C++ linker, so pass the flag for that case too; otherwise GCC /
        # MinGW gives "undefined reference to GOMP_parallel, omp_get_num_threads".
        if(NOT MSVC AND OpenMP_C_FLAGS)
            target_link_options(radar_core INTERFACE
                "$<$<LINK_LANGUAGE:CXX>:${OpenMP_C_FLAGS}>")
        endif()
    else()
        message(WARNING "OpenMP not found: the receiver runs single-threaded")
    endif()
endif()
if(NOT OpenMP_C_FOUND AND NOT MSVC)
    target_compile_options(radar_core PRIVATE -Wno-unknown-pragmas)   # the omp pragmas are ignored
endif()
find_package(Threads REQUIRED)
target_link_libraries(radar_core PUBLIC Threads::Threads)
if(WIN32)
    target_link_libraries(radar_core PUBLIC ws2_32)
else()
    target_link_libraries(radar_core PUBLIC m)
endif()

# ---------------------------------------------------------------------------
# Optional Qt PPI
# ---------------------------------------------------------------------------
set(RADAR_QT_FOUND OFF)
if(RADAR_WITH_QT)
    # Qt's CMake packages need C++ enabled before they are searched
    include(CheckLanguage)
    check_language(CXX)
    if(CMAKE_CXX_COMPILER)
        enable_language(CXX)
        set(CMAKE_CXX_STANDARD 17)
        set(CMAKE_CXX_STANDARD_REQUIRED ON)
        find_package(QT NAMES Qt6 Qt5 COMPONENTS Widgets QUIET)
    else()
        message(WARNING "No C++ compiler: building without the PPI display")
    endif()
    if(QT_FOUND AND QT_VERSION VERSION_LESS 5.12)
        message(WARNING "Qt ${QT_VERSION} is too old (5.12 or newer is needed): building without the PPI display")
        set(QT_FOUND OFF)
    endif()
    if(QT_FOUND)
        find_package(Qt${QT_VERSION_MAJOR} COMPONENTS Widgets REQUIRED)
        find_package(Threads REQUIRED)
        add_library(radar_ppi STATIC src/ppi_qt.cpp)
        set_target_properties(radar_ppi PROPERTIES AUTOMOC ON)
        target_include_directories(radar_ppi PUBLIC src)
        target_link_libraries(radar_ppi PUBLIC Qt${QT_VERSION_MAJOR}::Widgets Threads::Threads)
        set(RADAR_QT_FOUND ON)
        message(STATUS "PPI display: Qt ${QT_VERSION}")
    else()
        message(WARNING "Qt Widgets not found: building without the PPI display")
    endif()
endif()

# ---------------------------------------------------------------------------
# Receiver program
# ---------------------------------------------------------------------------
add_executable(radar_receiver src/main.c)
target_link_libraries(radar_receiver PRIVATE radar_core)
if(RADAR_QT_FOUND)
    target_compile_definitions(radar_receiver PRIVATE RADAR_WITH_QT)
    target_link_libraries(radar_receiver PRIVATE radar_ppi)
    # The C main program links with the C++ PPI library.
    set_target_properties(radar_receiver PROPERTIES LINKER_LANGUAGE CXX)
endif()

# Reading the DDC files of the simulator (ddc_<n>.mat). zlib (optional) is
# needed only for compressed MAT-files (save -v7); -v6 files need nothing.
add_library(radar_ddc_files STATIC tools/mat_file.c tools/ddc_files.c)
target_include_directories(radar_ddc_files PUBLIC tools)
find_package(ZLIB QUIET)
if(ZLIB_FOUND)
    target_compile_definitions(radar_ddc_files PRIVATE RADAR_HAVE_ZLIB)
    target_link_libraries(radar_ddc_files PRIVATE ZLIB::ZLIB)
    message(STATUS "DDC files: -v6 and compressed -v7 MAT-files (zlib)")
else()
    message(STATUS "zlib not found: the DDC files must be saved with -v6 (cfg.io.matFormat)")
endif()
if(NOT WIN32)
    target_link_libraries(radar_ddc_files PUBLIC m)
endif()

# Sends the DDC files (ddc_<n>.mat) or a recorded stream (--record) over UDP
# with the real PRI timing
add_executable(radar_stream_send tools/stream_send.c)
target_link_libraries(radar_stream_send PRIVATE radar_core radar_ddc_files)

# Windows: copy the Qt DLLs and plugins next to a program that uses the PPI,
# so it starts without Qt's bin folder on PATH.
function(radar_deploy_qt target)
    if(WIN32 AND RADAR_QT_FOUND)
        get_target_property(qmake Qt${QT_VERSION_MAJOR}::qmake IMPORTED_LOCATION)
        get_filename_component(qtBin "${qmake}" DIRECTORY)
        find_program(WINDEPLOYQT windeployqt HINTS "${qtBin}")
        if(WINDEPLOYQT)
            add_custom_command(TARGET ${target} POST_BUILD
                COMMAND "${WINDEPLOYQT}" --no-translations --no-compiler-runtime "$<TARGET_FILE:${target}>"
                COMMENT "windeployqt ${target}")
        else()
            message(WARNING "windeployqt not found: add Qt's bin folder to PATH to run ${target}")
        endif()
    endif()
endfunction()
radar_deploy_qt(radar_receiver)

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(RADAR_BUILD_TESTS)
    enable_testing()

    add_executable(test_core tests/test_core.c)
    target_link_libraries(test_core PRIVATE radar_core)
    add_test(NAME core COMMAND test_core ${CMAKE_CURRENT_SOURCE_DIR}/config/radar_config.bin)
    # phase codes with mismatch filters (export_receiver_config.m, see README)
    add_test(NAME core_phasecode COMMAND test_core ${CMAKE_CURRENT_SOURCE_DIR}/config/radar_config_phasecode.bin)

    # the CFAR against direct sums (white-box: the test includes processor.c)
    add_executable(test_cfar tests/test_cfar.c src/config.c src/fft.c src/extractor.c src/tracker.c src/platform.c)
    target_include_directories(test_cfar PRIVATE src)
    target_link_libraries(test_cfar PRIVATE Threads::Threads)
    if(OpenMP_C_FOUND)
        target_link_libraries(test_cfar PRIVATE OpenMP::OpenMP_C)
    endif()
    if(NOT WIN32)
        target_link_libraries(test_cfar PRIVATE m)
    endif()
    add_test(NAME cfar COMMAND test_cfar ${CMAKE_CURRENT_SOURCE_DIR}/config/radar_config.bin)

    add_executable(test_tracker tests/test_tracker.c)
    target_link_libraries(test_tracker PRIVATE radar_core)
    add_test(NAME tracker COMMAND test_tracker ${CMAKE_CURRENT_SOURCE_DIR}/config/radar_config.bin)

    add_executable(test_ddc_files tests/test_ddc_files.c)
    target_link_libraries(test_ddc_files PRIVATE radar_core radar_ddc_files)
    if(ZLIB_FOUND)
        target_compile_definitions(test_ddc_files PRIVATE RADAR_HAVE_ZLIB)   # writes compressed files
        target_link_libraries(test_ddc_files PRIVATE ZLIB::ZLIB)
    endif()
    add_test(NAME ddc_files COMMAND test_ddc_files)

    add_executable(test_udp_input tests/test_udp_input.c)
    target_link_libraries(test_udp_input PRIVATE radar_core)
    add_test(NAME udp_input COMMAND test_udp_input)

    if(NOT WIN32)
        add_executable(test_serial tests/test_serial.c)
        target_link_libraries(test_serial PRIVATE radar_core)
        find_library(UTIL_LIBRARY util)
        if(UTIL_LIBRARY)
            target_link_libraries(test_serial PRIVATE ${UTIL_LIBRARY})
        endif()
        add_test(NAME serial COMMAND test_serial)
    endif()

    if(RADAR_QT_FOUND)
        add_executable(test_ppi tests/test_ppi.c)
        target_link_libraries(test_ppi PRIVATE radar_ppi)
        set_target_properties(test_ppi PROPERTIES LINKER_LANGUAGE CXX)
        add_test(NAME ppi COMMAND test_ppi ${CMAKE_CURRENT_BINARY_DIR}/ppi_test.png)
        set_tests_properties(ppi PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen")

        # interaction (mouse, wheel, keys): the test includes ppi_qt.cpp
        add_executable(test_ppi_interactive tests/test_ppi_interactive.cpp)
        target_include_directories(test_ppi_interactive PRIVATE src)
        target_link_libraries(test_ppi_interactive PRIVATE Qt${QT_VERSION_MAJOR}::Widgets Threads::Threads)
        add_test(NAME ppi_interactive COMMAND test_ppi_interactive ${CMAKE_CURRENT_BINARY_DIR})
        set_tests_properties(ppi_interactive PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen")
        radar_deploy_qt(test_ppi_interactive)
        radar_deploy_qt(test_ppi)
    endif()
endif()
