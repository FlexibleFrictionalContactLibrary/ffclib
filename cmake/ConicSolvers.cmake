# The conic solvers, each behind FFCLIB_WITH_<NAME>: its source joins the library, which then defines
# FFCLIB_WITH_<NAME> for its users. Clarabel, SCS, QOCO and ConicXX are fetched (a local checkout:
# -DFETCHCONTENT_SOURCE_DIR_<NAME>=<path>), cuOpt is found in an installation (FindCUOPT.cmake, CUOPT_ROOT).
include(FetchContent)

if(FFCLIB_WITH_CLARABEL)
    # Needs Rust/cargo. Clarabel.cpp builds its Rust core in cargo's release profile only for CMAKE_BUILD_TYPE
    # Release; the other optimized build types are presented as Release (cargo's debug profile is 15-60x slower).
    FetchContent_Declare(clarabel GIT_REPOSITORY https://github.com/oxfordcontrol/Clarabel.cpp.git GIT_TAG v0.11.1)
    set(_ffclib_build_type "${CMAKE_BUILD_TYPE}")
    if(CMAKE_BUILD_TYPE MATCHES "^(RelWithDebInfo|MinSizeRel)$")
        set(CMAKE_BUILD_TYPE Release)
    endif()
    FetchContent_MakeAvailable(clarabel)
    set(CMAKE_BUILD_TYPE "${_ffclib_build_type}")
    target_sources(ffclib PRIVATE src/ffclib/solver/conic/clarabel.cpp)
    target_link_libraries(ffclib PRIVATE libclarabel_c_shared)
    target_compile_definitions(ffclib PUBLIC FFCLIB_WITH_CLARABEL)
endif()

# SCS with its direct linear-system solver only (scsdir: AMD + QDLDL), serial (every CPU solver runs on one
# thread). SCS vendors AMD and QDLDL under their original symbol names, which other solvers bring into the same
# process too (QOCO, ConicXX); scsdir is therefore a shared library that exports only the scs_* API.
if(FFCLIB_WITH_SCS)
    find_package(LAPACK QUIET)  # SCS compiles its Anderson acceleration only with LAPACK
    set(USE_LAPACK ${LAPACK_FOUND})
    set(BUILD_SHARED_LIBS ON)
    set(DLONG OFF)
    set(SFLOAT OFF)
    set(BUILD_TESTING OFF)  # SCS's test programs need the symbols that the version script hides
    set(_ffclib_openmp "${OpenMP_C_FOUND}")
    set(OpenMP_C_FOUND FALSE)
    FetchContent_Declare(scs GIT_REPOSITORY https://github.com/cvxgrp/scs.git GIT_TAG 3.3.1 GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(scs)
    set(OpenMP_C_FOUND "${_ffclib_openmp}")
    unset(USE_LAPACK)
    unset(BUILD_SHARED_LIBS)
    unset(DLONG)
    unset(SFLOAT)
    unset(BUILD_TESTING)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_options(scsdir PRIVATE "LINKER:--version-script=${CMAKE_CURRENT_LIST_DIR}/scs_exports.map")
    endif()
    target_sources(ffclib PRIVATE src/ffclib/solver/conic/scs.cpp)
    target_link_libraries(ffclib PRIVATE scsdir)
    target_compile_definitions(ffclib PUBLIC FFCLIB_WITH_SCS)
endif()

# QOCO, statically linked (its AMD and QDLDL do not clash: SCS hides its copies).
if(FFCLIB_WITH_QOCO)
    FetchContent_Declare(qoco GIT_REPOSITORY https://github.com/qoco-org/qoco.git GIT_TAG v0.3.2)
    FetchContent_MakeAvailable(qoco)
    target_sources(ffclib PRIVATE src/ffclib/solver/conic/qoco.cpp)
    target_link_libraries(ffclib PRIVATE qocostatic)
    target_compile_definitions(ffclib PUBLIC FFCLIB_WITH_QOCO)
endif()

# ConicXX (C++, Eigen), the commit that CardilloCxx uses.
if(FFCLIB_WITH_CONICXX)
    set(CONICXX_BUILD_TESTS OFF)
    set(CONICXX_BUILD_EXAMPLES OFF)
    set(CONICXX_BUILD_BENCHMARKS OFF)
    set(CONICXX_INSTALL OFF)
    FetchContent_Declare(conicxx GIT_REPOSITORY https://github.com/ConicXX/ConicXX.git GIT_TAG 05b3b42a3a75eae79b90d034f0f02020a8ffdddb)
    FetchContent_MakeAvailable(conicxx)
    target_sources(ffclib PRIVATE src/ffclib/solver/conic/conicxx.cpp)
    target_link_libraries(ffclib PRIVATE conicxx::conicxx)
    target_compile_definitions(ffclib PUBLIC FFCLIB_WITH_CONICXX)
endif()

# cuOpt: the C API of its mathematical-optimization component; its header needs no CUDA toolkit.
if(FFCLIB_WITH_CUOPT)
    find_package(CUOPT 26.10.1 REQUIRED)
    message(STATUS "ffclib: cuOpt ${CUOPT_VERSION} (${CUOPT_LIBRARY})")
    target_sources(ffclib PRIVATE src/ffclib/solver/conic/cuopt.cpp)
    target_link_libraries(ffclib PRIVATE CUOPT::CUOPT)
    target_compile_definitions(ffclib PUBLIC FFCLIB_WITH_CUOPT)
endif()

# The solvers' headers are not ours to fix.
foreach(_target libclarabel_c_shared scsdir qocostatic conicxx)
    if(TARGET ${_target})
        get_target_property(_dirs ${_target} INTERFACE_INCLUDE_DIRECTORIES)
        if(_dirs)
            set_target_properties(${_target} PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_dirs}")
        endif()
    endif()
endforeach()
