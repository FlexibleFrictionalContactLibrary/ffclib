# The solvers of the benchmark, each behind CFC_WITH_<NAME>. Clarabel, QOCO, ConicXX and SCS are
# fetched (a local checkout: -DFETCHCONTENT_SOURCE_DIR_<NAME>=<path>), cuOpt is found in an installation.
# Versions as in CardilloCxx (cmake/ConicBackends.cmake).
include(FetchContent)

if(CFC_WITH_CLARABEL)
    # Needs Rust/cargo. Clarabel.cpp passes --release to cargo only for CMAKE_BUILD_TYPE=Release; the
    # other optimized build types are presented as Release so its core is not built in cargo's debug profile.
    FetchContent_Declare(clarabel GIT_REPOSITORY https://github.com/oxfordcontrol/Clarabel.cpp.git GIT_TAG v0.11.1)
    set(_cfcbench_build_type "${CMAKE_BUILD_TYPE}")
    if(CMAKE_BUILD_TYPE MATCHES "^(RelWithDebInfo|MinSizeRel)$")
        set(CMAKE_BUILD_TYPE Release)
    endif()
    FetchContent_MakeAvailable(clarabel)
    set(CMAKE_BUILD_TYPE "${_cfcbench_build_type}")
endif()

if(CFC_WITH_QOCO)
    FetchContent_Declare(qoco GIT_REPOSITORY https://github.com/qoco-org/qoco.git GIT_TAG v0.3.2)
    FetchContent_MakeAvailable(qoco)
endif()

if(CFC_WITH_CONICXX)
    set(CONICXX_BUILD_TESTS OFF)
    set(CONICXX_BUILD_EXAMPLES OFF)
    set(CONICXX_BUILD_BENCHMARKS OFF)
    set(CONICXX_INSTALL OFF)
    FetchContent_Declare(conicxx GIT_REPOSITORY https://github.com/ConicXX/ConicXX.git GIT_TAG main)
    FetchContent_MakeAvailable(conicxx)
endif()

# SCS: only its direct linear-system backend (scsdir: AMD + QDLDL). SCS vendors AMD and QDLDL under their
# original symbol names, which QOCO and ConicXX bring into the same process; scsdir is therefore a shared
# library that exports only the public scs_* API (scs_exports.map).
if(CFC_WITH_SCS)
    find_package(LAPACK QUIET)  # Anderson acceleration is compiled only with LAPACK
    set(USE_LAPACK ${LAPACK_FOUND})
    set(BUILD_SHARED_LIBS ON)
    set(DLONG OFF)
    set(SFLOAT OFF)
    set(BUILD_TESTING OFF)  # its test programs need the symbols the version script hides
    set(_cfcbench_omp "${OpenMP_C_FOUND}")
    set(OpenMP_C_FOUND FALSE)  # serial SCS: the benchmark runs every CPU solver on one thread
    FetchContent_Declare(scs GIT_REPOSITORY https://github.com/cvxgrp/scs.git GIT_TAG 3.3.1 GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(scs)
    set(OpenMP_C_FOUND "${_cfcbench_omp}")
    unset(BUILD_SHARED_LIBS)
    unset(USE_LAPACK)
    unset(DLONG)
    unset(SFLOAT)
    unset(BUILD_TESTING)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_options(scsdir PRIVATE "LINKER:--version-script=${CMAKE_CURRENT_LIST_DIR}/scs_exports.map")
    endif()
endif()

# cuOpt: the C API of its mathematical-optimization component; its header needs no CUDA toolkit.
if(CFC_WITH_CUOPT)
    find_package(CUOPT 26.10.1 REQUIRED)
    message(STATUS "cfc_solve: cuOpt ${CUOPT_VERSION} (${CUOPT_LIBRARY})")
endif()

# The solvers' headers are not ours to fix.
foreach(_target libclarabel_c_shared conicxx scsdir qoco)
    if(TARGET ${_target})
        get_target_property(_dirs ${_target} INTERFACE_INCLUDE_DIRECTORIES)
        if(_dirs)
            set_target_properties(${_target} PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_dirs}")
        endif()
    endif()
endforeach()
