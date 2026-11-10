cmake_minimum_required(VERSION 3.20)
project(cppquant_options VERSION 0.1.0 LANGUAGES CXX)

option(CPPQUANT_BUILD_TESTS "Build the unit tests" ON)
option(CPPQUANT_BUILD_BENCHMARKS "Build the benchmark" ON)
option(CPPQUANT_BUILD_EXAMPLES "Build the examples" ON)
option(CPPQUANT_SANITIZE "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(CPPQUANT_WERROR "Treat compiler warnings as errors" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

function(cppquant_set_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4)
        if(CPPQUANT_WERROR)
            target_compile_options(${target} PRIVATE /WX)
        endif()
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Wshadow)
        if(CPPQUANT_WERROR)
            target_compile_options(${target} PRIVATE -Werror)
        endif()
    endif()
endfunction()

if(CPPQUANT_SANITIZE)
    if(MSVC)
        message(FATAL_ERROR "CPPQUANT_SANITIZE requires Clang or GCC")
    endif()
    add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address,undefined)
endif()

add_library(cppquant_options
    src/black_scholes.cpp
    src/implied_volatility.cpp
    src/binomial.cpp
    src/monte_carlo.cpp
    src/svi.cpp
    src/surface.cpp
    src/american_mc.cpp
    src/pde.cpp
)
add_library(cppquant::options ALIAS cppquant_options)
target_include_directories(cppquant_options PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(cppquant_options PUBLIC cxx_std_20)
set_target_properties(cppquant_options PROPERTIES EXPORT_NAME options)
cppquant_set_warnings(cppquant_options)

if(CPPQUANT_BUILD_EXAMPLES)
    add_executable(option_pricer
        examples/option_pricer.cpp
    )
    target_link_libraries(option_pricer PRIVATE cppquant::options)
    cppquant_set_warnings(option_pricer)
endif()

if(CPPQUANT_BUILD_BENCHMARKS)
    add_executable(options_bench
        bench/bench_options.cpp
    )
    target_link_libraries(options_bench PRIVATE cppquant::options)
    cppquant_set_warnings(options_bench)
endif()

if(CPPQUANT_BUILD_TESTS)
    enable_testing()

    add_executable(options_tests
        tests/main.cpp
        tests/test_options.cpp
        tests/test_svi.cpp
        tests/test_surface.cpp
        tests/test_american_mc.cpp
        tests/test_pde.cpp
    )
    target_link_libraries(options_tests PRIVATE cppquant::options)
    target_include_directories(options_tests PRIVATE tests)
    cppquant_set_warnings(options_tests)

    add_test(NAME options_tests COMMAND options_tests)
    set_tests_properties(options_tests PROPERTIES TIMEOUT 180)

    if(TARGET option_pricer)
        add_test(NAME option_pricer_american
            COMMAND option_pricer --type put --style american --steps 200 --no-mc
        )
        set_tests_properties(option_pricer_american PROPERTIES TIMEOUT 60)

        add_test(NAME option_pricer_european
            COMMAND option_pricer --type call --spot 100 --strike 100 --time 1
                    --vol 0.2 --rate 0.05 --market 10.450583572186
                    --mc-paths 20000
        )
        set_tests_properties(option_pricer_european PROPERTIES TIMEOUT 60)
    endif()

    if(TARGET options_bench)
        add_test(NAME options_bench_smoke COMMAND options_bench --quick)
        set_tests_properties(options_bench_smoke PROPERTIES TIMEOUT 120)
    endif()
endif()

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS cppquant_options
    EXPORT cppquantOptionsTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT cppquantOptionsTargets
    FILE cppquantOptionsTargets.cmake
    NAMESPACE cppquant::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cppquant_options
)

configure_package_config_file(
    cmake/cppquant_optionsConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cppquant_optionsConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cppquant_options
)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/cppquant_optionsConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/cppquant_optionsConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/cppquant_optionsConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cppquant_options
)
