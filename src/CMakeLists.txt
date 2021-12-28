set(SOURCE
    point.cpp
    airfoil.cpp
    panel.cpp
    profiledeserializer.cpp
    main.cpp
    )
add_executable(integral ${SOURCE})
target_compile_features(integral PUBLIC cxx_std_17)
find_package(Eigen3 3.3.9 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)
target_link_libraries(integral PRIVATE Eigen3::Eigen OpenMP::OpenMP_CXX GSL::gsl)
if(MSVC)
  target_compile_options(integral PRIVATE /W4 /WX)
else()
  target_compile_options(integral PRIVATE -Wall -Wextra -Wpedantic -Werror -O3)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      target_compile_options(integral PRIVATE -fuse-ld=/usr/local/lib/mold)
  endif()
endif()
