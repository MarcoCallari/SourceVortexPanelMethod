set(SOURCE
    point.cpp
    airfoil.cpp
    panel.cpp
    profiledeserializer.cpp
    main.cpp
    systemsolver.hpp
    )
add_executable(integral ${SOURCE})
target_compile_features(integral PUBLIC cxx_std_17)
find_package(Eigen3 3.3.9 REQUIRED)
target_link_libraries(integral PRIVATE Eigen3::Eigen) #matplot)
if(MSVC)
  #target_compile_options(integral PRIVATE /W4 /WX)
else()
  #target_compile_options(integral PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()
