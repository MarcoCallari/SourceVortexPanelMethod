set(SOURCE
    point.cpp
    airfoil.cpp
    panel.cpp
    main.cpp
   )
add_executable(integral ${SOURCE})
find_package(Eigen3 3.4 REQUIRED)
target_link_libraries(integral PRIVATE Eigen3::Eigen matplot)
