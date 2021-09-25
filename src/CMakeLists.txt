set(SOURCE
    point.cpp
    airfoil.cpp
    panel.cpp
    profiledeserializer.cpp
    main.cpp
   )
add_executable(integral ${SOURCE})
target_compile_features(integral PUBLIC cxx_std_17)
find_package(Eigen3 3.4 REQUIRED)
target_link_libraries(integral PRIVATE Eigen3::Eigen) #matplot)
