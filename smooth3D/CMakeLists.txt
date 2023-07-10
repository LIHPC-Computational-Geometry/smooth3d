#==============================================================================
# smooth3D CORE MODULE
#==============================================================================
#project(smooth3D CXX)

#add_definitions(-DDLLEXPORT)

#==============================================================================
# CODE STRUCTURE
#==============================================================================

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/inc")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src")

set(SMOOTH3D_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(SMOOTH3D_HEADER_FILES
	inc/smooth3D/smooth3d.h
        )


set(SMOOTH3D_SRC_FILES
        inc/smooth3D/smooth.h
	src/bridgeMesquite/ConditionNumber.cpp
	src/bridgeMesquite/InverseMeanRatio.cpp
	src/bridgeMesquite/Smooth3DMesh.cpp
	src/bridgeMesquite/Smooth3DMesh.h
	src/genuine/laplace.c
	src/genuine/Tipton.c
	src/getme/GetMe2D.cpp
	src/getme/GetMe3D.cpp
	src/math/BFGS.h
	src/math/condition.h
	src/math/ConjugateGradientPr.h
	src/math/LagrangeMapping.h
	src/math/LineSearch.h
	src/math/MeanRatio.h
	src/math/Optimize.h
	src/math/orthogonal.h
	src/math/Real3.h
	src/opt2D/IObjective3D.h
	src/opt2D/ObjectiveCondition2D.cpp
	src/opt2D/ObjectiveCondition2D.h
	src/opt2D/ObjectiveLag2D.cpp
	src/opt2D/ObjectiveLag2D.h
	src/opt2D/ObjectiveOrthogonal2D.cpp
	src/opt2D/ObjectiveOrthogonal2D.h
	src/opt2D/s3_conditionnumber2D.cpp
	src/opt2D/s3_orthogonal2D.cpp
	src/tiptonJun/jun2d.cpp
	src/tiptonJun/Jun.cpp
	src/s3mem.c
	src/s3_3dgeom.c
        )


#==============================================================================
# DEPENDENCIES
#==============================================================================
find_package(MachineTypes REQUIRED)

find_package(Lima REQUIRED)

find_package(gmds REQUIRED)

include ("${CMAKE_CURRENT_SOURCE_DIR}/FindMesquite.cmake")

#==============================================================================
# TEST MODULE
#==============================================================================

#set(GMDS_TEST ON CACHE BOOL "GMDS Kernel test module")
#if (GMDS_TEST)
#  add_subdirectory(test)
#endif ()

#==============================================================================
# LIBRARY BUILDING
#==============================================================================

add_library(smooth3D ${SMOOTH3D_HEADER_FILES} ${SMOOTH3D_SRC_FILES})

#set_target_properties(gmds PROPERTIES PUBLIC_HEADER "${GMDS_HEADER_FILES}")


#target_compile_features(gmds PUBLIC cxx_std_11)


target_link_libraries(smooth3D PRIVATE gmds)
target_link_libraries(smooth3D PRIVATE Lima::Lima Mesquite::Mesquite)
target_link_libraries(smooth3D PUBLIC MachineTypes::MachineTypes)

#add_custom_command(
#   TARGET gmds
#   POST_BUILD
#   COMMAND ${CMAKE_COMMAND} ARGS -E copy_directory  $<TARGET_FILE:gmds> ${CMAKE_BINARY_DIR}/gen/
#   COMMENT "Copy gmds lib"
#   VERBATIM
#   )

#target_include_directories(gmds INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>/generated>)
target_include_directories(smooth3D PUBLIC "$<INSTALL_INTERFACE:include>")

install(TARGETS smooth3D
        EXPORT smooth3D DESTINATION
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include
        )

install(FILES
        ${SMOOTH3D_HEADER_FILES} DESTINATION include/smooth3D)


set(INCLUDE_INSTALL_DIR include/ )
set(LIB_INSTALL_DIR lib/ )
set(SYSCONFIG_INSTALL_DIR etc/gmds/ )

include(CMakePackageConfigHelpers)
configure_package_config_file(
        smooth3DConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/smooth3DConfig.cmake
        INSTALL_DESTINATION ${LIB_INSTALL_DIR}/smooth3D/cmake
        PATH_VARS INCLUDE_INSTALL_DIR SYSCONFIG_INSTALL_DIR
)

write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/smooth3DConfigVersion.cmake
        VERSION "${smooth3D_VERSION}"
        COMPATIBILITY SameMajorVersion
)

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/smooth3DConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/smooth3DConfigVersion.cmake
        DESTINATION ${LIB_INSTALL_DIR}/smooth3D/cmake
        )


install(EXPORT smooth3D DESTINATION lib/smooth3D)
