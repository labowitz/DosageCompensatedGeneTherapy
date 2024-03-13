cmake_minimum_required(VERSION 3.20)

## make sure to use --triplet=x64-windows when doing vcpckg install
project(Analysis)
set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_TOOLCHAIN_FILE "${PROJECT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake")
set(VCPKG_TARGET_TRIPLET "x64-windows" CACHE STRING "VCPKG Target Triplet to use")
set(CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/vcpkg/packages")


set(OpenCV_DIR "${PROJECT_SOURCE_DIR}/vcpkg/packages/opencv4_x64-windows/share/opencv4/")
set(TIFF_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/vcpkg/packages/tiff_x64-windows/include")
set(TIFF_LIBRARY "${PROJECT_SOURCE_DIR}/vcpkg/packages/tiff_x64-windows/lib/tiff.lib")


#find_library(TIFF "${PROJECT_SOURCE_DIR}/vcpkg/packages/tiff_x64-windows/bin/tiff.dll")

find_package(OpenCV REQUIRED)
find_package(TIFF REQUIRED) 


set(EXECUTABLE_LIST analysisBrain2023-10-30 analysisBrain2023-11-04 analysisU2OSHCR2023-07-11_2)

FOREACH(EXECUTABLE ${EXECUTABLE_LIST})
	add_executable(${EXECUTABLE} ${EXECUTABLE}.cpp)


	target_link_libraries(${EXECUTABLE} PRIVATE
					${OpenCV_LIBS}
					${TIFF_LIBS}
					${TIFF_LIBRARY}
	)
	target_link_libraries(${EXECUTABLE} PUBLIC
					"${PROJECT_SOURCE_DIR}/LimFileStatic.lib"
					"${PROJECT_SOURCE_DIR}/Nd2ReadSdkStatic.lib"
					"${PROJECT_SOURCE_DIR}/LimFile.lib"
					"${PROJECT_SOURCE_DIR}/Nd2ReadSdk.lib"
					)

	target_include_directories(${EXECUTABLE} PUBLIC
						${PROJECT_BINARY_DIR}
						${OpenCV_INCLUDE_DIRS}
						${TIFF_INCLUDE_DIRS}
						) 

	set_property (TARGET ${EXECUTABLE}
		PROPERTY
		# Enable C++17 standard compliance
			CXX_STANDARD 17
	)
ENDFOREACH()
