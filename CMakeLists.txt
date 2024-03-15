cmake_minimum_required(VERSION 3.20)

project(Analysis)
set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")

find_package(OpenCV REQUIRED)
find_package(TIFF REQUIRED)
# find_package(ZLIB REQUIRED) 

set(EXECUTABLE_LIST analysisBrain2023-10-30 analysisBrain2023-11-04 analysisU2OSHCR2023-07-11_2)

FOREACH(EXECUTABLE ${EXECUTABLE_LIST})
	add_executable(${EXECUTABLE} ${EXECUTABLE}.cpp)


	target_link_libraries(${EXECUTABLE} PRIVATE
					${OpenCV_LIBS}
					${TIFF_LIBS}
					#${ZLIB_LIBS}
					#${TIFF_LIBRARY}
	)
	target_link_libraries(${EXECUTABLE} PRIVATE
					"${PROJECT_SOURCE_DIR}/limfile-shared.lib"
					"${PROJECT_SOURCE_DIR}/nd2readsdk-shared.lib"
					)

	target_include_directories(${EXECUTABLE} PUBLIC
						${PROJECT_BINARY_DIR}
						${OpenCV_INCLUDE_DIRS}
						${TIFF_INCLUDE_DIRS}
						#${ZLIB_INCLUDE_DIRS}
						) 

	set_property (TARGET ${EXECUTABLE}
		PROPERTY
		# Enable C++17 standard compliance
			CXX_STANDARD 17
	)
ENDFOREACH()
