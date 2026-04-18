cmake_minimum_required(VERSION 3.16)

project(ASEN_6044_Final_Project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

find_package(Eigen3 REQUIRED)
find_package(yaml-cpp REQUIRED)

file(
	GLOB_RECURSE PROJECT_SOURCES CONFIGURE_DEPENDS
	"${SRC_DIR}/*.cpp"
	"${SRC_DIR}/*.cc"
	"${SRC_DIR}/*.cxx"
)

if(NOT PROJECT_SOURCES)
	message(WARNING "No C++ source files found in ${SRC_DIR}. Add at least one .cpp file.")
else()
	add_executable(${PROJECT_NAME} ${PROJECT_SOURCES})

	target_include_directories(${PROJECT_NAME}
		PRIVATE
			${INCLUDE_DIR}
	)

    target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)
	target_link_libraries(${PROJECT_NAME} yaml-cpp::yaml-cpp)

	if(MSVC)
		target_compile_options(${PROJECT_NAME} PRIVATE /W4)
	else()
		target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic)
	endif()
endif()