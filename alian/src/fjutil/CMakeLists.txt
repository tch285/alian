set(MODULE_NAME fjutil)
set(MODULE ${PROJECT_NAME}_${MODULE_NAME})
set(NAME_LIB ${MODULE})

# Find Python
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
include_directories(${Python3_INCLUDE_DIRS})
link_directories(${Python3_LIBRARY_DIRS})

message(STATUS "${Green}PYTHON_LIBS: ${Python3_LIBRARIES}${ColourReset}")

message(STATUS "${Magenta}Module: ${MODULE}${ColourReset}")

include_directories(${NAME_LIB} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB SOURCES_LIB ${CMAKE_CURRENT_SOURCE_DIR}/*.cxx)
message(STATUS "SOURCES: ${SOURCES_LIB}")
string(REPLACE ".cxx" ".hh" HEADERS_LIB "${SOURCES_LIB}")

add_library(${NAME_LIB} SHARED ${SOURCES_LIB})
target_include_directories(${NAME_LIB} PUBLIC ${ROOT_INCLUDE_DIRS} ${FASTJET_DIR}/include ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(${NAME_LIB} PUBLIC ${FASTJET_LIBS} ${ROOT_LIBRARIES} ${Python3_LIBRARIES})

# Define NPY_NO_DEPRECATED_API to disable deprecated NumPy API warnings
target_compile_definitions(${NAME_LIB} PUBLIC -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)

# get the path to the root libraries using root-config --libdir
execute_process(COMMAND root-config --libdir OUTPUT_VARIABLE ROOT_LIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "${Green}FASTJET_DIR: ${FASTJET_DIR}${ColourReset}")
message(STATUS "${Green}ROOT_LIBDIR: ${ROOT_LIBDIR}${ColourReset}")

# Add the ROOT library directory to the link directories
target_link_directories(${NAME_LIB} PUBLIC ${FASTJET_DIR}/lib ${ROOT_LIBDIR})

# get the RPATH for PYTHON libraries
execute_process(COMMAND python3-config --ldflags OUTPUT_VARIABLE PYTHON_RPATH OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "${Green}PYTHON_RPATH: ${PYTHON_RPATH}${ColourReset}")

# Set the RPATH to include the ROOT and FastJet library directories
set_target_properties(${NAME_LIB} PROPERTIES
    INSTALL_RPATH "${FASTJET_DIR}/lib;${ROOT_LIBDIR};${Python_LIBDIR}"
    BUILD_WITH_INSTALL_RPATH TRUE
)

install_targets_custom()