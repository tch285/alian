set(MODULE_NAME demo)
set(MODULE ${PROJECT_NAME}_${MODULE_NAME})
set(NAME_LIB ${MODULE})

message(STATUS "${Magenta}Module: ${MODULE}${ColourReset}")

include_directories(${NAME_LIB} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB SOURCES_LIB ${CMAKE_CURRENT_SOURCE_DIR}/*.cxx)
message(STATUS "SOURCES: ${SOURCES_LIB}")
string(REPLACE ".cxx" ".hh" HEADERS_LIB "${SOURCES_LIB}")

add_library(${NAME_LIB} SHARED ${SOURCES_LIB})
target_include_directories(${NAME_LIB} PUBLIC ${ROOT_INCLUDE_DIRS} ${FASTJET_DIR}/include ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(${NAME_LIB} PUBLIC ${FASTJET_LIBS} ${ROOT_LIBRARIES})

# Get the path to the ROOT libraries using root-config --libdir
execute_process(COMMAND root-config --libdir OUTPUT_VARIABLE ROOT_LIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "${Green}FASTJET_DIR: ${FASTJET_DIR}${ColourReset}")
message(STATUS "${Green}ROOT_LIBDIR: ${ROOT_LIBDIR}${ColourReset}")

# Add the ROOT library directory to the link directories
target_link_directories(${NAME_LIB} PUBLIC ${FASTJET_DIR}/lib ${ROOT_LIBDIR})

# Set the RPATH to include the ROOT and FastJet library directories
set_target_properties(${NAME_LIB} PROPERTIES
    INSTALL_RPATH "${FASTJET_DIR}/lib;${ROOT_LIBDIR}"
    BUILD_WITH_INSTALL_RPATH TRUE
)

install_targets_custom()