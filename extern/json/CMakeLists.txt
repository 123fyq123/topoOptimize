set(LIBRARY_NAME json)

file(GLOB SOURCE_FILES *.cpp)
file(GLOB INCLUDE_FILES *.hpp)

add_library(${LIBRARY_NAME} STATIC ${INCLUDE_FILES} ${SOURCE_FILES})

#target_link_libraries(
#        ${LIBRARY_NAME}
#        PUBLIC helper_opencascade
#)