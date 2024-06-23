set(LIBRARY_NAME tinyXML)

file(GLOB SOURCE_FILES *.*)

add_library(${LIBRARY_NAME} STATIC ${SOURCE_FILES})
#target_link_libraries(
#        ${LIBRARY_NAME}
#        PUBLIC Geom
#        PUBLIC MeshKernel
#)