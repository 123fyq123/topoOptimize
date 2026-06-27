get_target_property(_qmake_executable Qt5::qmake IMPORTED_LOCATION)
message(STATUS "qmake executable:--------------------------------------------------------------------------------------- ${_qmake_executable}")
get_filename_component(_qt_bin_dir "${_qmake_executable}" DIRECTORY)
find_program(DEPLOYQT_EXECUTABLE NAMES windeployqt macdeployqt HINTS "${_qt_bin_dir}")

add_subdirectory(Poster)
add_subdirectory(SolverGPU)
add_subdirectory(DevGPU)
add_subdirectory(TestApp0)