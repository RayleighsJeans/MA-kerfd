
include_directories( ../../h )

SET( SOURCES_COMMON solver.cpp output_fem.cpp
										four_face_current.cpp seven_face_current.cpp ten_face_current.cpp)

ADD_LIBRARY( fem STATIC ${SOURCES_COMMON} )
