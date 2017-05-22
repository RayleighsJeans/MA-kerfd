cmake_minimum_required(VERSION 2.8)

##########compile_options##########
option ( BUILD_DEBUG "build in debug mode?" off )
option ( USE_PGI "use the pgi compiler" off )
option ( RUN_CPPCHECK "runs cppcheck before the main application is build" off )

option ( NEUTRALS_ONLY "Run the code with neutrals?" off )
option ( NEUTRALS_CONSTANT "Run the code with constant neutral background?" off )
option ( TESTING "Compile Tests" off )
option ( USE_MPI "Switch MPI on or off" on ) 
option ( USE_HDF5 "Switch HDF5 on or off" on ) 
option ( USE_UNARY_COLLISIONS "does not allow more than one neutral-neutral collision for each neutral" off ) 
option ( REDUCE_NTRL_COLL "does not allow more neutral-neutral collision than numbers of neutrals in the cell" off ) 

# RF options
option ( SIE "secondary ion emission" off ) 
option ( USE_FEEDGAS_OXYGEN  "Switch between Oxygen (on) and Argon (off) as feedgas (RF_EN/RF)" off ) 

# HEMP options
option ( USE_TWOPLUS_IONS  "Double charged ions" off ) 

# diagnosis options
option ( DIAG_FLUX "Flux diagnostics at target boundary" off ) 
option ( DIAG_DIST "Diagnostics for velocity distributions" off ) 
option ( DIAG_VELZ "Diagnostics for averaged velocity per cell" off ) 
option ( DIAG_TEMPZ "Diagnostics for averaged energy and temperature per cell" off ) 
option ( DIAG_FLUID "Diagnostics for fluid lines as integration of averaged velocities" off ) 
option ( DIAG_COLL "Collisions diagnostics" off ) 
option ( DIAG_PR "Phase resolved diagnostics" off ) 
option ( DIAG_ANGULAR_CURRENT "Angular ion current and energy diagnostic" off ) 
option ( DIAG_EMISSION "emission diagnostic for one specific characteristic line" off ) 
option ( USE_PARTICLE_IDS "ids to track particles" off ) 
option ( USE_FEM_SOLVER "use fem methods for current, density and field calc" off ) 
option ( USE_PARTICLE_ORIGIN "particles can store the position at which they were generated" off ) 
option ( USE_PARTICLE_ORIGIN_TYPE "particles can store the type of creation" off ) 
option ( USE_PARTICLE_BIRTH_STEP "particles can store the time they were generated" off ) 
option ( USE_PARTICLE_TRACING "particles can be traced." off ) 

# debug options
option ( DEBUG_PARTICLE_VELOCITIES "checks particle speeds" off ) 
option ( ASAN "address sanitizer" off )




#########compile settings##########
SET( CMAKE_EXPORT_COMPILE_COMMANDS ON )
SET( CMAKE_POSITION_INDEPENDENT_CODE ON ) # needed for static libraries to link 

#make debug work
if ( ${BUILD_DEBUG} ) 
  SET( CMAKE_BUILD_TYPE "DEBUG" )
else( )
  SET( CMAKE_BUILD_TYPE "RELEASE" )
endif( )

if ( ${USE_PGI} )
    SET( CMAKE_CXX_COMPILER pgc++ )
    #SET( CMAKE_CXX_LINK_EXECUTABLE g++ )
    SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS " ")

    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc -Minfo=accel -ta=nvidia,cc35,cuda7.0 " )

    LINK_DIRECTORIES( /opt/superlu/4.3/lib )
    message( "Linker flags: ${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS}" )
endif( )

SET ( GEOM_TYPE "NONE" CACHE TYPE STRING)
SET ( SOLVER_TYPE "SLU" CACHE TYPE STRING)


if ( ${GEOM_TYPE} STREQUAL "NONE" ) 
  message( "you have to set GEOM_TYPE" )
endif( )

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGEOM_TYPE_${GEOM_TYPE}" )

if ( ${USE_PARTICLE_TRACING} AND NOT ${USE_PARTICLE_IDS} )
  message( "Particle tracing not working without Particle IDs. USE_PARTICLE_IDS=1 is set now." )
  SET( USE_PARTICLE_IDS on )
endif( )

if ( ${DIAG_EMISSION} )
  if ( ${NTRL_ONLY} )
      message( "Emission diagnostic is not working with NTRL_ONLY. DIAG_EMISSION=0 is set now." )
      SET( DIAG_EMISSION off )
  endif( )
endif( )
  
if ( ${DIAG_EMISSION} AND NOT ${DIAG_DIST} )
  message( "Emission diagnostic is not working without DIAG_DIST. DIAG_DIST=1 is set now." )
    SET( DIAG_DIST on )
endif( )

macro( add_flag_if_true FLAG FLAG_NAME ) 
  if ( ${FLAG} ) 
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D${FLAG_NAME}=1" )
  else()
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D${FLAG_NAME}=0" )
  endif()
endmacro()

# DIAGNOSTIC FLAGS
add_flag_if_true( ${DIAG_FLUX} "DIAG_FLUX" )
add_flag_if_true( ${DIAG_DIST} "DIAG_DIST" )
add_flag_if_true( ${DIAG_VELZ} "DIAG_VELZ" )
add_flag_if_true( ${DIAG_TEMPZ} "DIAG_TEMPZ" )
add_flag_if_true( ${DIAG_FLUID} "DIAG_FLUID" )
add_flag_if_true( ${DIAG_COLL} "DIAG_COLL" )
add_flag_if_true( ${DIAG_PR} "DIAG_PR" )
add_flag_if_true( ${DIAG_ANGULAR_CURRENT} "DIAG_ANGULAR_CURRENT" )
add_flag_if_true( ${DIAG_EMISSION} "DIAG_EMISSION" )
add_flag_if_true( ${DEBUG_PARTICLE_VELOCITIES} "DEBUG_PARTICLE_VELOCITIES" )
add_flag_if_true( ${NEUTRALS_ONLY} "NTRL_ONLY" )
add_flag_if_true( ${NEUTRALS_CONSTANT} "NTRL_CONST" )
add_flag_if_true( ${USE_UNARY_COLLISIONS} "USE_UNARY_COLLISIONS" )
add_flag_if_true( ${REDUCE_NTRL_COLL} "REDUCE_NTRL_COLL" )
add_flag_if_true( ${SIE} "SIE" )
add_flag_if_true( ${USE_FEEDGAS_OXYGEN} "USE_FEEDGAS_OXYGEN" )
add_flag_if_true( ${USE_TWOPLUS_IONS} "USE_TWOPLUS_IONS" )
add_flag_if_true( ${USE_HDF5} "USE_HDF5" )
add_flag_if_true( ${USE_PARTICLE_IDS}  "USE_PARTICLE_IDS" ) 
add_flag_if_true( ${USE_FEM_SOLVER}  "USE_FEM_SOLVER" ) 
add_flag_if_true( ${USE_PARTICLE_ORIGIN} "USE_PARTICLE_ORIGIN" ) 
add_flag_if_true( ${USE_PARTICLE_ORIGIN_TYPE} "USE_PARTICLE_ORIGIN_TYPE" ) 
add_flag_if_true( ${USE_PARTICLE_BIRTH_STEP}  "USE_PARTICLE_BIRTH_STEP" ) 
add_flag_if_true( ${USE_PARTICLE_TRACING}  "USE_PARTICLE_TRACING" ) 



SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fno-rtti -fopenmp -fno-exceptions -Wno-literal-suffix -Wno-unused-result" )
SET( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -g -pedantic -Wall" )
SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-exceptions -Og" )

if ( ${ASAN} ) 
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address" )
endif()

if ( ${USE_MPI} )
  SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_MPI=1" )
  SET( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DUSE_MPI=1" )
  SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DUSE_MPI=1" )
  FIND_PACKAGE( MPI REQUIRED )
else()
  SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_MPI=0" )
  SET( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DUSE_MPI=0" )
  SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DUSE_MPI=0" )
endif()



##########build files##########
include_directories( libpush )
add_subdirectory( libpush )
include_directories( libinjection )
add_subdirectory( libinjection )
include_directories( libsolver )
add_subdirectory( libsolver )
include_directories( libcollision )
add_subdirectory( libcollision )
include_directories( libfem )
add_subdirectory( libfem )

SET( SOURCES_COMMON main.cpp aver.cpp density.cpp domain_init.cpp efield.cpp emission.cpp 
    engy.cpp memory.cpp init.cpp rng64.cpp
    moms.cpp output.cpp init_diag.cpp fill_rhs.cpp 
    domain_decomposition.cpp bdd.cpp backup_and_restore.cpp restart.cpp particle_tracing.cpp
    read_input_file.cpp h5save.cpp boundary_check.cpp angular_current_diag.cpp )

SET( SOURCES_EADS "${SOURCES_COMMON}"  init_eads.cpp bfield_eads.cpp)

SET( SOURCES_HEMP "${SOURCES_COMMON}" init_hemp.cpp bfield_hemp.cpp)

SET( SOURCES_RF "${SOURCES_COMMON}"  init_rf.cpp phase_resolved.cpp bfield_rf.cpp)

INCLUDE_DIRECTORIES( ../h )

SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ../wd )  


function ( add_cppcheck_target_for_target TARGET )

    get_target_property(cppcheckSourceFiles ${TARGET} SOURCES)
    message( STATUS "SOURCES ${cppcheckSourceFiles}" )

    # get include dirs of the target
    get_target_property(cppcheckIncludeDirs ${TARGET} INCLUDE_DIRECTORIES)
    set( cppcheckIncludes )
    foreach(includeDir ${cppcheckIncludeDirs} )
	    LIST (APPEND cppcheckIncludes "-I${includeDir}")
    endforeach()

    message( STATUS "INCLUDES ${cppcheckIncludes}" )

    # get compile definitions
    
    IF(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "")
	    STRING(TOUPPER ${CMAKE_BUILD_TYPE} buildTypeUpper)
    ENDIF()
    IF ("${buildTypeUpper}" STREQUAL "")
	    SET(buildTypeUpper "DEBUG")
    ENDIF()
    

    get_target_property(cppcheckDefinesRaw ${TARGET} COMPILE_DEFINITIONS_${buildTypeUpper})
    set( cppcheckDefines )
    foreach(cppcheckDefine ${cppcheckDefinesRaw} )
	    LIST (APPEND cppcheckDefines "-D${cppcheckDefine}")
    endforeach()

    message( STATUS "DEFINITIONS ${cppcheckDefines}" )

    # TODO implement a FindCppCheck to find the executable 
    SET( CPPCHECK_EXECUTABLE cppcheck )
    message( STATUS "EXE ${CPPCHECK_EXECUTABLE}" )

    file( RELATIVE_PATH RPATH ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} )

    add_custom_command(TARGET ${TARGET}
	    PRE_BUILD
	    COMMAND ${CPPCHECK_EXECUTABLE}
	    ARGS -q --inline-suppr --enable=all
	    --suppress=incorrectStringBooleanError --suppress=unusedFunction
	    --template='${RPATH}/{file} {line}
	    \"{message}\" {id}' 
	    ${cppcheckDefines} ${cppcheckSourceFiles} 2> cpp_check_result.txt
	    )
endfunction()

function( pic_add_executable EXE SOURCES )

  # build a list of libraries that will be used 
  SET( LIBRARIES collision injection fem )

  # add the solver library
  if ( ${SOLVER_TYPE} STREQUAL "SLU" )
    SET ( LIBRARIES ${LIBRARIES} solver_slu )
    if ( NOT NEUTRALS_ONLY )
    SET ( LIBRARIES ${LIBRARIES} superlu )
    endif()
  endif()

  if ( ${SOLVER_TYPE} STREQUAL "SOR" )
    SET ( LIBRARIES ${LIBRARIES} solver_sor )
    if ( NOT NEUTRALS_ONLY )
    SET ( LIBRARIES ${LIBRARIES} superlu )
    endif()
  endif()

  if ( ${SOLVER_TYPE} STREQUAL "PETSC" )
    SET ( LIBRARIES ${LIBRARIES} solver_petsc )
    if ( NOT NEUTRALS_ONLY )
    SET ( LIBRARIES ${LIBRARIES} petsc )
    endif()
  endif()

  SET( LIBRARIES ${LIBRARIES} m blas push)
  if ( USE_HDF5 )
    SET( LIBRARIES ${LIBRARIES} hdf5 )
  endif()

  if ( ${USE_PGI} ) 
    SET( LIBRARIES ${LIBRARIES} pgftnrtl rt )
  endif()

  IF( ${USE_MPI} ) 
    IF( MPI_FOUND )
      include_directories(${MPI_INCLUDE_PATH})
      SET( LIBRARIES ${LIBRARIES} ${MPI_LIBRARIES} )
    ENDIF()
  ENDIF()

  # build executable
  add_executable( ${EXE} ${SOURCES} )
  # link to the executable
  target_link_libraries( ${EXE} ${LIBRARIES} )
  
  # build library
  add_library( pic_lib SHARED ${SOURCES} )
  # link to the library
  target_link_libraries( pic_lib ${LIBRARIES} )

  # build the converter
  add_executable( pic_converter pic_converter.cpp save.cpp read.cpp  )
  # link the pic_lib to converter
  target_link_libraries( pic_converter pic_lib )

  if( ${TESTING} ) 
    add_subdirectory( test )
  endif()

  if ( ${RUN_CPPCHECK} )
    # add cppcheck for this target
    add_cppcheck_target_for_target( ${EXE} )
  endif()

endfunction( pic_add_executable )

SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DGEOM_TYPE_${GEOM_TYPE}=1" )
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGEOM_TYPE_${GEOM_TYPE}=1" )

if ( ${GEOM_TYPE} STREQUAL "EADS" )
  pic_add_executable( eads2Dsav "${SOURCES_EADS}" )
endif()

if ( ${GEOM_TYPE} STREQUAL "HEMP" ) 
  pic_add_executable( hemp2Dsav "${SOURCES_HEMP}" )
endif()

if ( ${GEOM_TYPE} STREQUAL "RF" ) 
  pic_add_executable( rf2Dsav "${SOURCES_RF}" )
endif()


