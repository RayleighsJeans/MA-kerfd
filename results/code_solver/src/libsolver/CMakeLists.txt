
include_directories( ../../h )

add_library( solver_slu STATIC solver_slu.cpp solver.cpp )
add_library( solver_sor STATIC solver_sor.cpp solver.cpp )
#add_library( solver_petsc STATIC solver_petsc.cpp )
