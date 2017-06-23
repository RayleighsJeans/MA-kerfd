/***********************************************************
2d3v Axially Symmetric PIC-MCC ( ASPIC )
24.01.10
************************************************************/

#if USE_FEM_SOLVER
    if (nstep == 1) {
      // very important for fem
      // needed to relax efield before first step
      // to null everything
      reset_fem ( ); store_old_field ( );
      // grab efield from matrix solver
      efield_matrix2fem ( E_grid );
      // dump to oldfield used in first fem solve
      store_old_field ( ); reset_fem ( );
    }

    if (nstep >= 2) reset_fem ( );
    area_density ( );
    cellfacecurrent ( );
    efield_fem ( );
    output_fem ( E_grid );
    store_old_field ( );
#endif
  
