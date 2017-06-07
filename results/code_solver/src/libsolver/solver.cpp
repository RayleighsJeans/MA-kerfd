#include "solver.h"
#include "var_slu.h"
#include "boundary_check.h"
#include "output.h"
#if !NTRL_ONLY
#include <superlu/slu_ddefs.h>
#endif

#include "benchmark.hpp"




/**************************************************************************
Initializing Super LU structures, filling the LHS matrix, triangulating.

Only L, U and B are declared globaly and are visible in the main()
************************************************************************/

void init_slu(std::string file_name, std::vector<std::vector<double>> Eps) 
{
#if !NTRL_ONLY

  SuperMatrix *A, *AC;    
  double   *aNR, *aNC;         
  int      *colind, *rowptr,  *rowind, *colptr;    
  int      *etree, lwork=0;  // elimination tree

  int       panel_size, permc_spec, relax;
  superlu_options_t options;
  
  double east, west, south, north, bingo;
  double e_iP1_jP1, e_iP1_jM1, e_iM1_jP1, e_iM1_jM1;
  int i, j, nnz, inz, i_row;

  double factorize_time = 0.0; 
  auto begin_time = Benchmark::start();
 
  int NR=global_grid.mesh_r_dim;
  int NZ=global_grid.mesh_z_dim;
  int NG=global_grid.mesh_dim;
  
  if ( !(aNR    = doubleMalloc(5*NG)) ) ABORT("Malloc fails for aNR[].");
  if ( !(colind = intMalloc(5*NG)) )    ABORT("Malloc fails for colind[].");
  if ( !(rowptr = intMalloc(NG+1)) )    ABORT("Malloc fails for rowptr[].");
 
  inz = 0;
  for (i = 0; i < NR; ++i) 
  {
   for (j = 0; j < NZ; ++j) 
   {
     i_row = j + i*NZ;		 
       
     if ( Eps[i][j] == -100. ) // Fixed potential on metal
     {
       rowptr[i_row] = inz;
       colind[inz] =  i_row;	  // diagonal element	 
       aNR[inz] =  1.; 
       inz++;
     } 
     else if ( Eps[i][j] == -200. ) // Fixed gradient rhs
     {
       rowptr[i_row] = inz;
       colind[inz] =  i_row - 1;	  // c - west 	 
       aNR[inz] =  -1.;    
        
       colind[inz+1] =  i_row;	     // e - bingo 	 
       aNR[inz+1] =  1.; 
       inz+=2;
     } 
      else if ( Eps[i][j] == -300. ) // Fixed gradient lhs 
     {
       rowptr[i_row] = inz;
       colind[inz] =  i_row;	     // e - bingo 	 
       aNR[inz] =  1.; 

       colind[inz+1] =  i_row + 1;	  // east   	 
       aNR[inz+1] = -1.;    
       inz+=2;
     } 
     else if ( Eps[i][j] == -400. ) // Fixed gradient upper side  
     {
       rowptr[i_row] = inz;
       colind[inz] =  i_row-NZ;	     // south	 
       aNR[inz] =  -1.; 

       colind[inz+1] =  i_row;	  // e - bingo   	 
       aNR[inz+1] = 1.;    
       inz+=2;
     } 
     else if ( i == 0)          // Axis 
     {
       e_iP1_jM1 = Eps[i+1][j-1];
       if ( e_iP1_jM1 <= 0.)   e_iP1_jM1 = Eps[i][j];
       
       e_iP1_jP1 = Eps[i+1][j+1];
       if ( e_iP1_jP1 <= 0.)   e_iP1_jP1 =  Eps[i][j];  
       
       west = e_iP1_jM1;                                                 // a,b,c,d,e    *= dz*dz
       east = e_iP1_jP1;
       north = 2.*(e_iP1_jM1 + e_iP1_jP1);
       bingo = -3.*(e_iP1_jM1 + e_iP1_jP1);    
       
       rowptr[i_row] = inz;

       colind[inz] =  i_row - 1;	  // c - west 	 
       aNR[inz] =  west; 
         
       colind[inz+1] =  i_row;	  // e - bingo 	 
       aNR[inz+1] =  bingo; 
         
       colind[inz+2] =  i_row + 1;	  // d - east 	 
       aNR[inz+2] =  east; 
            
       colind[inz+3] =  i_row + NZ;	  // b - north 	 
       aNR[inz+3] =  north;      
         		      
       inz +=4;  
     }
     else              // the rest: channell, dielectric, boundaries, etc           
     {
       e_iP1_jM1 =  Eps[i+1][j-1];
       if ( e_iP1_jM1 <= 0.)   e_iP1_jM1 = Eps[i][j];
       if ( e_iP1_jM1 ==0 )
       { 
	 if (  Eps[i+1][j]  > 0 ) e_iP1_jM1 = Eps[i+1][j];
         else                           e_iP1_jM1 = Eps[i][j-1];   
       }
                   
       e_iP1_jP1 = Eps[i+1][j+1];
       if ( e_iP1_jP1 <= 0.)   e_iP1_jP1 = Eps[i][j];  
       if ( e_iP1_jP1 ==0 ) 
       {
	 if (  Eps[i+1][j]  > 0 )   e_iP1_jP1 = Eps[i+1][j];
         else      e_iP1_jP1 = Eps[i][j+1];   
       }
      
       e_iM1_jM1 =  Eps[i-1][j-1];
       if ( e_iM1_jM1 <= 0.)   e_iM1_jM1 = Eps[i][j];
       if (e_iM1_jM1 ==0 ) 
       {
	 if ( Eps[i-1][j] > 0 ) e_iM1_jM1 = Eps[i-1][j];
         else                         e_iM1_jM1 = Eps[i][j-1];  
       }      
                       
       e_iM1_jP1 = Eps[i-1][j+1];
       if ( e_iM1_jP1 <= 0.)   e_iM1_jP1 = Eps[i][j];
       if (e_iM1_jP1 ==0 )
       { 
	 if (  Eps[i-1][j] > 0 )  e_iM1_jP1 =  Eps[i-1][j];
         else                           e_iM1_jP1 = Eps[i][j+1];   
       }
             
       north = 0.25*(e_iP1_jM1 + e_iP1_jP1)/i + 0.5*(e_iP1_jM1 + e_iP1_jP1);
       south = 0.5*(e_iM1_jM1 + e_iM1_jP1) -  0.25*(e_iM1_jM1 + e_iM1_jP1)/i;
       bingo = 0.25*(e_iM1_jM1 + e_iM1_jP1 - e_iP1_jM1 - e_iP1_jP1)/i - (e_iP1_jM1 + e_iP1_jP1 + e_iM1_jM1 + e_iM1_jP1 );
       west = 0.5*(e_iP1_jM1 + e_iM1_jM1 );
       east = 0.5*(e_iP1_jP1 + e_iM1_jP1 ); 
          
       rowptr[i_row] = inz;
         
       colind[inz] =  i_row - NZ;	  // a - south 	 
       aNR[inz] =  south; 
          
       colind[inz+1] =  i_row - 1;	  // c - west 	 
       aNR[inz+1] =  west;    
       
       colind[inz+2] =  i_row;	     // e - bingo 	 
       aNR[inz+2] =  bingo; 
       
       colind[inz+3] =  i_row + 1;	     // d - east 	 
       aNR[inz+3] =  east; 
          
       colind[inz+4] =  i_row + NZ;	 // b - north 	 
       aNR[inz+4] =  north;      
       		      
        inz +=5;  
 
     }
   }
  }
  nnz = inz;    
  rowptr[NG] = nnz;

  // Converting matrix data from CompressedRow to CompressedColumn format
  dCompRow_to_CompCol(NG, NG,  nnz, aNR, colind, rowptr, &aNC, &rowind, &colptr);
 
  // releasing unnecessary storage 
  SUPERLU_FREE (aNR);  
  SUPERLU_FREE (colind);  
  SUPERLU_FREE (rowptr);     
  
  // Allocate the L, U, B (global) Matrixes 
  L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
  L->Store = NULL;
  if ( !L ) ABORT("Malloc fails for L");
  U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );	
  U->Store = NULL;
  if ( !U ) ABORT("Malloc fails for U");
  B = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
  B->Store = NULL;
  if ( !B ) ABORT("Malloc fails for B");
	
  if ( !(perm_r = intMalloc(NG)) ) ABORT("Malloc fails for perm_r[].");
  if ( !(perm_c = intMalloc(NG)) ) ABORT("Malloc fails for perm_c[].");
    
  if ( !(etree = intMalloc(NG)) ) ABORT("Malloc fails for etree[].");
    	
  AC = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
  A = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );

	
  // Create discretization matrix A in the format expected by SuperLU. 
  dCreate_CompCol_Matrix(A, NG, NG, nnz, aNC, rowind, colptr, SLU_NC, SLU_D, SLU_GE);

  // Set the options for factorization 
  set_default_options(&options);  
  options.ColPerm =  MMD_AT_PLUS_A;       //NATURAL; //
  options.SymmetricMode = YES;
  options.DiagPivotThresh = 0.01;

  // Initialize the statistics variables. 
  StatInit(&stat_SLU);
  permc_spec = options.ColPerm;
  get_perm_c(permc_spec, A, perm_c);
  sp_preorder(&options, A, perm_c, etree, AC);
  panel_size = sp_ienv(1);
  relax = sp_ienv(2);

  //    LU factorization       
  dgstrf(&options, AC, relax, panel_size, etree, NULL, lwork, perm_c, perm_r, L, U, &stat_SLU, &info);
  assert( U->Store != NULL && "store did not change" );

  // De-allocate unnecessary storage 
  Destroy_CompCol_Permuted(AC);       
  Destroy_CompCol_Matrix(A);
  
  SUPERLU_FREE(A);
  SUPERLU_FREE(AC);
  
  SUPERLU_FREE(etree);
  
  factorize_time += Benchmark::stop(begin_time);


//create phi=0 as initial parameter
  
  // Create right-hand side matrix B. 
  dCreate_Dense_Matrix(B, NG, 1, phi, NG, SLU_DN, SLU_D, SLU_GE);          
   
#if 0 
  printf("\nCompCol matrix %s:\n", "U");
  printf("Stype %d, Dtype %d, Mtype %d\n", U->Stype,U->Dtype,U->Mtype);
  printf("nrow %d, ncol %d, nnz %d\n", U->nrow,U->ncol,((NCformat *) U->Store)->nnz);

  printf("\nSuperNode matrix %s:\n", "L");
  printf("Stype %d, Dtype %d, Mtype %d\n", L->Stype,L->Dtype,L->Mtype);
  printf("nrow %d, ncol %d, nnz %d, nsuper %d\n", L->nrow,L->ncol, ((SCformat *) L->Store)->nnz, ((SCformat *) L->Store)->nsuper);

  printf("  Factorize= %11.3e s  \n \n",   factorize_time);
#endif

#endif
}



void destroy_slu()
{
//free everything from mymem_slu.c
 /* De-allocate storage */    
#if !NTRL_ONLY

    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    
    Destroy_SuperNode_Matrix(L);
    Destroy_CompCol_Matrix(U);    
    Destroy_SuperMatrix_Store(B);
        
    SUPERLU_FREE(B);
    SUPERLU_FREE(L);
    SUPERLU_FREE(U);
    
    // release phi memory 
    // SUPERLU_FREE (phi);

    StatPrint(&stat_SLU);
    StatFree(&stat_SLU);
#endif
}
