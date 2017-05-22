#pragma once

/**********************************************************************
                          Basic types definition
 **********************************************************************/ 

#include "diag_flags.h"
#include "use_flags.h"
#include <string>

typedef struct  { double r, z; }  	Field;
typedef struct  { double r, t, z; }  	Vec3d;

// mandatory part
#define PHASE_SPACE_VARS \
  double r;\
  double z;\
  double vr;\
  double vt;\
  double vz;

// optional part
#if DIAG_COLL // since i know it does not compile like this i leave it commented
#define DiagCollVars \
  double coll_r;\
  double coll_t;
#else 
#define DiagCollVars
#endif

#if USE_PARTICLE_IDS
#define PARTICLE_ID unsigned long particle_id;
#else
#define PARTICLE_ID 
#endif

#if USE_PARTICLE_ORIGIN
#define PARTICLE_ORIGIN int origin_r; int origin_z;
#else 
#define PARTICLE_ORIGIN 
#endif

#if USE_PARTICLE_WEIGHTING
#define PARTICLE_WEIGHT double w;
#else
#define PARTICLE_WEIGHT 
#endif

/// NOTE: if you want to add things here: always APPEND new types !!!!!!!!!! 
///       ( otherwise it will make old data incompatible )
enum PARTICLE_ORIGIN_TYPE_ENUM {
  UNKNOWN_ORIGIN = 0,
  INJECTION_ORIGIN,
  RECYCLING_ORIGIN,
  IONIZATION_ORIGIN,
  IONIZATION_TWOPLUS_ORIGIN
};


	//    |             $              |    
	//    | r+1,z       $              | r+1,z+1
	// ---|-------------$--------------|--------
	//		|             $              |
	//    |             $              |
	//    |    TOP_LEFT $    TOP_RIGHT |
	//    |             $              |
	//    |             $              |
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//    |             $              |
	//    |             $              |
	//    |             $              |
	//    | BOTTOM_LEFT $ BOTTOM_RIGHT |
	//    |             $              |
	//    | r,z         $              | r,z+1
	// ---|-------------$--------------|-------
	//    | r-1,z       $              | r-1,z+1
	//    |             $              |
	//
	// CELLPARTS according to 
	// fem_solver method


#if USE_FEM_SOLVER
	enum cellpart {
		TOP_LEFT,
		TOP_RIGHT,
		BOTTOM_LEFT,
		BOTTOM_RIGHT
	};
// fem solver variables for 
// diagnostic purposes
	#define FEM_SOLVER_VARS \
		double r_old;\
		double z_old;\
		cellpart oldpart;\
		cellpart newpart;
#else
	#define FEM_SOLVER_VARS
#endif

#if USE_PARTICLE_ORIGIN_TYPE
	#define PARTICLE_ORIGIN_TYPE PARTICLE_ORIGIN_TYPE_ENUM origin_type;
#else 
#define PARTICLE_ORIGIN_TYPE
#endif

#if USE_PARTICLE_BIRTH_STEP
#define PARTICLE_BIRTH_STEP int birth_step;
#else
#define PARTICLE_BIRTH_STEP
#endif


/// in the absence of inheritance the particle is 
/// composed by the preprocessor 
struct Particle {

  /// can be called from outside without assigning new ids
  /// and to reset the counter to zero 
  static unsigned long new_id ( bool set_start_id = false, long start_id = 0 ) {
#if USE_PARTICLE_IDS
    // meyers singleton
    static unsigned long id = 0;
    if ( set_start_id ) {
      //iprintf("starting new ids at %d\n", start_id);
      id = start_id;
      return id;
    }
    return id++;
#endif
  }

  /// set new id to the particle and return its id so it can be used for tracing
  unsigned long assign_new_id(){
#if USE_PARTICLE_IDS
    return particle_id = new_id();
#else
    // no particle ids used so no id
    return 0;
#endif
  }

  Particle(){
    // TODO This wont work for multiple processors
#if USE_PARTICLE_IDS
    assign_new_id();
#endif
  }
  PHASE_SPACE_VARS
	FEM_SOLVER_VARS
  PARTICLE_WEIGHT
  DiagCollVars
  PARTICLE_ID
  PARTICLE_ORIGIN
  PARTICLE_ORIGIN_TYPE
  PARTICLE_BIRTH_STEP
};

typedef struct  { double  ur, ut, uz, tr, tt, tz; int n; }	Moments;
