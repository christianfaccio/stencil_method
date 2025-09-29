#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <float.h>
#include <math.h>


#define NORTH 0
#define SOUTH 1
#define EAST  2
#define WEST  3

#define SEND 0
#define RECV 1

#define OLD 0
#define NEW 1

#define _x_ 0
#define _y_ 1

#define uint unsigned int 

// ============================================================

int initialize (    int         argc,			        // argc from cmd line
		            char	    **argv,			        // argv from cmd line
		            int	        *S,			            // two-uint array defining x,y dims of grid
		            int	        *periodic,		        // boolean flag for boundary conditions - if 1, 
                                                        // the plate behaves as infinite with periodic 
                                                        // boundaries; if 0, fixed boundaries
		            int 	    *Niterations,		    // n. iterations
		            int   	    *Nsources,		        // n. heat sources randomly placed on the grid 
                                                        // (default = 1)
		            int  	    **Sources,              // 2D array storing heat source positions;
                                                        // sources[2*i] and Sources[2*i+1] are x,y 
                                                        // coordinates of source i
		            double      *energy_per_source,	    // amount of energy each heat source injects 
                                                        // per injection event
		            double      **planes,		        // two 2D arrays for double-buffering;
                                                        // one holds current heat values, the other 
                                                        // stores updated values (swapped each iteration)
                    int   	    *output_energy_at_steps,// debug flag - if 1, prints energy statistics 
                                                        // and dumps grid data at each time step
                    int   	    *injection_frequency    // how often to inject energy
		        );  

int update_plane (const int     periodic, 
                  const int     size[2],
			      const double  *old,
                  double        *new
                );   

int dump (          const double  *data, 
                    const uint    size[2], 
                    const char    *filename, 
                    double        *min, 
                    double        *max
        );


// ============================================================

static inline int memory_allocate ( const int       size[2],
		                            double	        **planes_ptr
                                )
/*
 * allocate the memory for the planes
 * we need 2 planes: the first contains the
 * current data, the second the updated data
 *
 * in the integration loop then the roles are
 * swapped at every iteration
 *
 */
{
  if (planes_ptr == NULL )  
    return 1;

  /*
  use +2 for ghost cells:
  Ghost cells:    G G G G G G
                  G x x x x G  <- Real grid points (x)
                  G x x x x G     surrounded by ghost cells (G)
                  G x x x x G
                  G x x x x G
                  G G G G G G
  */ 
  unsigned int bytes = (size[_x_]+2)*(size[_y_]+2);

  /*
  planes_ptr ──→ [planes_ptr[OLD]] ──→ [allocated memory for OLD plane]
                 [planes_ptr[NEW]] ──→ [allocated memory for NEW plane]
  */
  planes_ptr[OLD] = (double*)calloc( 2*bytes, sizeof(double));  // Allocate one big block using calloc 
                                                                // such that values are initialized to 0
  planes_ptr[NEW] = planes_ptr[OLD] + bytes; // points to second half
      
  return 0;
}

static inline int initialize_sources(   uint        size[2],
			                            int    	    Nsources,
			                            int  	    **Sources
                                    )
/*
 * randomly spread heat suources
 * NEED ENHANCEMENT FOR PARALLELISM (also lrand48)
 */
{
  // Sources ──→ [*Sources] ──→ [sources array memory]
  *Sources = (int*)malloc( Nsources * 2 *sizeof(uint) );
  for ( int s = 0; s < Nsources; s++ )
    {
      (*Sources)[s*2] = 1+ lrand48() % size[_x_]; // coordinate in range [1, size]
      (*Sources)[s*2+1] = 1+ lrand48() % size[_y_]; // coordinate in range [1, size]
    }

  return 0;
}

static inline int inject_energy (   const int       periodic,
                                    const int       Nsources,
			                        const int       *Sources,
			                        const double    energy,
			                        const int       mysize[2],
                                    double          *plane 
                                )
/*
Periodic function to inject energy into the grid at the positions specified in Sources array.
If periodic is true, it also updates the ghost cells accordingly.

Grid with periodic boundaries and source at (1,3):

[G][G][G][G][G]
[G][S][x][x][G] <- Source at (1,3) injects here
[G][x][x][x][G]
[G][x][x][x][S] <- AND also here at wrapped position
[G][G][G][G][G]
*/
{
   #define IDX( i, j ) ( (j)*(mysize[_x_]+2) + (i) )    // replaced at compile time
                                                        // IDX is a macro that converts 
                                                        // 2D coordinates (i, j) into a
                                                        // 1D array index for row-major
                                                        // storage.
    for (int s = 0; s < Nsources; s++) {
        
        int x = Sources[2*s];
        int y = Sources[2*s+1];
        plane[IDX(x, y)] += energy;
        
        // handle periodic boundaries
        // - inner cells -> [1, mysize[_x_]] x [1, mysize[_y_]]
        // - ghost cells -> [0, mysize[_x_]+1] x [0, mysize[_y_]+1]
        if ( periodic )
            {
                if ( x == 1 )
                    plane[IDX(mysize[_x_]+1, y)] += energy; // wrap around left edge to right edge
                if ( x == mysize[_x_] )
                    plane[IDX(0, y)] += energy; // wrap around right edge to left edge
                if ( y == 1 )
                    plane[IDX(x, mysize[_y_]+1)] += energy; // wrap around bottom edge to top edge
                if ( y == mysize[_y_] )
                    plane[IDX(x, 0)] += energy; // wrap around top edge to bottom edge
            }
    }
   #undef IDX
    
    return 0;
}
 
static inline int get_total_energy( const int       size[2],
                                    const double    *plane,
                                    double          *energy
                                )
/*
 * NOTE: this routine a good candiadate for openmp
 *       parallelization
 */
{

    const int register xsize = size[_x_];
    
   #define IDX( i, j ) ( (j)*(xsize+2) + (i) )

   #if defined(LONG_ACCURACY)    
    long double totenergy = 0;
   #else
    double totenergy = 0;    
   #endif

    // HINT: you may attempt to
    //       (i)  manually unroll the loop
    //       (ii) ask the compiler to do it
    // for instance
    // #pragma GCC unroll 4
    for ( int j = 1; j <= size[_y_]; j++ )
        for ( int i = 1; i <= size[_x_]; i++ )
            totenergy += plane[ IDX(i, j) ];
    
   #undef IDX

    *energy = (double)totenergy;
    return 0;
}
                            
static inline int memory_release (  double *data, 
                                    int *sources
                                )
  
{
  if( data != NULL )
    free( data );

  if( sources != NULL )
    free( sources );

  
  
  return 0;
}