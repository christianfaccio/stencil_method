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
//
// function prototypes

int initialize (    int         argc,			        // argc from cmd line
		            char	    **argv,			        // argv from cmd line
		            int	        *S,			            // two-uint array defining x,y dims of grid
		            int	        *periodic,		        // periodic-boundary tag
		            int 	    *Niterations,		    // n. iterations
		            int   	    *Nsources,		        // n. heat sources
		            int  	    **Sources,              // the heat sources
		            double      *energy_per_source,	    // how much heat per source
		            double      **planes,		        // the two planes
                    int   	    *output_energy_at_steps,// whether to output the energy at every step
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

static inline int initialize_sources(   uint        size[2],
			                            int    	    Nsources,
			                            int  	    **Sources
                                    )
/*
 * randomly spread heat suources
 *
 */
{
  *Sources = (int*)malloc( Nsources * 2 *sizeof(uint) );
  for ( int s = 0; s < Nsources; s++ )
    {
      (*Sources)[s*2] = 1+ lrand48() % size[_x_];
      (*Sources)[s*2+1] = 1+ lrand48() % size[_y_];
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
{
   #define IDX( i, j ) ( (j)*(mysize[_x_]+2) + (i) )
    for (int s = 0; s < Nsources; s++) {
        
        int x = Sources[2*s];
        int y = Sources[2*s+1];
        plane[IDX(x, y)] += energy;

        if ( periodic )
            {
                if ( x == 1 )
                    plane[IDX(mysize[_x_]+1, y)] += energy;
                if ( x == mysize[_x_] )
                    plane[IDX(0, y)] += energy;
                if ( y == 1 )
                    plane[IDX(x, mysize[_y_]+1)] += energy;
                if ( y == mysize[_y_] )
                    plane[IDX(x, 0)] += energy;
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
    // an invalid pointer has been passed
    return 1;

  unsigned int bytes = (size[_x_]+2)*(size[_y_]+2);

  planes_ptr[OLD] = (double*)malloc( 2*bytes*sizeof(double) );
  memset ( planes_ptr[OLD], 0, 2*bytes*sizeof(double) );
  planes_ptr[NEW] = planes_ptr[OLD] + bytes;
      
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