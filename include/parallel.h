#define _XOPEN_SOURCE 600
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>


#define NORTH 	0
#define SOUTH 	1
#define EAST  	2
#define WEST  	3

#define SEND 	0
#define RECV 	1

#define OLD 	0
#define NEW 	1

#define _x_ 	0
#define _y_ 	1

#define IDX(i,j) ((j)*fsize + (i))
#define LONG_ACCURACY 0

typedef unsigned int uint;
typedef uint    vec2_t[2]; 		// for (x,y) coordinates
typedef double 	*restrict buffers_t[4]; // for MPI communication buffers (NORTH, SOUTH, EAST, WEST)

typedef struct {
    double	*restrict data;
    vec2_t  	size;
} plane_t; 				// for the 2D grid (with ghost cells)

// ============================================================

int initialize (	MPI_Comm 	*Comm,
                    	int       	*Me,			// my rank
                    	int       	*Ntasks,               	// total number of MPI ranks
                    	int         	argc,		 	// argc from cmd line
                    	char	    	**argv,	  		// argv from cmd line
                    	vec2_t	   	*S,        	   	// two-uint array defining (x,y) dims of grid
                    	vec2_t   	*N,                     // two-uint array defining the MPI tasks' grid 
								// dims (x,y)
                    	int	        *periodic,		// boolean flag for boundary conditions - if 1, 
                                                            	// the plate behaves as infinite with periodic 
                                                            	// boundaries; if 0, fixed boundaries
                    	int       	*output_energy_stat,   	// debug flag - if 1, prints energy statistics 
                                                        	// and dumps grid data at each time step
                    	int       	*neighbours,           	// four-int array that gives back the neighbours 
								// of the calling task 
                    	int 	    	*Niterations,		// n. iterations
                    	int   	    	*Nsources,		// n. heat sources randomly placed on the grid 
                                                            	// (default = 1)
                    	int   	    	*Nsources_local,       	// n. heat sources for this MPI task
                    	vec2_t  	**Sources_local,        // 2D array storing heat source positions;
                                                            	// sources[2*i] and Sources[2*i+1] are x,y 
                                                            	// coordinates of source i
                    	double      	*energy_per_source,	// amount of energy each heat source injects 
                                                            	// per injection event
                    	plane_t      	*planes,		// two planes: planes[OLD] contains the current data,
                                                            	// planes[NEW] will contain the updated data
                                                            	// the two planes swap their roles at every iteration
                    	buffers_t    	*buffers,            	// communication buffers for the four directions
                    	int   	*injection_frequency    	// how often to inject energy (every how many iters)
		);  

int update_plane (	const int       periodic,
			const plane_t   *oldplane,
                        plane_t   	*newplane
               	);   

// ============================================================

static inline int memory_allocate ( 	const int	*neighbours,
		                        buffers_t       *buffers_ptr,
		                        plane_t         *planes_ptr
                               	)
/*
  here you allocate the memory buffers that you need to
  (i)  hold the results of your computation
  (ii) communicate with your neighbours

  The memory layout that I propose to you is as follows:

  (i) --- calculations
  you need 2 memory regions: the "OLD" one that contains the
  results for the step (i-1)th, and the "NEW" one that will contain
  the updated results from the step ith.

  Then, the "NEW" will be treated as "OLD" and viceversa.

  These two memory regions are indexed by *planes_ptr:

  planes_ptr[0] ==> the "OLD" region
  planes_ptr[1] ==> the "NEW" region


  (ii) --- communications

  you may need two buffers (one for sending and one for receiving)
  for each one of your neighbours, that are at most 4:
  north, south, east and west.

  To them you need to communicate at most mysizex or mysizey
  double data.

  These buffers are indexed by the buffers_ptr pointer so
  that

  (*buffers_ptr)[SEND][ {NORTH,...,WEST} ] = .. some memory regions
  (*buffers_ptr)[RECV][ {NORTH,...,WEST} ] = .. some memory regions
  
  --->> Of course you can change this layout as you prefer
  */
{
  if (planes_ptr == NULL )  
    return 1;
  if (buffers_ptr == NULL )
    return 2;

  /*
  use +2 for ghost cells:
  Ghost cells:    G G G G G G
                  G x x x x G  <- Real grid points (x)
                  G x x x x G     surrounded by ghost cells (G)
                  G x x x x G
                  G x x x x G
                  G G G G G G
  */ 
  unsigned int frame_size = (planes_ptr[OLD].size[_x_]+2) * (planes_ptr[OLD].size[_y_]+2);

  planes_ptr[OLD].data = (double*)calloc( 2*frame_size, sizeof(double));// Allocate one big block using calloc
                                                                        // such that values are initialized to 0
  if ( planes_ptr[OLD].data == NULL )
    return 3;

  planes_ptr[NEW].data = planes_ptr[OLD].data + frame_size; 		// points to second half
  if ( planes_ptr[NEW].data == NULL )
    return 4;

  // ··················································
  // buffers for north and south communication 
  // are not really needed
  //
  // in fact, they are already contiguous, just the
  // first and last line of every rank's plane
  //
  // you may just make some pointers pointing to the
  // correct positions
  //

  // or, if you prefer, just go on and allocate buffers
  // also for north and south communications

  // ··················································
  // allocate buffers
  //
  for ( int sr=0; sr < 2; sr++)
  	for ( int dir = 0; dir < 4; dir++ )
    	{
      		if ( neighbours[dir] != MPI_PROC_NULL )
        	{
          		unsigned int bufsize = (dir < 2? planes_ptr[OLD].size[_x_] : planes_ptr[OLD].size[_y_]);
          		buffers_ptr[sr][dir] = (double*)malloc( bufsize * sizeof(double) );
          		if ( buffers_ptr[sr][dir] == NULL )
            			return 5;
        	}
      		else
        	{
          		buffers_ptr[sr][dir] = NULL;          
        	}
    	}

  // ··················································
  return 0;
}

static inline int initialize_sources(   int     	Me,
                                        int     	Ntasks,
                                        MPI_Comm 	*Comm,
                                        vec2_t  	mysize,
			                int    		Nsources,
                                        int     	*Nsources_local,
			                vec2_t  	**Sources
                                    )
{
  srand48(time(NULL) ^ Me);	// seed different for every task
  int *tasks_with_sources = (int*)malloc( Nsources * sizeof(int) );
  
  if ( Me == 0 )
    {
      for ( int i = 0; i < Nsources; i++ )
	      tasks_with_sources[i] = (int)lrand48() % Ntasks;
    }
  
  MPI_Bcast( tasks_with_sources, Nsources, MPI_INT, 0, *Comm );
  
  int nlocal = 0;
  for ( int i = 0; i < Nsources; i++ )
    nlocal += (tasks_with_sources[i] == Me);
  *Nsources_local = nlocal;
  
  if ( nlocal > 0 )
    {
      vec2_t *restrict helper = (vec2_t*)malloc( nlocal * sizeof(vec2_t) );      
      for ( int s = 0; s < nlocal; s++ )
        {
          helper[s][_x_] = 1 + lrand48() % mysize[_x_];
          helper[s][_y_] = 1 + lrand48() % mysize[_y_];
        }
      *Sources = helper;
    }
  
  free( tasks_with_sources );

  return 0;
}

static inline int inject_energy (   const int       	periodic,
                                    const int       	Nsources,
			            const vec2_t   	*Sources,
			            const double    	energy,
			            const vec2_t       	N,
                                    plane_t      	*plane
                                )
/*
Function to inject energy into the grid at the positions specified in Sources array.
If periodic is true, it also updates the ghost cells accordingly.

Grid with periodic boundaries and source at (1,3):

[G][G][G][G][G]
[G][S][x][x][S] <- Source at (1,3) injects here
[G][x][x][x][G]
[G][x][x][x][G]
[G][S][G][G][G]
*/
{
  register const uint fsize = plane->size[_x_]+2;
  double *restrict data = plane->data;
  
  for (int s = 0; s < Nsources; s++)
  {
      int x = Sources[s][_x_];
      int y = Sources[s][_y_];
      data[IDX(x, y)] += energy;

      // Only update ghost cells if periodic AND single process in that dimension
      // (otherwise MPI will handle the periodic communication)
      if (periodic) {
	  // If only 1 process in x-direction, handle x-wrapping locally
	  if ( N[_x_] == 1 ) {
	      if ( x == 1 )
		  data[IDX(plane->size[_x_]+1, y)] += energy;
	      if ( x == plane->size[_x_] )
		  data[IDX(0, y)] += energy;
	  }

	  // If only 1 process in y-direction, handle y-wrapping locally
	  if ( N[_y_] == 1 ) {
	      if ( y == 1 )
		  data[IDX(x, plane->size[_y_]+1)] += energy;
	      if ( y == plane->size[_y_] )
		  data[IDX(x, 0)] += energy;
	  }
      }
  }

  return 0;
}
 
static inline uint simple_factorization(  uint 	A,
                                          int 	*Nfactors, 
                                          uint 	**factors 
					)
/*
 * rought factorization;
 * assumes that A is small, of the order of <~ 10^5 max,
 * since it represents the number of tasks
 */
{
  int N = 0; 	// counter for the number of factors found
  int f = 2;	// factor currently being tested
  uint _A_ = A;	// used to preserve original A since we do two passes

  while ( f*f <= _A_ )
    {
      while( _A_ % f == 0 ) {
	N++;
	_A_ /= f; }

      f++;
    }
  if ( _A_ > 1 )
	  N++;

  *Nfactors = N;
  uint *_factors_ = (uint*)malloc( N * sizeof(uint) );

  N   = 0;
  f   = 2;
  _A_ = A;

  while ( f*f <= _A_ )
    {
      while( _A_ % f == 0 ) {
	_factors_[N++] = f;
	_A_ /= f; }
      f++;
    }
  if ( _A_ > 1 )
	  _factors_[N++] = _A_;

  *factors = _factors_;
  return 0;
}

static inline int get_total_energy( plane_t *plane,
                                    double  *energy
                                )
{	
    	register const int xsize = plane->size[_x_];
    	register const int ysize = plane->size[_y_];
    	register const int fsize = xsize+2;

    	double * restrict data = plane->data;
    	
   	#if defined(LONG_ACCURACY)    
    	long double totenergy = 0;
   	#else
    	double totenergy = 0;    
   	#endif

	#pragma omp parallel for schedule(static) reduction(+:totenergy)
	for ( int j = 1; j <= ysize; j++ )
		#pragma GCC unroll 4
		for ( int i = 1; i <= xsize; i++ )
			totenergy += data[ IDX(i, j) ];

	*energy = (double)totenergy;
	return 0;
}
                            
static inline int output_energy_stat (  int 		step, 
                          		plane_t 	*plane, 
                          		double 		budget, 
                          		int 		Me, 
                          		MPI_Comm 	*Comm 
					)
{

  double system_energy = 0;
  double tot_system_energy = 0;
  get_total_energy ( plane, &system_energy );
  
  MPI_Reduce ( &system_energy, &tot_system_energy, 1, MPI_DOUBLE, MPI_SUM, 0, *Comm );
  
  if ( Me == 0 )
    {
      if ( step >= 0 )
	printf(" [ step %4d ] ", step );
      fflush(stdout);

      
      printf( "total injected energy is %g, "
	      "system energy is %g "
	      "( in avg %g per grid point)\n",
	      budget,
	      tot_system_energy,
	      tot_system_energy / (plane->size[_x_]*plane->size[_y_]) );
    }
  
  return 0;
}

static inline int memory_release (  	buffers_t	*buffers,
					plane_t   	*planes,
                                    	vec2_t   	*sources
                                )
{
  if (buffers != NULL) {
    for ( int sr = 0; sr < 2; sr++ )
      for ( int dir = 0; dir < 4; dir++ )
        if ( buffers[sr][dir] != NULL )
          free( buffers[sr][dir] );
  }
  if( planes != NULL )
    	free( planes[OLD].data );
  if( sources != NULL )
    	free( sources );
  return 0;
}
