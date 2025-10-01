#include "serial.h"

// ------------------------------------------------------------------

int main(int argc, char **argv)
{
  /* MPI Init */
  int Rank, Nprocs;
  int mpi_provided_thread_level;
  MPI_Comm STENCIL_WORLD;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &mpi_provided_thread_level); // !!! to be set for multithreading for Leonardo
  if (provided < MPI_THREAD_SINGLE) {
    printf("Aborting...")
    MPI_Abort(STENCIL_WORLD, 1)
  }
  
  MPI_Comm_dup(MPI_COMM_WORLD, &STENCIL_WORLD); // separate context

  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
  MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);

  MPI_Barrier(STENCIL_WORLD);

  /* init */
  int  Niterations;
  int  periodic;
  int  S[2];
  
  int     Nsources;
  int    *Sources;
  double  energy_per_source;

  double *planes[2];
  
  double injected_heat = 0;

  int injection_frequency;
  int output_energy_at_steps = 0;
   
  initialize ( argc, argv, &S[0], &periodic, &Niterations,
	       &Nsources, &Sources, &energy_per_source, &planes[0],
	       &output_energy_at_steps, &injection_frequency );
  
  int current = OLD;

  /* main loop */
  for (int iter = 0; iter < Niterations; iter++)
    {
      /* (opt) inject energy */
      if ( iter % injection_frequency == 0 )
        {
          inject_energy( periodic, Nsources, Sources, energy_per_source, S, planes[current] );
          injected_heat += Nsources*energy_per_source;
        }
                  
      /* update grid points */
      update_plane(periodic, S, planes[current], planes[!current] );

      if ( output_energy_at_steps )
        {
          double system_heat;
          get_total_energy( S, planes[!current], &system_heat);
                  
          printf("step %d :: injected energy is %g, updated system energy is %g\n", iter, 
          injected_heat, system_heat );

          char filename[100];
          sprintf( filename, "plane_%05d.bin", iter );
          dump( planes[!current], S, filename, NULL, NULL );
            
        }

      /* swap planes for the new iteration */
      current = !current;
      
    }
  
  /* get final heat in the system */
  double system_heat;
  get_total_energy( S, planes[current], &system_heat);

  printf("injected energy is %g, system energy is %g\n",
	 injected_heat, system_heat );
  
  memory_release( planes[OLD], Sources );

  MPI_Finalize();
  return 0;
}

// ------------------------------------------------------------------

int initialize (  MPI_Comm  *Comm,
                  int       *Me,
                  int       *Ntasks,
                  int		    argc,                   
		              char   	  **argv,                	
                  vec2_t    *S,                   
                  vec2_t    *N,	
                  int     	*periodic,   
                  int       *output_energy_stat,  
                  int       *neighbours,       	
                  int     	*Niterations,         	
                  int     	*Nsources,   
                  int       *Nsources_local         	
                  vec2_t 	  **Sources_local,              
                  double  	*energy_per_source,   	
                  plane_t 	*planes,          
                  buffers_t *buffers,   	
                  //int     	*output_energy_at_steps,
                  //int     	*injection_frequency    
		            )
{
  int ret;  // return value
  
  // ··································································
  // set default values

  S[_x_]            = 1000;
  S[_y_]            = 1000;
  *periodic         = 0;
  *Nsources         = 1;
  *Niterations      = 99;
  *output_energy_at_steps = 0;
  *energy_per_source = 1.0;
  *injection_frequency = *Niterations;

  double freq = 0; 
  
  // ··································································
  // process the commadn line
  // 
  while ( 1 )
  {
    int opt;
    while((opt = getopt(argc, argv, ":x:y:e:E:f:n:p:o:")) != -1)
      {
	switch( opt )
	  {
	  case 'x': S[_x_] = (uint)atoi(optarg);
	    break;

	  case 'y': S[_y_] = (uint)atoi(optarg);
	    break;

	  case 'e': *Nsources = atoi(optarg);
	    break;

	  case 'E': *energy_per_source = atof(optarg);
	    break;

	  case 'n': *Niterations = atoi(optarg);
	    break;

	  case 'p': *periodic = (atoi(optarg) > 0);
	    break;

	  case 'o': *output_energy_at_steps = (atoi(optarg) > 0);
	    break;

	  case 'f': freq = atof(optarg);
	    break;
	    
	  case 'h': printf( "valid options are ( values btw [] are the default values ):\n"
			    "-x    x size of the plate [1000]\n"
			    "-y    y size of the plate [1000]\n"
			    "-e    how many energy sources on the plate [1]\n"
			    "-E    how many energy sources on the plate [1.0]\n"
			    "-f    the frequency of energy injection [0.0]\n"
			    "-n    how many iterations [100]\n"
			    "-p    whether periodic boundaries applies  [0 = false]\n"
			    "-o    whether to print the energy budgest at every step [0 = false]\n"
			    );
	    break;
	    

	  case ':': printf( "option -%c requires an argument\n", optopt);
	    break;

	  case '?': printf(" -------- help unavailable ----------\n");
	    break;
	  }
      }

    if ( opt == -1 )
      break;
  }

  if ( freq == 0 )
    *injection_frequency = *Niterations;
  else
    {
      *injection_frequency = *Niterations / (int)freq;
    }

  // ··································································
  /*
   * we should check for all the parms being meaningful
   *
   */

  if ( S[_x_] < 1 || S[_y_] < 1 )
    {
      printf("error: the size of the plate must be positive\n");
      exit(1);
    }
  if ( *Niterations < 1 )
    {
      printf("error: the number of iterations must be positive\n");
      exit(1);
    }
  if ( *Nsources < 1 )
    {
      printf("error: the number of heat sources must be positive\n");
      exit(1);
    }
  if ( *energy_per_source <= 0 )
    {
      printf("error: the energy per source must be positive\n");
      exit(1);
    }
  if ( *injection_frequency < 1 || *injection_frequency > *Niterations )
    {
      printf("error: the injection frequency must be in [1,%d]\n", *Niterations);
      exit(1);
    }
  

  // ··································································
  // allocae the needed memory
  //
  ret = memory_allocate( S, planes );
  

  // ··································································
  // allocae the heat sources
  //
  ret = initialize_sources( S, *Nsources, Sources );
  
  
  return 0;  
}

int update_plane (const int     periodic, 
                  const vec2_t  N,
			            const plane_t *old,
                        plane_t *new
                  )
/*
 * calculate the new energy values
 * the old plane contains the current data, the new plane
 * will store the updated data
 *
 * NOTE: in parallel, every MPI task will perform the
 *       calculation for its patch
 *
 */
{
  uint register fxsize = oldplane->size[_x_]+2;
  uint register fysize = oldplane->size[_y_]+2;
  
  uint register xsize = oldplane->size[_x_];
  uint register ysize = oldplane->size[_y_];
    
  #define IDX( i, j ) ( (j)*fxsize + (i) )
  
  // HINT: you may attempt to
  //       (i)  manually unroll the loop
  //       (ii) ask the compiler to do it
  // for instance
  // #pragma GCC unroll 4
  //
  // HINT: in any case, this loop is a good candidate
  //       for openmp parallelization

  double * restrict old = oldplane->data;
  double * restrict new = newplane->data;
  
  for (uint j = 1; j <= ysize; j++)
      for ( uint i = 1; i <= xsize; i++)
          {
              
              // NOTE: (i-1,j), (i+1,j), (i,j-1) and (i,j+1) always exist even
              //       if this patch is at some border without periodic conditions;
              //       in that case it is assumed that the +-1 points are outside the
              //       plate and always have a value of 0, i.e. they are an
              //       "infinite sink" of heat
              
              // five-points stencil formula
              //
              // HINT : check the serial version for some optimization
              //
              new[ IDX(i,j) ] =
                  old[ IDX(i,j) ] / 2.0 + ( old[IDX(i-1, j)] + old[IDX(i+1, j)] +
                                            old[IDX(i, j-1)] + old[IDX(i, j+1)] ) /4.0 / 2.0;
              
          }

  if ( periodic )
      {
          if ( N[_x_] == 1 )
              {
                  // propagate the boundaries as needed
                  // check the serial version
              }

          if ( N[_y_] == 1 ) 
              {
                  // propagate the boundaries as needed
                  // check the serial version
              }
      }

  
  #undef IDX
  return 0;
}

int dump (        const double  *data, 
                  const uint    size[2], 
                  const char    *filename, 
                  double        *min, 
                  double        *max
          )
/* dump the data in a binary file, in single precision */
{
  if ( (filename != NULL) && (filename[0] != '\0') )
    {
      FILE *outfile = fopen( filename, "w" );
      if ( outfile == NULL )
	return 2;
      
      float *array = (float*)malloc( size[0] * sizeof(float) );
      
      double _min_ = DBL_MAX;
      double _max_ = 0;

      for ( int j = 0; j < size[1]; j++ )
	{
	  /*
	  float y = (float)j / size[1];
	  fwrite ( &y, sizeof(float), 1, outfile );
	  */
	  
	  const double * restrict line = data + j*size[0];
	  for ( int i = 0; i < size[0]; i++ ) {
	    array[i] = (float)line[i];
	    _min_ = ( line[i] < _min_? line[i] : _min_ );
	    _max_ = ( line[i] > _max_? line[i] : _max_ ); }
	  
	  fwrite( array, sizeof(float), size[0], outfile );
	}
      
      free( array );
      
      fclose( outfile );

      if ( min != NULL )
	*min = _min_;
      if ( max != NULL )
	*max = _max_;
    }

  else return 1;
  
}
