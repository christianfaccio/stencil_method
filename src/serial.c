#include "serial.h"

// ------------------------------------------------------------------

int main(int argc, char **argv)
{
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

int initialize (  int		    argc,                   
		              char   	  **argv,                	
                  int     	*S,                   	
                  int     	*periodic,            	
                  int     	*Niterations,         	
                  int     	*Nsources,            	
                  int    	  **Sources,              
                  double  	*energy_per_source,   	
                  double 	  **planes,             	
                  int     	*output_energy_at_steps,
                  int     	*injection_frequency    
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
                  const int     size[2],
			            const double  *old,
                  double        *new
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
    const int register fxsize = size[_x_]+2;
    const int register fysize = size[_y_]+2;
    const int register xsize = size[_x_];
    const int register ysize = size[_y_];
    
   #define IDX( i, j ) ( (j)*fxsize + (i) )

    // HINT: you may attempt to
    //       (i)  manually unroll the loop
    //       (ii) ask the compiler to do it
    // for instance
    // #pragma GCC unroll 4
    //
    // HINT: in any case, this loop is a good candidate
    //       for openmp parallelization
    for (int j = 1; j <= ysize; j++)
        for ( int i = 1; i <= xsize; i++)
            {
                //
                // five-points stencil formula
                //

                
                // simpler stencil with no explicit diffusivity
                // always conserve the smoohed quantity
                // alpha here mimics how much "easily" the heat
                // travels
                
                double alpha = 0.6;
                double result = old[ IDX(i,j) ] *alpha;
                double sum_i  = (old[IDX(i-1, j)] + old[IDX(i+1, j)]) / 4.0 * (1-alpha);
                double sum_j  = (old[IDX(i, j-1)] + old[IDX(i, j+1)]) / 4.0 * (1-alpha);
                result += (sum_i + sum_j);
                

                /*

                  // implentation from the derivation of
                  // 3-points 2nd order derivatives
                  // however, that should depends on an adaptive
                  // time-stepping so that given a diffusivity
                  // coefficient the amount of energy diffused is
                  // "small"
                  // however the imlic methods are not stable
                  
               #define alpha_guess 0.5     // mimic the heat diffusivity

                double alpha = alpha_guess;
                double sum = old[IDX(i,j)];
                
                int   done = 0;
                do
                    {                
                        double sum_i = alpha * (old[IDX(i-1, j)] + old[IDX(i+1, j)] - 2*sum);
                        double sum_j = alpha * (old[IDX(i, j-1)] + old[IDX(i, j+1)] - 2*sum);
                        result = sum + ( sum_i + sum_j);
                        double ratio = fabs((result-sum)/(sum!=0? sum : 1.0));
                        done = ( (ratio < 2.0) && (result >= 0) );    // not too fast diffusion and
                                                                     // not so fast that the (i,j)
                                                                     // goes below zero energy
                        alpha /= 2;
                    }
                while ( !done );
                */

                new[ IDX(i,j) ] = result;
                
            }

    if ( periodic )
        /*
         * propagate boundaries if they are periodic
         *
         * NOTE: when is that needed in distributed memory, if any?
         */
        {
            for ( int i = 1; i <= xsize; i++ )
                {
                    new[ i ] = new[ IDX(i, ysize) ];
                    new[ IDX(i, ysize+1) ] = new[ i ];
                }
            for ( int j = 1; j <= ysize; j++ )
                {
                    new[ IDX( 0, j) ] = new[ IDX(xsize, j) ];
                    new[ IDX( xsize+1, j) ] = new[ IDX(1, j) ];
                }
        }
    
    return 0;

   #undef IDX
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
