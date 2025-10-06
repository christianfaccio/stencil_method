#include "parallel.h"

// ------------------------------------------------------------------

int main(int argc, char **argv)
{
  /* init */
  int  Niterations;
  int  periodic;
  vec2_t S,N;
  int Nsources, Nsources_local;
  vec2_t *Sources_local;
  double  energy_per_source;
  plane_t planes[2];
  buffers_t buffers[2];
  int output_energy_stat_perstep;
  int injection_frequency;

  /* MPI Init */
  int Rank, Nprocs;
  int neighbours[4];
  int mpi_provided_thread_level;
  MPI_Comm STENCIL_WORLD;
  
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_provided_thread_level);
  if (mpi_provided_thread_level < MPI_THREAD_FUNNELED) {
    printf("Aborting...\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MPI_Comm_dup(MPI_COMM_WORLD, &STENCIL_WORLD); // separate context

  MPI_Comm_rank(STENCIL_WORLD, &Rank);
  MPI_Comm_size(STENCIL_WORLD, &Nprocs);

  /* initial setting */
  int ret = initialize (  &STENCIL_WORLD, &Rank, &Nprocs, argc, argv, &S, &N, &periodic, &output_energy_stat_perstep,
                          neighbours, &Niterations,
                          &Nsources, &Nsources_local, &Sources_local, &energy_per_source,
                          planes, buffers, &injection_frequency );
  if (ret)
  {
    printf("task %d is opting out with termination code %d\n", Rank, ret );
    MPI_Finalize();
    return 0;
  }

  int current = OLD;
  double t0 = MPI_Wtime(); // for wall-clock time

  // Timing variables
  double total_comm_time = 0.0;
  double total_comp_time = 0.0;
  double t_start, t_end;

  /* main loop */
  for (int iter = 0; iter < Niterations; iter++)
    {

      MPI_Request reqs[8]; // requests for non-blocking comms

      // COMPUTATION: Energy injection
      t_start = MPI_Wtime();
      if ( iter % injection_frequency == 0)
      {
      	inject_energy( periodic, Nsources_local, Sources_local, energy_per_source, N, &planes[current] );
      }
      t_end = MPI_Wtime();
      total_comp_time += (t_end - t_start);

      /* -------------------------------------- */
      // [A] fill the buffers
      t_start = MPI_Wtime();
      double *data = planes[current].data;
      uint sizex = planes[current].size[_x_];
      uint sizey = planes[current].size[_y_];
      uint fsize = sizex + 2;  // frame size including ghost cells

      // NORTH: send first row (j=1)
      if (neighbours[NORTH] != MPI_PROC_NULL) {
          #pragma omp parallel for
          for (uint i = 0; i < sizex; i++)
              buffers[SEND][NORTH][i] = data[IDX(i+1, 1)];
      }

      // SOUTH: send last row (j=sizey)
      if (neighbours[SOUTH] != MPI_PROC_NULL) {
          #pragma omp parallel for
          for (uint i = 0; i < sizex; i++)
              buffers[SEND][SOUTH][i] = data[IDX(i+1, sizey)];
      }

      // EAST: send rightmost column (i=sizex)
      if (neighbours[EAST] != MPI_PROC_NULL) {
          #pragma omp parallel for
          for (uint j = 0; j < sizey; j++)
              buffers[SEND][EAST][j] = data[IDX(sizex, j+1)];
      }

      // WEST: send leftmost column (i=1)
      if (neighbours[WEST] != MPI_PROC_NULL) {
          #pragma omp parallel for
          for (uint j = 0; j < sizey; j++)
              buffers[SEND][WEST][j] = data[IDX(1, j+1)];
      }
      t_end = MPI_Wtime();
      total_comm_time += (t_end - t_start);

      // [B] perform the halo communications
      t_start = MPI_Wtime();
      //     (1) use Send / Recv
      //     (2) use Isend / Irecv
      //         --> can you overlap communication and compution in this way?
      
      int req_count = 0;
      for (int dir=0; dir<4; dir++){
	      if (neighbours[dir] != MPI_PROC_NULL){
		      uint bufsize = (dir < 2 ? sizex : sizey);
		      int send_tag = dir;
		      int recv_tag = (dir < 2 ? 1-dir : 5-dir); // opposite direction
		      MPI_Isend(	buffers[SEND][dir], bufsize, MPI_DOUBLE,
				      	neighbours[dir], send_tag, STENCIL_WORLD, &reqs[req_count++]
				);
		      MPI_Irecv(	buffers[RECV][dir], bufsize, MPI_DOUBLE,
				     	neighbours[dir], recv_tag, STENCIL_WORLD, &reqs[req_count++]
				);
	      }
      }
      MPI_Waitall(req_count, reqs, MPI_STATUSES_IGNORE);
      t_end = MPI_Wtime();
      total_comm_time += (t_end - t_start);

      // [C] copy the haloes data
      t_start = MPI_Wtime();

      // NORTH: copy received data to ghost row j=0
      if (neighbours[NORTH] != MPI_PROC_NULL) {
          #pragma omp parallel for
          for (uint i = 0; i < sizex; i++)
              data[IDX(i+1, 0)] = buffers[RECV][NORTH][i];
      }

      // SOUTH: copy received data to ghost row j=sizey+1
      if (neighbours[SOUTH] != MPI_PROC_NULL) {
          #pragma omp parallel for
          for (uint i = 0; i < sizex; i++)
               data[IDX(i+1, sizey+1)] = buffers[RECV][SOUTH][i];
      }

      // EAST: copy received data to ghost column i=sizex+1
      if (neighbours[EAST] != MPI_PROC_NULL) {
          #pragma omp parallel for
          for (uint j = 0; j < sizey; j++)
              data[IDX(sizex+1, j+1)] = buffers[RECV][EAST][j];
      }

      // WEST: copy received data to ghost column i=0
      if (neighbours[WEST] != MPI_PROC_NULL) {
          #pragma omp parallel for
          for (uint j = 0; j < sizey; j++)
              data[IDX(0, j+1)] = buffers[RECV][WEST][j];
      }
      t_end = MPI_Wtime();
      total_comm_time += (t_end - t_start);

      /* --------------------------------------  */
      /* update grid points */
      t_start = MPI_Wtime();
      update_plane( periodic, &planes[current], &planes[!current] );
      t_end = MPI_Wtime();
      total_comp_time += (t_end - t_start);

      /* output if needed */
      if ( output_energy_stat_perstep )
 	      output_energy_stat ( iter, &planes[!current], (iter+1)*Nsources*energy_per_source, Rank, Nprocs, S, N, &STENCIL_WORLD );
	
      /* swap plane indexes for the new iteration */
      current = !current;

    }

  double t1 = MPI_Wtime() - t0;

  output_energy_stat ( -1, &planes[!current], (Niterations/injection_frequency) * Nsources * energy_per_source, Rank, Nprocs, S, N, &STENCIL_WORLD );

  // Report timing statistics
  if ( Rank == 0 ) {
    printf("Elapsed time: %g seconds\n", t1);
    printf("\nTiming breakdown (Rank 0):\n");
    printf("  Computation time: %g seconds (%.2f%%)\n", total_comp_time, 100.0*total_comp_time/t1);
    printf("  Communication time: %g seconds (%.2f%%)\n", total_comm_time, 100.0*total_comm_time/t1);
    printf("  Other time: %g seconds (%.2f%%)\n", t1-total_comp_time-total_comm_time,
           100.0*(t1-total_comp_time-total_comm_time)/t1);
  }

  memory_release( buffers, planes, Sources_local );
  
  
  MPI_Finalize();
  return 0;
}

// ------------------------------------------------------------------

int initialize (  MPI_Comm  	*Comm,
                  int       	*Me,
                  int       	*Ntasks,
                  int		argc,
		  char   	**argv,
                  vec2_t    	*S,
                  vec2_t    	*N,
                  int     	*periodic,
                  int       	*output_energy_stat,
                  int       	*neighbours,
                  int     	*Niterations,
                  int     	*Nsources,
                  int       	*Nsources_local,
                  vec2_t 	**Sources_local,
                  double  	*energy_per_source,
                  plane_t 	*planes,
                  buffers_t 	*buffers,
                  int     	*injection_frequency
		)
{
  int halt = 0;     // halt flag
  int ret;          // return value
  int verbose = 0;  // verbosity level
  
  // ····························�····
  // set default values

  (*S)[_x_]         = 10000;
  (*S)[_y_]         = 10000;
  *periodic         = 0;
  *Nsources         = 1;
  *Nsources_local   = 0;
  *Sources_local    = NULL;
  *Niterations      = 100;
  *energy_per_source = 1.0;
  *output_energy_stat = 0;
  *injection_frequency = 10;
  
  planes[OLD].size[0] = 0;
  planes[OLD].size[1] = 0;
  planes[NEW].size[0] = 0;
  planes[NEW].size[1] = 0;
  planes[OLD].data = NULL;
  planes[NEW].data = NULL;
  
  neighbours[0] = MPI_PROC_NULL;
  neighbours[1] = MPI_PROC_NULL;
  neighbours[2] = MPI_PROC_NULL;
  neighbours[3] = MPI_PROC_NULL;

  buffers[0][0] = NULL;
  buffers[0][1] = NULL; 
  buffers[0][2] = NULL; 
  buffers[0][3] = NULL; 
  buffers[1][0] = NULL; 
  buffers[1][1] = NULL; 
  buffers[1][2] = NULL; 
  buffers[1][3] = NULL; 

  // ································
  // process the commadn line
  
  while ( 1 )
  {
    int opt;
    while((opt = getopt(argc, argv, ":hx:y:e:E:n:o:p:f:v:")) != -1)
      {
	switch( opt )
	  {
	  case 'x': (*S)[_x_] = (uint)atoi(optarg);
	    break;

	  case 'y': (*S)[_y_] = (uint)atoi(optarg);
	    break;

	  case 'e': *Nsources = atoi(optarg);
	    break;

	  case 'E': *energy_per_source = atof(optarg);
	    break;

	  case 'n': *Niterations = atoi(optarg);
	    break;

	  case 'o': *output_energy_stat = (atoi(optarg) > 0);
	    break;

	  case 'p': *periodic = (atoi(optarg) > 0);
	    break;

	  case 'f': *injection_frequency = atoi(optarg);
	    break;

	  case 'v': verbose = atoi(optarg);
	    break;

	  case 'h': {
	    if ( *Me == 0 )
	      printf( "\nvalid options are ( values btw [] are the default values ):\n"
		      "-x    x size of the plate [10000]\n"
		      "-y    y size of the plate [10000]\n"
		      "-e    how many energy sources on the plate [4]\n"
		      "-E    how many energy sources on the plate [1.0]\n"
		      "-n    how many iterations [100]\n"
		      "-p    whether periodic boundaries applies  [0 = false]\n"
		      "-f    injection frequency, i.e. every how many iterations to inject [10]\n\n"
		      );
	    halt = 1; }
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

  if ( halt )
    return 1;
  
  
  // ·······························��
  /*
   * here we should check for all the parms being meaningful
   */

  if ( (*S)[_x_] < 1 || (*S)[_y_] < 1 )
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
  if ( *output_energy_stat != 0 && *output_energy_stat != 1 )
    {
      printf("error: the output energy stat flag must be either 0 or 1\n");
      exit(1);
    }
  if ( *periodic != 0 && *periodic != 1 )
    {
      printf("error: the periodic flag must be either 0 or 1\n");
      exit(1);
    }
  if ( *injection_frequency < 1 || *injection_frequency > *Niterations )
    {
      printf("error: the injection frequency must be in [1,%d]\n", *Niterations);
      exit(1);
    }
  
  // ····························�···
  /*
   * find a suitable domain decomposition
   * very simple algorithm, you may want to
   * substitute it with a better one
   *
   * the plane Sx x Sy will be solved with a grid
   * of Nx x Ny MPI tasks
   */

  vec2_t Grid;
  double formfactor = ((*S)[_x_] >= (*S)[_y_] ? (double)(*S)[_x_]/(*S)[_y_] : (double)(*S)[_y_]/(*S)[_x_] );
  int    dimensions = 2 - (*Ntasks <= ((int)formfactor+1) ); 	// heuristics to avoid creating thin, inefficient
								// domain patches when having few procs
  
  if ( dimensions == 1 )
    {
      if ( (*S)[_x_] >= (*S)[_y_] )
	Grid[_x_] = *Ntasks, Grid[_y_] = 1;
      else
	Grid[_x_] = 1, Grid[_y_] = *Ntasks;
    }
  else
    {
      int   Nf;		// number of the factors (e.g. 3)
      uint  *factors; 	// array of factors (e.g. [2,2,3])
      uint  first = 1;
      ret = simple_factorization( *Ntasks, &Nf, &factors );
      
      for ( int i = 0; (i < Nf) && ((*Ntasks/first)/first > formfactor); i++ )
	first *= factors[i];

      if ( (*S)[_x_] > (*S)[_y_] )
	Grid[_x_] = *Ntasks/first, Grid[_y_] = first;
      else
	Grid[_x_] = first, Grid[_y_] = *Ntasks/first;
    }

  (*N)[_x_] = Grid[_x_];
  (*N)[_y_] = Grid[_y_];

  // ····························�···
  // my cooridnates in the grid of processors
  
  int X = *Me % Grid[_x_];
  int Y = *Me / Grid[_x_];

  // ····························�···
  // find my neighbours

  if ( Grid[_x_] > 1 )
    {  
      if ( *periodic ) {       
	neighbours[EAST]  = Y*Grid[_x_] + (*Me + 1 ) % Grid[_x_];
	neighbours[WEST]  = (X%Grid[_x_] > 0 ? *Me-1 : (Y+1)*Grid[_x_]-1); }
      
      else {
	neighbours[EAST]  = ( X < Grid[_x_]-1 ? *Me+1 : MPI_PROC_NULL );
	neighbours[WEST]  = ( X > 0 ? (*Me-1)%*Ntasks : MPI_PROC_NULL ); }  
    }

  if ( Grid[_y_] > 1 )
    {
      if ( *periodic ) {      
	neighbours[NORTH] = (*Ntasks + *Me - Grid[_x_]) % *Ntasks;
	neighbours[SOUTH] = (*Ntasks + *Me + Grid[_x_]) % *Ntasks; }

      else {    
	neighbours[NORTH] = ( Y > 0 ? *Me - Grid[_x_]: MPI_PROC_NULL );
	neighbours[SOUTH] = ( Y < Grid[_y_]-1 ? *Me + Grid[_x_] : MPI_PROC_NULL ); }
    }
  
  // ····························�····
  // the size of my patch

  /*
   * every MPI task determines the size sx x sy of its own domain
   * REMIND: the computational domain will be embedded into a frame
   *         that is (sx+2) x (sy+2)
   *         the outern frame will be used for halo communication
   */
  
  vec2_t mysize;
  uint s = (*S)[_x_] / Grid[_x_];
  uint r = (*S)[_x_] % Grid[_x_];
  mysize[_x_] = s + (X < r);
  s = (*S)[_y_] / Grid[_y_];
  r = (*S)[_y_] % Grid[_y_];
  mysize[_y_] = s + (Y < r);

  planes[OLD].size[0] = mysize[0];
  planes[OLD].size[1] = mysize[1];
  planes[NEW].size[0] = mysize[0];
  planes[NEW].size[1] = mysize[1];
  

  if ( verbose > 0 )
    { 
      if (*Me == 0)
      {
		printf("(DEBUG) -- Tasks are decomposed in a grid %d x %d\n\n", Grid[_y_], Grid[_x_] );
	     	fflush(stdout);
      }

      for ( int t = 0; t < *Ntasks; t++ )
      {   
	  MPI_Barrier(*Comm);
	  if ( t == *Me )
	  {
	      printf("(DEBUG) -- Task %d:\n"
		     "(DEBUG) -- grid coordinates : %d,%d\n"
		     "(DEBUG) -- grid size : %d,%d\n"
		     "(DEBUG) -- neighbours: N->%d, E->%d, S->%d, W->%d\n",
		     *Me, Y, X, mysize[1], mysize[0],
		     neighbours[NORTH], neighbours[EAST],
		     neighbours[SOUTH], neighbours[WEST] );
	      fflush(stdout);
	  }
      }
    }

  // ····························�···
  // allocae the needed memory

  ret = memory_allocate(neighbours, buffers, planes);
  
  // ····························�···
  // allocate the heat sources

  ret = initialize_sources( *Me, *Ntasks, Comm, mysize, *Nsources, Nsources_local, Sources_local );
 
  
  return 0; 
}

int update_plane (	const int     	periodic,
			const plane_t 	*oldplane,
                        plane_t 	*newplane
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
  register uint fsize = oldplane->size[_x_]+2;

  register uint xsize = oldplane->size[_x_];
  register uint ysize = oldplane->size[_y_];

  // HINT: you may attempt to
  //       (i)  manually unroll the loop
  //       (ii) ask the compiler to do it
  // for instance
  // #pragma GCC unroll 4
  //
  // HINT: in any case, this loop is a good candidate
  //       for openmp parallelization

  double *restrict old_data = oldplane->data;
  double *restrict new_data = newplane->data;
  const double alpha = 0.25;

  #pragma omp parallel for schedule(static)
  for (uint j = 1; j <= ysize; j++)
	{
	      #pragma GCC unroll 4
	      for ( uint i = 1; i <= xsize; i++)
		  {
		      // five-points stencil formula
		      double result = old_data[IDX(i,j)] * (1.0 - 4*alpha);
		      result += alpha * (old_data[IDX(i-1,j)] + old_data[IDX(i+1,j)] + old_data[IDX(i,j-1)] + old_data[IDX(i,j+1)]);

		      new_data[ IDX(i,j) ] = result;
		  }
	}

  if ( periodic )
      {
	  #pragma omp parallel for
	  for ( int i=1; i<=xsize; i++)
	  {
		new_data[ IDX(i, 0) ] = new_data[ IDX(i, ysize) ];
		new_data[ IDX(i, ysize+1) ] = new_data[ IDX(i, 1) ];
	  }

	  #pragma omp parallel for
	  for (int j=1; j<=ysize; j++)
	  {
		new_data[ IDX( 0, j) ] = new_data[ IDX(xsize, j) ];
		new_data[ IDX( xsize+1, j) ] = new_data[ IDX(1, j) ];
	  }
      }
  
  return 0;
}


