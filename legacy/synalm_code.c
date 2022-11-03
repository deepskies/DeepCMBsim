fcomplex **fs_synalm(int nx,int ny,flouble lx,flouble ly,int nmaps,
		     nmt_k_function **cells,nmt_k_function **beam,int seed)/* function which returns float double (?) that takes in: # of pixels in x direction, # of pixels in y direction, degrees in x direction, degrees in y direction, # of fields to be generated, cls, beam patters, random seed*/

{
  int imap; //map to be looped over
  fcomplex **alms; //defining alms array
  alms=my_malloc(nmaps*sizeof(fcomplex *)); //allocate memory for alms array
  for(imap=0;imap<nmaps;imap++) //allocate memory alms data in said array
    alms[imap]=dftw_malloc(ny*(nx/2+1)*sizeof(fcomplex));

  //Switch off error handler for Cholesky decomposition //I don't know what that means
  gsl_error_handler_t *geh=gsl_set_error_handler_off();

  int numthr=0;

#pragma omp parallel default(none)			\ // I think this is for parallel processors?
  shared(nx,ny,lx,ly,nmaps,cells,beam,seed,alms,numthr)
  {
    //This is to avoid using the omp.h library
    int ithr;
#pragma omp critical
    {
      ithr=numthr;
      numthr++;
    }

    int iy;
    double dkx=2*M_PI/lx,dky=2*M_PI/ly;
    double inv_dkvol=1./(dkx*dky);
    gsl_vector *rv1=gsl_vector_alloc(nmaps);
    gsl_vector *iv1=gsl_vector_alloc(nmaps);
    gsl_vector *rv2=gsl_vector_alloc(nmaps);
    gsl_vector *iv2=gsl_vector_alloc(nmaps);
    gsl_matrix *clmat=gsl_matrix_calloc(nmaps,nmaps); 
    gsl_vector *eval =gsl_vector_alloc(nmaps);
    gsl_matrix *evec =gsl_matrix_alloc(nmaps,nmaps); 
    gsl_eigen_symmv_workspace *wsym=gsl_eigen_symmv_alloc(nmaps);
    unsigned int seed_thr=(unsigned int)(seed+ithr);
    gsl_rng *rng=init_rng(seed_thr);
    gsl_interp_accel *intacc_cells=gsl_interp_accel_alloc();
    gsl_interp_accel *intacc_beam=gsl_interp_accel_alloc();

#pragma omp for
    for(iy=0;iy<ny;iy++) { //looping over each column
      int ix; 
      flouble ky;
      if(2*iy<=ny)
	ky=iy*dky; //defining associated y-frequency for each column, positive on left side
      else
	ky=-(ny-iy)*dky; //negative frequency on right side of map
      for(ix=0;ix<=nx/2;ix++) { //looping over each pixel in each column
	int imp1,imp2;
	flouble kx=ix*dkx; //x-frequency for each pixel
	long index=ix+(nx/2+1)*iy;
	flouble kmod=sqrt(kx*kx+ky*ky); //overall frequency for pixel
	if(kmod<0) {
	  for(imp1=0;imp1<nmaps;imp1++)
	    alms[imp1][index]=0; //0 amplitude if the overall frequency is negative
	}
	else {
	  //Get power spectrum
	  int icl=0; //C_l being iterated over?
	  for(imp1=0;imp1<nmaps;imp1++) { //iterating over maps
	    for(imp2=imp1;imp2<nmaps;imp2++) {
	      //Fill up only lower triangular part //not sure what they mean by that
	      
	      flouble cl=0.5*inv_dkvol*nmt_k_function_eval(cells[icl],kmod,intacc_cells);
	      //evaluation function really just splines the data and gives an input
	      
	      gsl_matrix_set(clmat,imp1,imp2,cl); //setting data point in matrix
	      if(imp2!=imp1)
		gsl_matrix_set(clmat,imp2,imp1,cl); //mirroring point
	      icl++;
	    }
	  }

	  //Take square root
	  gsl_eigen_symmv(clmat,eval,evec,wsym); //calculates eigenvalues
	  for(imp1=0;imp1<nmaps;imp1++) {
	    double dr,di; //At the same time get white random numbers
	    rng_gauss(rng,&dr,&di); //generating some random numbers
	    gsl_vector_set(rv1,imp1,dr); //setting vector ranges?
	    gsl_vector_set(iv1,imp1,di);
	    for(imp2=0;imp2<nmaps;imp2++) { //looping over maps
	      double oij=gsl_matrix_get(evec,imp1,imp2); //setting matrix for randomness?
	      double lambda=gsl_vector_get(eval,imp2); //setting vector to add some randomness?
	      if(lambda<=0) lambda=0;
	      else lambda=sqrt(lambda);
	      gsl_matrix_set(clmat,imp1,imp2,oij*lambda); //adding some randomnes to Cls matrix?
	    }
	  }

	  //Get correlate random numbers
	  gsl_blas_dgemv(CblasNoTrans,1.,clmat,rv1,0,rv2); //I don't know what this does
	  gsl_blas_dgemv(CblasNoTrans,1.,clmat,iv1,0,iv2);
	  for(imp1=0;imp1<nmaps;imp1++) {
	    flouble bm=nmt_k_function_eval(beam[imp1],kmod,intacc_beam);//getting beam for map
	    flouble a_re=bm*gsl_vector_get(rv2,imp1);//real and imaginary parts of random numbers?
	    flouble a_im=bm*gsl_vector_get(iv2,imp1);//seems like amplitude is the beam against the vector containing the randomly generated number
	    if(ix==0) {
	      if(iy>ny/2)
		continue;
	      else {
		if(iy==0)
		  alms[imp1][index]=(fcomplex)(M_SQRT2*a_re+I*0*a_im);//0,0 pixel must be real
		else {
		  int iyy=ny-iy;
		  alms[imp1][index]=(fcomplex)(a_re+I*a_im);//setting equal to value
		  alms[imp1][ix+(nx/2+1)*iyy]=(fcomplex)(a_re-I*a_im);//flipping about axis (symmetry)
		}
	      }
	    }
	    else
	      alms[imp1][index]=(fcomplex)(a_re+I*a_im);
	  }
	}
      }
    } //omp end for
    gsl_vector_free(rv1);
    gsl_vector_free(iv1);
    gsl_vector_free(rv2);
    gsl_vector_free(iv2);
    gsl_matrix_free(clmat);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_eigen_symmv_free(wsym);
    end_rng(rng);
    gsl_interp_accel_free(intacc_cells);
    gsl_interp_accel_free(intacc_beam);
  } //omp end parallel

  //Restore error handler
  gsl_set_error_handler(geh);

  return alms;
}
