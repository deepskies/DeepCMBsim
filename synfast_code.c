flouble **nmt_synfast_flat(int nx,int ny,flouble lx,flouble ly,int nfields,int *spin_arr,
			   int nl_beam,flouble *l_beam,flouble **beam_fields,
			   int nl_cell,flouble *l_cell,flouble **cell_fields,
			   int seed) /* function which returns float double (?)
that takes in: # of pixels in x direction, # of pixels in y direction, degrees in x direction, degrees in y direction, # of fields to be generated, spins of these fields, number of beam patterns for different l modes (???), corresponding beam patterns for these l modes (???), fields associated with these beam patters (???), number of cell patterns for different l modes (???), correspoinding cell patterns (???) */
{
  int ifield,imap; //represents current field or map being looped over
  int nmaps=0,ncls=0; //# of maps, # of spectra
  long npix=nx*ny; //# of pixels in map
  nmt_k_function **beam,**cell; //creating namaster function pointer
  flouble **maps; //float double
  fcomplex **alms; //float doubles for alms
  nmt_flatsky_info *fs=nmt_flatsky_info_alloc(nx,ny,lx,ly); //function to get info about map params?
  for(ifield=0;ifield<nfields;ifield++) { // looping over different fields
    int nmp=1;
    if(spin_arr[ifield]) nmp=2;
    nmaps+=nmp; //increase the number of maps for the number of fields?
  }

  imap=0; //map being looped over
  beam=my_malloc(nmaps*sizeof(nmt_k_function *)); //allocating storage for beam patterns?
  maps=my_malloc(nmaps*sizeof(flouble *)); //allocating storage for maps?
  for(ifield=0;ifield<nfields;ifield++) { //looping over fields (again? why again?)
    int imp,nmp=1; //still not sure what these variables are for...
    if(spin_arr[ifield]) nmp=2;
    for(imp=0;imp<nmp;imp++) {
      beam[imap+imp]=nmt_k_function_alloc(nl_beam,l_beam,beam_fields[ifield],beam_fields[ifield][0],0.,0); //defining beam pattern for a given map?
      maps[imap+imp]=dftw_malloc(npix*sizeof(flouble));//allocating storage for maps of size npix
    }
    imap+=nmp;
  }

  ncls=(nmaps*(nmaps+1))/2; //# of spectra (aka C_l's/cls, given the # of maps)
  cell=my_malloc(ncls*sizeof(nmt_k_function *)); //allocating data for spectra/cls
  for(imap=0;imap<ncls;imap++) // looping over maps
    cell[imap]=nmt_k_function_alloc(nl_cell,l_cell,cell_fields[imap],cell_fields[imap][0],0.,0);//allocating storage for spectrum/cls for the map being looped over

  alms=fs_synalm(nx,ny,lx,ly,nmaps,cell,beam,seed); //calculating the amplitude spectra for all the maps given number of pixels on each side, degrees on each side, maps, spectra/cls, beam pattern, and random seed

  for(imap=0;imap<nmaps;imap++) 
    nmt_k_function_free(beam[imap]);//I don't know what this does
  free(beam);//does this function de-allocate "beam"?
  for(imap=0;imap<ncls;imap++)
    nmt_k_function_free(cell[imap]);//I don't know what this does
  free(cell);//does this function de-allocate "cell"?

  imap=0;
  for(ifield=0;ifield<nfields;ifield++) {//looping over maps 
    int imp,nmp=1;
    if(spin_arr[ifield]) nmp=2;
    fs_alm2map(fs,1,spin_arr[ifield],&(maps[imap]),&(alms[imap]));//calculating the maps based on amplitude spectra calculated above
    for(imp=0;imp<nmp;imp++)
      dftw_free(alms[imap+imp]);
    imap+=nmp;
  }
  free(alms);
  nmt_flatsky_info_free(fs);

  return maps;
}
