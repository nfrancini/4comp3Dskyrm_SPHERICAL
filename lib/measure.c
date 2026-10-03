// Includo l'header
#include "../head/header.h"

// Funzione che calcola l'energia di una singola replica
double get_energy(chain_t *single_chain, param_t *myParams){
  double ene=0.0;
  double heis_term=0.0;
  double DM_term=0.0;
  // Energy with mean_field
  // Energy with PBC
  #if BC==0
    for(int iSite=0;iSite<myParams->N;iSite++){
      for(int mu=0;mu<NUM_NEAREST_NEIGHBOR;mu++){
        heis_term+=spin_product(single_chain->lattice[iSite],single_chain->lattice[nn[iSite][mu]]);
        DM_term+=spin_product(single_chain->lattice[iSite],e4)*spin_product(single_chain->lattice[nn[iSite][mu]],DM_vectors[mu])-spin_product(single_chain->lattice[iSite],DM_vectors[mu])*spin_product(single_chain->lattice[nn[iSite][mu]],e4);
      }
      ene+=(-myParams->heis)*heis_term/2.0 + (-myParams->DM)*DM_term/2.0 -spin_product(myParams->ext_field,single_chain->lattice[iSite]);
      heis_term=0;
      DM_term=0;
    }
  #elif BC==1
    for(int iSite=0;iSite<myParams->N;iSite++){
      // Bulk spin
      if((iSite%myParams->size!=0)&&((iSite/myParams->size)%myParams->size!=0)&&((iSite/((int)pow(myParams->size,2)))%myParams->size!=0)&&(iSite%myParams->size!=myParams->size-1)&&((iSite/myParams->size)%myParams->size!=myParams->size-1)&&((iSite/((int)pow(myParams->size,2)))%myParams->size!=myParams->size-1)){
        for(int mu=0;mu<NUM_NEAREST_NEIGHBOR;mu++){
          heis_term+=spin_product(single_chain->lattice[iSite],single_chain->lattice[nn[iSite][mu]]);
          DM_term+=spin_product(single_chain->lattice[iSite],e4)*spin_product(single_chain->lattice[nn[iSite][mu]],DM_vectors[mu])-spin_product(single_chain->lattice[iSite],DM_vectors[mu])*spin_product(single_chain->lattice[nn[iSite][mu]],e4);
        }
        ene+=(-myParams->heis)*heis_term/2.0 + (-myParams->DM)*DM_term/2.0 -spin_product(myParams->ext_field,single_chain->lattice[iSite]);
        heis_term=0;
        DM_term=0;
      }
      // Spin on x=0, y=/={0,L-1}, z=/={0,L-1}
      // Calculuate along +x
      else if((iSite%myParams->size==0)&&((iSite/myParams->size)%myParams->size!=0)&&((iSite/((int)pow(myParams->size,2)))%myParams->size!=0)&&((iSite/myParams->size)%myParams->size!=myParams->size-1)&&((iSite/((int)pow(myParams->size,2)))%myParams->size!=myParams->size-1)){
        heis_term+=spin_product(single_chain->lattice[iSite],single_chain->lattice[nn[iSite][mu]]);
        DM_term+=spin_product(single_chain->lattice[iSite],e4)*spin_product(single_chain->lattice[nn[iSite][mu]],DM_vectors[mu])-spin_product(single_chain->lattice[iSite],DM_vectors[mu])*spin_product(single_chain->lattice[nn[iSite][mu]],e4);
        ene+=(-myParams->heis)*heis_term/2.0 + (-myParams->DM)*DM_term/2.0;
        heis_term=0;
        DM_term=0;
      }
      // Spin on x=L-1, y=/={0,L-1}, z=/={0,L-1}
      // Calculate along -x
      else if((iSite%myParams->size==myParams->size-1)&&((iSite/myParams->size)%myParams->size!=0)&&((iSite/((int)pow(myParams->size,2)))%myParams->size!=0)&&((iSite/myParams->size)%myParams->size!=myParams->size-1)&&((iSite/((int)pow(myParams->size,2)))%myParams->size!=myParams->size-1)){
        heis_term+=spin_product(single_chain->lattice[iSite],single_chain->lattice[nn[iSite][mu]]);
        DM_term+=spin_product(single_chain->lattice[iSite],e4)*spin_product(single_chain->lattice[nn[iSite][mu]],DM_vectors[mu])-spin_product(single_chain->lattice[iSite],DM_vectors[mu])*spin_product(single_chain->lattice[nn[iSite][mu]],e4);
        ene+=(-myParams->heis)*heis_term/2.0 + (-myParams->DM)*DM_term/2.0;
        heis_term=0;
        DM_term=0;
      }
      // Spin on y=0, x=/={0,L-1}, z=/={0,L-1}
      // Calculate along +y
      else if((iSite%myParams->size!=0)&&((iSite/myParams->size)%myParams->size==0)&&((iSite/((int)pow(myParams->size,2)))%myParams->size!=0)&&(iSite%myParams->size!=myParams->size-1)&&((iSite/((int)pow(myParams->size,2)))%myParams->size!=myParams->size-1)){
        heis_term+=spin_product(single_chain->lattice[iSite],single_chain->lattice[nn[iSite][mu]]);
        DM_term+=spin_product(single_chain->lattice[iSite],e4)*spin_product(single_chain->lattice[nn[iSite][mu]],DM_vectors[mu])-spin_product(single_chain->lattice[iSite],DM_vectors[mu])*spin_product(single_chain->lattice[nn[iSite][mu]],e4);
        ene+=(-myParams->heis)*heis_term/2.0 + (-myParams->DM)*DM_term/2.0;
        heis_term=0;
        DM_term=0;
      }
      // Spin on y=L-1, x=/={0,L-1}, z=/={0,L-1}
      // Calculate along -y
      else if((iSite%myParams->size!=0)&&((iSite/myParams->size)%myParams->size==myParams->size-1)&&((iSite/((int)pow(myParams->size,2)))%myParams->size!=0)&&(iSite%myParams->size!=myParams->size-1)&&((iSite/((int)pow(myParams->size,2)))%myParams->size!=myParams->size-1)){
        heis_term+=spin_product(single_chain->lattice[iSite],single_chain->lattice[nn[iSite][mu]]);
        DM_term+=spin_product(single_chain->lattice[iSite],e4)*spin_product(single_chain->lattice[nn[iSite][mu]],DM_vectors[mu])-spin_product(single_chain->lattice[iSite],DM_vectors[mu])*spin_product(single_chain->lattice[nn[iSite][mu]],e4);
        ene+=(-myParams->heis)*heis_term/2.0 + (-myParams->DM)*DM_term/2.0;
        heis_term=0;
        DM_term=0;
      }
      // Spin on z=0, x=/={0,L-1}, y=/={0,L-1}
      // Calculate along +z
      else if((iSite%myParams->size!=0)&&((iSite/myParams->size)%myParams->size!=0)&&((iSite/((int)pow(myParams->size,2)))%myParams->size==0)&&(iSite%myParams->size!=myParams->size-1)&&((iSite/myParams->size)%myParams->size!=myParams->size-1)){
        heis_term+=spin_product(single_chain->lattice[iSite],single_chain->lattice[nn[iSite][mu]]);
        DM_term+=spin_product(single_chain->lattice[iSite],e4)*spin_product(single_chain->lattice[nn[iSite][mu]],DM_vectors[mu])-spin_product(single_chain->lattice[iSite],DM_vectors[mu])*spin_product(single_chain->lattice[nn[iSite][mu]],e4);
        ene+=(-myParams->heis)*heis_term/2.0 + (-myParams->DM)*DM_term/2.0;
        heis_term=0;
        DM_term=0;
      }
      // Spin on z=L-1, x=/={0,L-1}, y=/={0,L-1}
      // Calculate along -z
      else if((iSite%myParams->size!=0)&&((iSite/myParams->size)%myParams->size!=0)&&((iSite/((int)pow(myParams->size,2)))%myParams->size==myParams->size-1)&&(iSite%myParams->size!=myParams->size-1)&&((iSite/myParams->size)%myParams->size!=myParams->size-1)){
        heis_term+=spin_product(single_chain->lattice[iSite],single_chain->lattice[nn[iSite][mu]]);
        DM_term+=spin_product(single_chain->lattice[iSite],e4)*spin_product(single_chain->lattice[nn[iSite][mu]],DM_vectors[mu])-spin_product(single_chain->lattice[iSite],DM_vectors[mu])*spin_product(single_chain->lattice[nn[iSite][mu]],e4);
        ene+=(-myParams->heis)*heis_term/2.0 + (-myParams->DM)*DM_term/2.0;
        heis_term=0;
        DM_term=0;
      }
    }
  #endif
  return ene;
}

// Procedura che misura le varie magnetizzazioni del reticolo
// Utilizzo una supercella da 8 siti, quindi devo tenerne conto nella
// misura delle varie magnetizzazioni
void get_obs(chain_t *single_chain, param_t *myParams, obs_t *single_obs){
  double ene=0.0;
  double heis_term=0.0;
  double DM_term=0.0;
  double aux_spin[NUM_COMPS];
  double aux_topo[NUM_COMPS];
  double aux_topo_S2[NUM_COMPS];
  //Azzero gli errori se DEBUG è definito
  #ifdef DEBUG
  resetErr();
  #endif
  // Azzero tutte le misure
  single_obs->ene_density=0;
  single_obs->top_charge=0;
  single_obs->top_charge_S2=0;
  for(int i=0;i<NUM_COMPS;i++){
    single_obs->magn[i]=0;
  }
  // Measure energy, magnetization and topological charge
  // Periodic BC
  #if BC==0
    for(int i=0;i<myParams->N;i++){
      for(int mu=0;mu<NUM_NEAREST_NEIGHBOR;mu++){
        heis_term+=spin_product(single_chain->lattice[i],single_chain->lattice[nn[i][mu]]);
        DM_term+=spin_product(single_chain->lattice[i],e4)*spin_product(single_chain->lattice[nn[i][mu]],DM_vectors[mu])-spin_product(single_chain->lattice[i],DM_vectors[mu])*spin_product(single_chain->lattice[nn[i][mu]],e4);
      }
      ene+=(-myParams->heis)*heis_term/2.0 + (-myParams->DM)*DM_term/2.0 -spin_product(myParams->ext_field,single_chain->lattice[i]);
      heis_term=0;
      DM_term=0;
      for(int j=0;j<NUM_COMPS;j++){
        single_obs->magn[j] += single_chain->lattice[i][j]/myParams->N;
      }
      // Measure also topological charge
      cross_product(single_chain->lattice[nn[i][0]],single_chain->lattice[nn[i][1]],single_chain->lattice[nn[i][2]],aux_topo);
      reduced_cross_product(single_chain->lattice[nn[i][0]],single_chain->lattice[nn[i][1]],aux_topo_S2);
      single_obs->top_charge+=(spin_product(single_chain->lattice[i],aux_topo))/(2*PI*PI);
      single_obs->top_charge_S2+=(spin_product(single_chain->lattice[i],aux_topo_S2))/(4*PI);
    }
    single_obs->ene_density=ene/(myParams->N);
  // Dirichlet BC
  #elif BC==1
    for(int i=0;i<myParams->N;i++){
      // Spin in the bulk
      if((i%myParams->size!=0)&&((i/myParams->size)%myParams->size!=0)&&((i/((int)pow(myParams->size,2)))%myParams->size!=0)&&(i%myParams->size!=myParams->size-1)&&((i/myParams->size)%myParams->size!=myParams->size-1)&&((i/((int)pow(myParams->size,2)))%myParams->size!=myParams->size-1)){
        for(int mu=0;mu<NUM_NEAREST_NEIGHBOR;mu++){
          heis_term+=spin_product(single_chain->lattice[i],single_chain->lattice[nn[i][mu]]);
          DM_term+=spin_product(single_chain->lattice[i],e4)*spin_product(single_chain->lattice[nn[i][mu]],DM_vectors[mu])-spin_product(single_chain->lattice[i],DM_vectors[mu])*spin_product(single_chain->lattice[nn[i][mu]],e4);
        }
        ene+=(-myParams->heis)*heis_term/2.0 + (-myParams->DM)*DM_term/2.0 -spin_product(myParams->ext_field,single_chain->lattice[i]);
        heis_term=0;
        DM_term=0;
        // Bulk magnetization
        for(int j=0;j<NUM_COMPS;j++){
          single_obs->magn[j] += single_chain->lattice[i][j]/(pow(myParams->size-2,3));
        }
        // Measure also topological charge
        cross_product(single_chain->lattice[nn[i][0]],single_chain->lattice[nn[i][1]],single_chain->lattice[nn[i][2]],aux_topo);
        reduced_cross_product(single_chain->lattice[nn[i][0]],single_chain->lattice[nn[i][1]],aux_topo_S2);
        single_obs->top_charge+=(spin_product(single_chain->lattice[i],aux_topo))/(2*PI*PI);
        single_obs->top_charge_S2+=(spin_product(single_chain->lattice[i],aux_topo_S2))/(4*PI);
      }
      // Spin on x=0, y=/={0,L-1}, z=/={0,L-1}
      // Calculuate along +x
      else if((i%myParams->size==0)&&((i/myParams->size)%myParams->size!=0)&&((i/((int)pow(myParams->size,2)))%myParams->size!=0)&&((i/myParams->size)%myParams->size!=myParams->size-1)&&((i/((int)pow(myParams->size,2)))%myParams->size!=myParams->size-1)){
        heis_term+=spin_product(single_chain->lattice[i],single_chain->lattice[nn[i][mu]]);
        DM_term+=spin_product(single_chain->lattice[i],e4)*spin_product(single_chain->lattice[nn[i][mu]],DM_vectors[mu])-spin_product(single_chain->lattice[i],DM_vectors[mu])*spin_product(single_chain->lattice[nn[i][mu]],e4);
        ene+=(-myParams->heis)*heis_term/2.0 + (-myParams->DM)*DM_term/2.0;
        heis_term=0;
        DM_term=0;
        // Measure also topological charge
        cross_product(single_chain->lattice[nn[i][0]],single_chain->lattice[nn[i][1]],single_chain->lattice[nn[i][2]],aux_topo);
        reduced_cross_product(single_chain->lattice[nn[i][0]],single_chain->lattice[nn[i][1]],aux_topo_S2);
        single_obs->top_charge+=(spin_product(single_chain->lattice[i],aux_topo))/(2*PI*PI);
        single_obs->top_charge_S2+=(spin_product(single_chain->lattice[i],aux_topo_S2))/(4*PI);
      }
      // Spin on x=L-1, y=/={0,L-1}, z=/={0,L-1}
      // Calculate along -x
      else if((i%myParams->size==myParams->size-1)&&((i/myParams->size)%myParams->size!=0)&&((i/((int)pow(myParams->size,2)))%myParams->size!=0)&&((i/myParams->size)%myParams->size!=myParams->size-1)&&((i/((int)pow(myParams->size,2)))%myParams->size!=myParams->size-1)){
        heis_term+=spin_product(single_chain->lattice[i],single_chain->lattice[nn[i][mu]]);
        DM_term+=spin_product(single_chain->lattice[i],e4)*spin_product(single_chain->lattice[nn[i][mu]],DM_vectors[mu])-spin_product(single_chain->lattice[i],DM_vectors[mu])*spin_product(single_chain->lattice[nn[i][mu]],e4);
        ene+=(-myParams->heis)*heis_term/2.0 + (-myParams->DM)*DM_term/2.0;
        heis_term=0;
        DM_term=0;
        // Measure also topological charge
        cross_product(single_chain->lattice[nn[i][0]],single_chain->lattice[nn[i][1]],single_chain->lattice[nn[i][2]],aux_topo);
        reduced_cross_product(single_chain->lattice[nn[i][0]],single_chain->lattice[nn[i][1]],aux_topo_S2);
        single_obs->top_charge+=(spin_product(single_chain->lattice[i],aux_topo))/(2*PI*PI);
        single_obs->top_charge_S2+=(spin_product(single_chain->lattice[i],aux_topo_S2))/(4*PI);
      }
      // Spin on y=0, x=/={0,L-1}, z=/={0,L-1}
      // Calculate along +y
      else if((i%myParams->size!=0)&&((i/myParams->size)%myParams->size==0)&&((i/((int)pow(myParams->size,2)))%myParams->size!=0)&&(i%myParams->size!=myParams->size-1)&&((i/((int)pow(myParams->size,2)))%myParams->size!=myParams->size-1)){
        heis_term+=spin_product(single_chain->lattice[i],single_chain->lattice[nn[i][mu]]);
        DM_term+=spin_product(single_chain->lattice[i],e4)*spin_product(single_chain->lattice[nn[i][mu]],DM_vectors[mu])-spin_product(single_chain->lattice[i],DM_vectors[mu])*spin_product(single_chain->lattice[nn[i][mu]],e4);
        ene+=(-myParams->heis)*heis_term/2.0 + (-myParams->DM)*DM_term/2.0;
        heis_term=0;
        DM_term=0;
        // Measure also topological charge
        cross_product(single_chain->lattice[nn[i][0]],single_chain->lattice[nn[i][1]],single_chain->lattice[nn[i][2]],aux_topo);
        reduced_cross_product(single_chain->lattice[nn[i][0]],single_chain->lattice[nn[i][1]],aux_topo_S2);
        single_obs->top_charge+=(spin_product(single_chain->lattice[i],aux_topo))/(2*PI*PI);
        single_obs->top_charge_S2+=(spin_product(single_chain->lattice[i],aux_topo_S2))/(4*PI);
      }
      // Spin on y=L-1, x=/={0,L-1}, z=/={0,L-1}
      // Calculate along -y
      else if((i%myParams->size!=0)&&((i/myParams->size)%myParams->size==myParams->size-1)&&((i/((int)pow(myParams->size,2)))%myParams->size!=0)&&(i%myParams->size!=myParams->size-1)&&((i/((int)pow(myParams->size,2)))%myParams->size!=myParams->size-1)){
        heis_term+=spin_product(single_chain->lattice[i],single_chain->lattice[nn[i][mu]]);
        DM_term+=spin_product(single_chain->lattice[i],e4)*spin_product(single_chain->lattice[nn[i][mu]],DM_vectors[mu])-spin_product(single_chain->lattice[i],DM_vectors[mu])*spin_product(single_chain->lattice[nn[i][mu]],e4);
        ene+=(-myParams->heis)*heis_term/2.0 + (-myParams->DM)*DM_term/2.0;
        heis_term=0;
        DM_term=0;
        // Measure also topological charge
        cross_product(single_chain->lattice[nn[i][0]],single_chain->lattice[nn[i][1]],single_chain->lattice[nn[i][2]],aux_topo);
        reduced_cross_product(single_chain->lattice[nn[i][0]],single_chain->lattice[nn[i][1]],aux_topo_S2);
        single_obs->top_charge+=(spin_product(single_chain->lattice[i],aux_topo))/(2*PI*PI);
        single_obs->top_charge_S2+=(spin_product(single_chain->lattice[i],aux_topo_S2))/(4*PI);
      }
      // Spin on z=0, x=/={0,L-1}, y=/={0,L-1}
      // Calculate along +z
      else if((i%myParams->size!=0)&&((i/myParams->size)%myParams->size!=0)&&((i/((int)pow(myParams->size,2)))%myParams->size==0)&&(i%myParams->size!=myParams->size-1)&&((i/myParams->size)%myParams->size!=myParams->size-1)){
        heis_term+=spin_product(single_chain->lattice[i],single_chain->lattice[nn[i][mu]]);
        DM_term+=spin_product(single_chain->lattice[i],e4)*spin_product(single_chain->lattice[nn[i][mu]],DM_vectors[mu])-spin_product(single_chain->lattice[i],DM_vectors[mu])*spin_product(single_chain->lattice[nn[i][mu]],e4);
        ene+=(-myParams->heis)*heis_term/2.0 + (-myParams->DM)*DM_term/2.0;
        heis_term=0;
        DM_term=0;
        // Measure also topological charge
        cross_product(single_chain->lattice[nn[i][0]],single_chain->lattice[nn[i][1]],single_chain->lattice[nn[i][2]],aux_topo);
        reduced_cross_product(single_chain->lattice[nn[i][0]],single_chain->lattice[nn[i][1]],aux_topo_S2);
        single_obs->top_charge+=(spin_product(single_chain->lattice[i],aux_topo))/(2*PI*PI);
        single_obs->top_charge_S2+=(spin_product(single_chain->lattice[i],aux_topo_S2))/(4*PI);
      }
      // Spin on z=L-1, x=/={0,L-1}, y=/={0,L-1}
      // Calculate along -z
      else if((i%myParams->size!=0)&&((i/myParams->size)%myParams->size!=0)&&((i/((int)pow(myParams->size,2)))%myParams->size==myParams->size-1)&&(i%myParams->size!=myParams->size-1)&&((i/myParams->size)%myParams->size!=myParams->size-1)){
        heis_term+=spin_product(single_chain->lattice[i],single_chain->lattice[nn[i][mu]]);
        DM_term+=spin_product(single_chain->lattice[i],e4)*spin_product(single_chain->lattice[nn[i][mu]],DM_vectors[mu])-spin_product(single_chain->lattice[i],DM_vectors[mu])*spin_product(single_chain->lattice[nn[i][mu]],e4);
        ene+=(-myParams->heis)*heis_term/2.0 + (-myParams->DM)*DM_term/2.0;
        heis_term=0;
        DM_term=0;
        // Measure also topological charge
        cross_product(single_chain->lattice[nn[i][0]],single_chain->lattice[nn[i][1]],single_chain->lattice[nn[i][2]],aux_topo);
        reduced_cross_product(single_chain->lattice[nn[i][0]],single_chain->lattice[nn[i][1]],aux_topo_S2);
        single_obs->top_charge+=(spin_product(single_chain->lattice[i],aux_topo))/(2*PI*PI);
        single_obs->top_charge_S2+=(spin_product(single_chain->lattice[i],aux_topo_S2))/(4*PI);
      }
    }
    // Bulk energy density
    // Always divide energy by number of bulk spins
    single_obs->ene_density=ene/(pow(myParams->size-2,3));
  #endif
}

// Procedura per le misure
void get_measures(chain_t *myChains, param_t *myParams, filenames_t *myFilenames, obs_t *myObs, int *tempToreplica, int *exchAcc, dsfmt_t *dsfmt){
  FILE *fOut[myParams->n_temp];
  int rank;
  int flag;
  // Se DEBUG è definito resetto gli errori
  #ifdef DEBUG
  resetErr();
  #endif
  // Reset delle accettanze degli swap
  for(int i=0;i<myParams->n_temp-1;i++){
    exchAcc[i]=0;
  }
  // Apro i file delle misure
  for(int i=0;i<myParams->n_temp;i++){
    fOut[i]=fopen((myFilenames + i)->data_file, "a");
    if(fOut[i]==NULL){
      perror("Errore in apertura in get_measures");
      exit(EXIT_FAILURE);
    }
  }
  // Ciclo di misurazioni
  for(int j=0;j<myParams->iMis;j++){
    for(int k=0;k<myParams->iDec;k++){
      // Parallelizzo per il Metropolis
      #pragma omp parallel
      {
        #pragma omp for
        for(int i=0;i<myParams->n_temp;i++){
          rank=omp_get_thread_num();
          update_metropolis(myParams, myChains+tempToreplica[i], dsfmt, rank);
        }
        #pragma omp barrier
      }
      // Overrelaxation
      for(int j=0;j<IOVERR;j++){
        #pragma omp parallel
        {
          #pragma omp for
          for(int i=0;i<myParams->n_temp;i++){
            rank=omp_get_thread_num();
            update_overrelax(myParams, myChains+tempToreplica[i], dsfmt, rank);
          }
          #pragma omp barrier
        }
      }
      swap_temperatures(myParams, myChains, dsfmt, rank, tempToreplica, exchAcc);
    }
    // Misure in parallelo
    #pragma omp parallel
    {
      #pragma omp for
      for(int i=0;i<myParams->n_temp;i++){
        rank=omp_get_thread_num();
        // Misuro le magnetizzazioni che mi interessano
        get_obs((myChains + tempToreplica[i]), myParams, (myObs+i));
        #ifdef DEBUG
        if(errno == EDOM){
          printf("Funzione get_magn, divisione\n");
          perror("    errno == EDOM");
        }
        if(errno == ERANGE){
          printf("Funzione get_magn, divisione\n");
          perror("    errno == ERANGE");
          if(fetestexcept(FE_OVERFLOW)){
            puts("    FE_OVERFLOW was raised");
            exit(EXIT_FAILURE);
          }
          else if(fetestexcept(FE_UNDERFLOW)){
            puts("    FE_UNDERFLOW was raised");
            exit(EXIT_FAILURE);
          }
        }
        if(fetestexcept(FE_INVALID)){
          puts("    FE_INVALID was raised");
          exit(EXIT_FAILURE);
        }
        if(fetestexcept(FE_DIVBYZERO)){
          puts("    FE_DIVBYZERO was raised");
          exit(EXIT_FAILURE);
        }
        #endif
      }
      #pragma omp barrier
    }
    // Scrivo le misure
    // Apro i file delle misure
    for(int i=0;i<myParams->n_temp;i++){
      fprintf(fOut[i], "%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\n", (myObs+i)->ene_density, (myObs+i)->magn[0], (myObs+i)->magn[1], (myObs+i)->magn[2], (myObs+i)->magn[3], (myObs+i)->top_charge,(myObs+i)->top_charge_S2);
    }
    // Rinormalizzo gli spin ogni 10 000 misure
    if((j+1)%INORM==0){
      #pragma omp parallel
      {
        #pragma omp for
        for(int i=0;i<myParams->n_temp;i++){
          renormalize_spin((myChains +tempToreplica[i]), myParams);
        }
        #pragma omp barrier
      }
    }
    // Salvo le configurazioni a T finita
    if((myParams->iSave!=0)&&((j+1)%myParams->iSave==0)){
      save_Configs(myChains, myParams, myFilenames, tempToreplica);
    }
    // Backup
    if((j+1)%myParams->iBackup==0){
      // Rimuovo il file contentente i nomi e lo riscrivo
      flag=remove(myParams->filenames_file);
      if(flag!=0){
        perror("Errore nel rimuovere filenames in get_measure\n");
        exit(EXIT_FAILURE);
      }
      write_Status(myParams, myChains, myFilenames, tempToreplica);
    }
  }
  // Chiudo i file delle misure
  for(int i=0;i<myParams->n_temp;i++){
    fclose(fOut[i]);
  }
}
