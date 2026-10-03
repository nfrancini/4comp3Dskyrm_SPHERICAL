// Includo l'header
#include "../head/header.h"

// PROCEDURA GEOMETRY PER MEMORIZZARE I PRIMI VICINI;
void geometry(param_t *myParams){
  // SFRUTTO L'ORDINAMENTO LESSICOGRAFICO i_x = x1 + x2*L + x3*L^2
  // IN QUESTO MODO
  // x1 = i_x % L
  // x2 = (i_x / L) % L
  // x3 = (i_x / L^2) % L

  // IMPLEMENTAZIONE NPP
  for(int i=0;i<(myParams->N);i++){
    for(int mu=0;mu<D;mu++){
      if((i / (int)pow(myParams->size,mu)) % (myParams->size) == ((myParams->size)-1)){
        nn[i][mu] = i - (int)pow(myParams->size,mu)*((myParams->size)-1);
      }
      else nn[i][mu] = i + (int)pow(myParams->size,mu);
    }
  }

  // IMPLEMENTAZIONE NMM
  for(int i=0;i<(myParams->N);i++){
    for(int mu=0;mu<D;mu++){
      if((i / (int)pow((myParams->size),mu)) % (myParams->size) == 0){
        nn[i][mu+D] = i + (int)pow((myParams->size),mu)*((myParams->size)-1);
      }
      else nn[i][mu+D] = i - (int)pow((myParams->size),mu);
    }
  }
}

// Implementazione inizializzazione array ausiliari
void initialize_auxiliary_array(param_t *myParams, int *tempToreplica, int *exchAcc){
  for (int i = 0; i < myParams->n_temp; i++) {
    tempToreplica[i]=i;
  }
  for (int i = 0; i < myParams->n_temp-1; i++) {
    exchAcc[i]=0;
  }
}

// Initialize DM vectors
// These are the delta vectors
// 6 vectors: DM[0]=(1,0,0)
//            DM[1]=(0,1,0)
//            DM[2]=(0,0,1)
//            DM[3]=(-1,0,0)
//            DM[4]=(0,-1,0)
//            DM[5]=(0,0,-1)
// Here i initialize also e4
void initialize_DM_vectors(){
  for(int nu=0;nu<NUM_NEAREST_NEIGHBOR;nu++){
    for(int mu=0;mu<NUM_COMPS;mu++){
      DM_vectors[nu][mu]=0.0;
      e4[mu]=0.0;
    }
  }
  DM_vectors[0][0]=1.0;
  DM_vectors[D][0]=-1.0;

  DM_vectors[1][1]=1.0;
  DM_vectors[D+1][1]=-1.0;

  DM_vectors[2][2]=1.0;
  DM_vectors[D+2][2]=-1.0;

  e4[3]=1.0;
}

// Libreria per implementare le funzioni di inizalizzazione
// Le funzioni di gestione input/output sono trattate in un altro file
dsfmt_t * initialize_system(chain_t **myChains, param_t *myParams, obs_t **myObs, filenames_t **myFilenames, int **tempToreplica, int **exchAcc, dsfmt_t *dsfmt, char const *finput, chain_t *zero_T_replica){
  // Se DEBUG è definito allora resetto gli errori
  #ifdef DEBUG
  resetErr();
  #endif
  // 1. Leggo il file di input
  read_from_input_Param(myParams, finput);
  // 2. Alloco la memoria
  nn = malloc(sizeof(int[myParams->N][NUM_NEAREST_NEIGHBOR]));
  zero_T_replica->lattice = malloc(sizeof(double[myParams->N][NUM_COMPS]));
  *tempToreplica = (int *)malloc((myParams->n_temp)*sizeof(int));
  *exchAcc = (int *)malloc((myParams->n_temp-1)*sizeof(int));
  *myObs = (obs_t *)malloc((myParams->n_temp)*sizeof(obs_t));
  *myChains = (chain_t *)malloc((myParams->n_temp)*sizeof(chain_t));
  *myFilenames = (filenames_t *)malloc((myParams->n_temp)*sizeof(filenames_t));
  for(int i=0;i<myParams->n_temp;i++){
    (*myChains +i)->lattice = malloc(sizeof(double[myParams->N][NUM_COMPS]));
    (*myObs +i)->ene_density=0.0;
    (*myObs +i)->top_charge=0.0;
    (*myObs +i)->top_charge_S2=0.0;
    for(int j=0;j<NUM_COMPS;j++){
      (*myObs +i)->magn[j]=0.0;
    }
  }
  initialize_auxiliary_array(myParams, *tempToreplica, *exchAcc);
  // 3. PBC
  geometry(myParams);
  // 4. Generatore random
  int rank_max = omp_get_max_threads();                     // Numero massimo di threads
  dsfmt = (dsfmt_t *)malloc(rank_max * sizeof(dsfmt_t));    // Alloco dinamicamente la memoria
  dsfmt_init_gen_rand(&dsfmt[0], myParams->dSFMT_seed);     // Inizializzo a partire dall'indice di thread 0
  for(int i=1;i<rank_max;i++){
    dsfmt[i]=dsfmt[i-1];
    dSFMT_jump(&dsfmt[i], jump);
  }
  // 5. Inizializzo i vettori DM
  initialize_DM_vectors();
  // 6. Inizializzo a seconda di iStart
  if(myParams->iStart==2){
    // Leggo da file precedente
    read_from_filenames(*myFilenames, myParams);
    read_Status(myParams, *myChains, *myFilenames, *tempToreplica);
  }
  else if(myParams->iStart==1){           // Partenza a freddo (ordinata) -> ALONG VACUUM
    for(int i=0;i<myParams->n_temp;i++){
      // Azzero acc, count e definisco le temperature di partenza equispaziate
      (*myChains + i)->acc=0;
      (*myChains + i)->count=0;
      (*myChains + i)->T=myParams->t_max * pow(myParams->t_min/myParams->t_max, (double)i/(myParams->n_temp-1));
      (*myChains + i)->acc_Parameter = 0.5;
      for(int j=0;j<myParams->N;j++){
        for(int k=0;k<NUM_COMPS;k++){
          // VACUUM ORDER
          if(k==NUM_COMPS-1){
            (*myChains + i)->lattice[j][k]=1.0;
          }
          else{
            (*myChains + i)->lattice[j][k]=0.0;
          }
        }
      }
    }
  }
  else if(myParams->iStart==0){           // Partenza a caldo (disordinata)
    for(int i=0;i<myParams->n_temp;i++){
      // Azzero acc, count e definisco le temperature di partenza equispaziate
      (*myChains + i)->acc=0;
      (*myChains + i)->count=0;
      (*myChains + i)->T=myParams->t_max * pow(myParams->t_min/myParams->t_max, (double)i/(myParams->n_temp-1));
      (*myChains + i)->acc_Parameter = PI/4;
      for(int j=0;j<myParams->N;j++){
        // PBC
        #if BC==0
          for(int k=0;k<NUM_COMPS;k++){     // Componenti random tra [-1;1) che dopo normalizzo
            (*myChains + i)->lattice[j][k]=2*dsfmt_genrand_close_open(&dsfmt[omp_get_num_threads()])-1.0;
          }
        #endif
        // FIXED BC
        #if BC==1
          // If x=0 or y=0 or z=0 then fix spin to vacuum
          if((j%myParams->size==0)||((j/myParams->size)%myParams->size==0)||((j/((int)pow(myParams->size,2)))%myParams->size==0)){
            for(int k=0;k<NUM_COMPS;k++){
              // Vacuum order
              if(k==NUM_COMPS-1){
                (*myChains + i)->lattice[j][k]=1.0;
              }
              else{
                (*myChains + i)->lattice[j][k]=0.0;
              }
            }
          }
          // If x=L-1 or y=L-1 or z=L-1 then fix spin to vacuum
          else if((j%myParams->size==myParams->size-1)||((j/myParams->size)%myParams->size==myParams->size-1)||((j/((int)pow(myParams->size,2)))%myParams->size==myParams->size-1)){
            for(int k=0;k<NUM_COMPS;k++){
              // Vacuum order
              if(k==NUM_COMPS-1){
                (*myChains + i)->lattice[j][k]=1.0;
              }
              else{
                (*myChains + i)->lattice[j][k]=0.0;
              }
            }
          }
          else{   // Random configuration for bulk spins
            for(int k=0;k<NUM_COMPS;k++){     // Componenti random tra [-1;1) che dopo normalizzo
              (*myChains + i)->lattice[j][k]=2*dsfmt_genrand_close_open(&dsfmt[omp_get_num_threads()])-1.0;
            }
          }
        #endif
        // APBC
      }
    }
    #pragma omp parallel
    {
      #pragma omp for
      for(int i=0;i<myParams->n_temp;i++){  // Rinormalizzo gli array random
        renormalize_spin((*myChains +i), myParams);
      }
      #pragma omp barrier
    }
  }
  return dsfmt;
}

// Procedura che controlla le accettanze Metropolis
void check_metro_acc(chain_t *single_chain, param_t *myParams){
  double acc;
  double step = PI/100;
  acc=(double)single_chain->acc/single_chain->count;
  if(acc>acc_max){
    single_chain->acc_Parameter=single_chain->acc_Parameter + step;
    if(single_chain->acc_Parameter>PI/2){
      single_chain->acc_Parameter=PI/2-PI/100;
    }
  }
  else if(acc<acc_min){
    single_chain->acc_Parameter=single_chain->acc_Parameter - step;
    if(single_chain->acc_Parameter<PI/100){
      single_chain->acc_Parameter=PI/50;
    }
  }
  single_chain->acc=0;
  single_chain->count=0;
}

// Procedura che aggiusta le temperature
// In questa procedura aggiusto anche i parametri di accettanza
void set_temperatures(chain_t *myChains, param_t *myParams, dsfmt_t *dsfmt, int *tempToreplica, int *exchAcc){
  // Array necessari alla procedura
  double g[myParams->n_temp], auxEnergies[myParams->n_temp];
  double beta_prime;
  int rank;
  // L'obietttivo finale è avere più o meno
  // un'accettanza piatta su tutto il sistema
  for(int t=0;t<myParams->iSet;t++){
    // Azzero il contatore degli swap
    for(int i=0;i<myParams->n_temp-1;i++){
      exchAcc[i]=0;
    }
    // Resetto gli array ausiliari
    for(int i=0;i<myParams->n_temp;i++){
      g[i]=0;
      auxEnergies[i]=0;
    }
    // Termalizzazione
    // Prima monte carlo indipendenti in parallelo
    // Poi swap
    for(int j=0;j<myParams->iPT;j++){
      #pragma omp parallel
      {
        #pragma omp for
        for(int i=0;i<myParams->n_temp;i++){
          rank=omp_get_thread_num();
          update_metropolis(myParams, (myChains + tempToreplica[i]), dsfmt, rank);
        }
        #pragma omp barrier
      }
      #pragma omp parallel
      {
        #pragma omp for
        for(int i=0;i<myParams->n_temp;i++){
          check_metro_acc((myChains + tempToreplica[i]), myParams);
        }
        #pragma omp barrier
      }
      // Scambio le temperature
      swap_temperatures(myParams, myChains, dsfmt, rank, tempToreplica, exchAcc);
      // Calcolo le energie medie
      // In parallelo
      #pragma omp parallel
      {
        #pragma omp for
        for(int i=0;i<myParams->n_temp;i++){
          auxEnergies[i]+=get_energy((myChains + tempToreplica[i]), myParams)/myParams->iPT;
        }
        #pragma omp barrier
      }
    }
    // Costruisco il vettore g
    for(int i=1;i<myParams->n_temp;i++){
      if(auxEnergies[i-1]==auxEnergies[i+1]){
        // Energie uguali
        g[i] = 1.0/(myChains+tempToreplica[i])->T;
      }
      else{
        g[i] = (auxEnergies[i-1]/(myChains+tempToreplica[i-1])->T - auxEnergies[i+1]/(myChains+tempToreplica[i+1])->T - auxEnergies[i]*(1.0/(myChains+tempToreplica[i-1])->T - 1.0/(myChains+tempToreplica[i+1])->T)) /(auxEnergies[i-1]-auxEnergies[i+1]);
      }
    }
    if((t%2)==0){
      // Cambio le temperature con indice pari tenendo le dispari fisse
      for(int i=2;i<myParams->n_temp-1;i=i+2){
        beta_prime = (1.0/((myChains+tempToreplica[i])->T) + g[i])/2.0;
        if(((1.0/beta_prime)<myParams->t_min)||((1.0/beta_prime)>myParams->t_max)){
          // DOES NOTHING
        }
        else{
          myChains[tempToreplica[i]].T = 1.0/beta_prime;
        }
      }
    }
    else if((t%2)==1){
      for(int i=1;i<myParams->n_temp-1;i=i+2){
        beta_prime = (1.0/((myChains + tempToreplica[i])->T) + g[i])/2.0;
        if(((1.0/beta_prime)<myParams->t_min)||((1.0/beta_prime)>myParams->t_max)){
          // DOES NOTHING
        }
        else{
          myChains[tempToreplica[i]].T = 1.0/beta_prime;
        }
      }
    }
    else{
      printf("Errore in set_temperatures\n");
      exit(EXIT_FAILURE);
    }
  }
  // CHEK IF EVERY TEMPERATURE IS IN GIVEN RANGE
  for(int i=1;i<myParams->n_temp-1;i++){
    if(((myChains+tempToreplica[i])->T>myParams->t_max)||((myChains+tempToreplica[i])->T<myParams->t_min)){
      printf("Errore in set_temperatures Tmin=%lf, Tmax=%lf , T=%lf\n", myParams->t_min, myParams->t_max, (myChains+tempToreplica[i])->T);
      exit(EXIT_FAILURE);
    }
  }
}

// Procedura per le accettanze
// Da chiamare una volta fissate le temperature
void set_acceptances(chain_t *myChains, param_t *myParams, dsfmt_t *dsfmt, int *tempToreplica){
  int rank;
  for(int j=0;j<myParams->iTerm;j++){
    #pragma omp parallel
    {
      #pragma omp for
      for(int i=0;i<myParams->n_temp;i++){
        rank=omp_get_thread_num();
        update_metropolis(myParams, (myChains + tempToreplica[i]), dsfmt, rank);
      }
      #pragma omp barrier
    }
    #pragma omp parallel
    {
      #pragma omp for
      for(int k=0;k<myParams->n_temp;k++){
        check_metro_acc((myChains + tempToreplica[k]), myParams);
      }
      #pragma omp barrier
    }
  }
}

// Procedura che crea i nomi dei file su cui scrivo
void set_filenames(chain_t *myChains, param_t *myParams, filenames_t *myFilenames, int *tempToreplica){
  for(int i=0;i<myParams->n_temp;i++){
    // Creo le stringhe che per i nomi dei file
    sprintf((myFilenames+i)->data_file, "data/data_%d_%lf.dat", myParams->size, (myChains+tempToreplica[i])->T);
    sprintf((myFilenames+i)->backup_file, "backup/backup_%d_%lf.bin", myParams->size, (myChains+tempToreplica[i])->T);
    sprintf((myFilenames+i)->aux_file, "aux/aux_%d_%lf.bin", myParams->size, (myChains+tempToreplica[i])->T);
    if(myParams->iSave!=0){
      sprintf((myFilenames+i)->conf_file, "conf/conf_%d_%lf.txt", myParams->size, (myChains+tempToreplica[i])->T);
    }
  }
  // Scrvo sui file dei dati e delle configurazioni l'header
  for(int i=0;i<myParams->n_temp;i++){
    FILE *fptr1;
    FILE *fptr2;
    fptr1=fopen((myFilenames+i)->data_file, "w");
    if(fptr1==NULL){
      perror("Errore in apertura dei file in set_filenames");
      exit(EXIT_FAILURE);
    }
    else{
      fprintf(fptr1, "%d\t%d\t%d\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\n", myParams->size, myParams->N, NUM_COMPS, (myChains+tempToreplica[i])->T, myParams->heis, myParams->DM, myParams->ext_field[0], myParams->ext_field[1], myParams->ext_field[2], myParams->ext_field[3]);
      fprintf(fptr1, "Ene_dens\t\tm1\tm2\tm3\tm4\tTop_Charge\tCharge_S2\n\n");
    }
    fclose(fptr1);
    if(myParams->iSave!=0){
      fptr2=fopen((myFilenames+i)->conf_file, "w");
      if(fptr2==NULL){
        perror("Errore in apertura dei file in set_filenames");
        exit(EXIT_FAILURE);
      }
      else {
        fprintf(fptr2, "%d\t%d\t%d\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\t%.13lf\n", myParams->size, myParams->N, NUM_COMPS, (myChains+tempToreplica[i])->T, myParams->heis, myParams->DM, myParams->ext_field[0], myParams->ext_field[1], myParams->ext_field[2], myParams->ext_field[3]);
        fprintf(fptr2, "iSite\tS1\tS2\tS3\tS4\n\n");
      }
      fclose(fptr2);
    }
  }
}

// Procedura per la prima termalizzazione
void burn_in(chain_t *myChains, param_t *myParams, dsfmt_t *dsfmt, int *tempToreplica, int *exchAcc){
  int rank;
  // Resetto le accettanze
  for(int i=0;i<myParams->n_temp-1;i++){
    exchAcc[i]=0;
  }
  // Cicli a vuoto
  for(int j=0;j<myParams->iTerm;j++){
    #pragma omp parallel
    {
      #pragma omp for
      for(int i=0;i<myParams->n_temp;i++){
        rank=omp_get_thread_num();
        update_metropolis(myParams, (myChains + tempToreplica[i]), dsfmt, rank);
      }
      #pragma omp barrier
    }
    swap_temperatures(myParams, myChains, dsfmt, rank, tempToreplica, exchAcc);
  }
}

// Procedura per la deallocazione della memoria
void deallocation(chain_t **myChains, param_t *myParams, obs_t **myObs, filenames_t **myFilenames, int **tempToreplica, int **exchAcc, dsfmt_t *dsfmt, chain_t *zero_T_replica){
  free(nn);
  free(*tempToreplica);
  free(*myObs);
  free(*exchAcc);
  free(dsfmt);
  free(*myFilenames);
  for(int i=0;i<myParams->n_temp;i++){
    free((*myChains +i)->lattice);
  }
  free(*myChains);
  free(zero_T_replica->lattice);
}
