// Includo l'header
#include "../head/header.h"

// Procedura per creare lo spin di prova
// La procedura prende in ingresso una singola catena
// il sito in cui viene creato lo spin di prova
// ed ovviamente i parametri di sistema ed un vettore trial
// dove salvo lo spin di prova.
// L'idea è di usare gli angoli di Eulero:
// Lo spin scelto può essere visto come combinazione di due
// rotazioni del vettore (001). Quindi prima genero un
// cono di punti uniformi attorno ad (001), poi lo ruoto sfruttando
// le rotazioni del vettore scelto. Così ottengo un cono
// di punti casuali attorno al punto scelto.
void spin_trial(chain_t *single_chain, dsfmt_t *dsfmt, int rank, int iSite, double *trial){
  // Metodo Angoli di Eulero
  double theta;
  int i,j;
  // Se DEBUG è definito resetto gli errori
  #ifdef DEBUG
  resetErr();
  #endif

  // Copio lo spin nel vettore trial
  spin_copy(single_chain->lattice[iSite],trial);
  // SELEZIONO DUE COMPONENTI RANDOM
  i = floor((NUM_COMPS)*dsfmt_genrand_close_open(&dsfmt[rank]));         // SCELGO DUE COMPONENTI DA MESCOLARE IN MODO RANDOM
  j = floor((NUM_COMPS)*dsfmt_genrand_close_open(&dsfmt[rank]));
  while (j == i){                // PASSAGGIO NECESSARIO PER IL FUNZIONAMENTO (DEVO AVERE i != j)
    j = floor((NUM_COMPS)*dsfmt_genrand_close_open(&dsfmt[rank]));
  }
  // ANGOLO UNIFORMEMENTE SCELTO TRA [0 E XI_MAX]
  theta=single_chain->acc_Parameter*(2*dsfmt_genrand_close_open(&dsfmt[rank])-1.0);
  // RUOTO SOLO DUE COMPONENTI
  trial[i]=cos(theta)*single_chain->lattice[iSite][i] - sin(theta)*single_chain->lattice[iSite][j];
  trial[j]=sin(theta)*single_chain->lattice[iSite][i] + cos(theta)*single_chain->lattice[iSite][j];

  #ifdef DEBUG
  if(errno == EDOM){
    printf("Funzione spin_trial, trial spin: cos or sin\n");
    perror("    errno == EDOM");
  }
  if(errno == ERANGE){
    printf("Funzione spin_trial, trial spin: cos or sin \n");
    perror("    errno == ERANGE");
  }
  if(fetestexcept(FE_INVALID)){
    puts("    FE_INVALID was raised");
    exit(EXIT_FAILURE);
  }
  #endif
}

// Procedura per creare il campo medio
// Prima calcolo la parte relativa al modello HKGG'
void mean_field(chain_t *single_chain, param_t *myParams, int iSite, double *mean_field){
  // Inizializzo a zero
  for(int i=0;i<NUM_COMPS;i++){
    mean_field[i]=-myParams->ext_field[i];
  }
  for(int j=0;j<NUM_NEAREST_NEIGHBOR;j++){
    for(int i=0;i<NUM_COMPS;i++){
      mean_field[i]+=(-myParams->heis)*(single_chain->lattice[nn[iSite][j]][i]) + (-myParams->DM)*(spin_product(single_chain->lattice[nn[iSite][j]],DM_vectors[j])*e4[i] - spin_product(single_chain->lattice[nn[iSite][j]],e4)*DM_vectors[j][i]);
    }
  }
}

// Update Metropolis
void update_metropolis(param_t *myParams, chain_t *single_chain, dsfmt_t *dsfmt, int rank){
  double trial[NUM_COMPS], diff[NUM_COMPS], mean_f[NUM_COMPS];
  double c=0;
  // Se DEBUG è definito resetto gli errori
  #ifdef DEBUG
  resetErr();
  #endif
  #if BC==0
  // Ciclo sui tutti i siti
  for(int iSite=0;iSite<myParams->N;iSite++){
    // Aumento il contatore
    single_chain->count+=1;
    // Creo lo spin di prova
    spin_trial(single_chain, dsfmt, rank, iSite, trial);
    // Differenza tra spin di prova e spin iniziale nel sito iSite
    spin_difference(trial, single_chain->lattice[iSite], diff);
    // Calcolo il campo medio nel sito iSite
    mean_field(single_chain, myParams, iSite, mean_f);
    // Differenza di energia
    c=-spin_product(diff, mean_f)/single_chain->T;
    #ifdef DEBUG
    double ene1;
    ene1=get_energy(single_chain, myParams);
    double dE=spin_product(diff, mean_f);
    #endif
    // Test dell'accettanza
    if((c>=0) || (exp(c)>dsfmt_genrand_close_open(&dsfmt[rank]))){
      #ifdef DEBUG
      if(errno == ERANGE){
        printf("Funzione update_metro\n");
        perror("    errno == ERANGE");
      }
      if(fetestexcept(FE_OVERFLOW)){
        printf("    FE_OVERFLOW was raised\n");
        exit(EXIT_FAILURE);
      }
      else if(fetestexcept(FE_UNDERFLOW)){
        printf("    FE_UNDERFLOW was raised\n");
        exit(EXIT_FAILURE);
      }
      else if(fetestexcept(FE_DIVBYZERO)){
        printf("    FE_DIVBYZERO was raised\n");
        exit(EXIT_FAILURE);
      }
      #endif
      // Accetto il passo e copio lo spin di prova
      spin_copy(trial, single_chain->lattice[iSite]);
      single_chain->acc+=1;
      #ifdef DEBUG
      double ene2;
      ene2=get_energy(single_chain, myParams);
      double delta_glob=ene2-ene1;
      if(fabs(delta_glob-dE)>1.0e-12){
        printf("Errore, ene_locale %.13lf\tene_globale %.13lf\n", dE, delta_glob);
        exit(EXIT_FAILURE);
      }
      #endif
    }
  }
  #elif BC==1
  // Ciclo sui tutti i siti
  for(int iSite=0;iSite<myParams->N;iSite++){
    if((iSite%myParams->size!=0)&&((iSite/myParams->size)%myParams->size!=0)&&((iSite/((int)pow(myParams->size,2)))%myParams->size!=0)&&(iSite%myParams->size!=myParams->size-1)&&((iSite/myParams->size)%myParams->size!=myParams->size-1)&&((iSite/((int)pow(myParams->size,2)))%myParams->size!=myParams->size-1)){
      // Aumento il contatore
      single_chain->count+=1;
      // Creo lo spin di prova
      spin_trial(single_chain, dsfmt, rank, iSite, trial);
      // Differenza tra spin di prova e spin iniziale nel sito iSite
      spin_difference(trial, single_chain->lattice[iSite], diff);
      // Calcolo il campo medio nel sito iSite
      mean_field(single_chain, myParams, iSite, mean_f);
      // Differenza di energia
      c=-spin_product(diff, mean_f)/single_chain->T;
      #ifdef DEBUG
      double ene1;
      ene1=get_energy(single_chain, myParams);
      double dE=spin_product(diff, mean_f);
      #endif
      // Test dell'accettanza
      if((c>=0) || (exp(c)>dsfmt_genrand_close_open(&dsfmt[rank]))){
        #ifdef DEBUG
        if(errno == ERANGE){
          printf("Funzione update_metro\n");
          perror("    errno == ERANGE");
        }
        if(fetestexcept(FE_OVERFLOW)){
          printf("    FE_OVERFLOW was raised\n");
          exit(EXIT_FAILURE);
        }
        else if(fetestexcept(FE_UNDERFLOW)){
          printf("    FE_UNDERFLOW was raised\n");
          exit(EXIT_FAILURE);
        }
        else if(fetestexcept(FE_DIVBYZERO)){
          printf("    FE_DIVBYZERO was raised\n");
          exit(EXIT_FAILURE);
        }
        #endif
        // Accetto il passo e copio lo spin di prova
        spin_copy(trial, single_chain->lattice[iSite]);
        single_chain->acc+=1;
        #ifdef DEBUG
        double ene2;
        ene2=get_energy(single_chain, myParams);
        double delta_glob=ene2-ene1;
        if(fabs(delta_glob-dE)>1.0e-12){
          printf("Errore, ene_locale %.13lf\tene_globale %.13lf\n", dE, delta_glob);
          exit(EXIT_FAILURE);
        }
        #endif
      }
    }
  }
  #endif
}

// Swap delle temperature
void swap_temperatures(param_t *myParams, chain_t *myChains, dsfmt_t *dsfmt, int rank, int *tempToreplica, int *exchAcc){
  double delta;
  double tempT;
  double tempAccparam;
  int tempindex;
  // Se DEBUG è definito resetto gli errori
  #ifdef DEBUG
  resetErr();
  #endif
  // Ciclo sulle temperature, m è un indice di temperatura
  for(int m=0; m<myParams->n_temp-1;m++){
    // Condizione di swapping tra T_m e T_{m+1}
    delta = (1.0/myChains[tempToreplica[m+1]].T - 1.0/myChains[tempToreplica[m]].T) * (get_energy((myChains + tempToreplica[m]), myParams) - get_energy((myChains + tempToreplica[m+1]), myParams));
    if((delta<0) || (exp(-delta)>dsfmt_genrand_close_open(&dsfmt[rank]))){
      #ifdef DEBUG
      if(errno == ERANGE){
        printf("Funzione swap_temperatures\n");
        perror("    errno == ERANGE");
      }
      if(fetestexcept(FE_OVERFLOW)){
        printf("    FE_OVERFLOW was raised\n");
        exit(EXIT_FAILURE);
      }
      else if(fetestexcept(FE_UNDERFLOW)){
        printf("    FE_UNDERFLOW was raised\n");
        exit(EXIT_FAILURE);
      }
      else if(fetestexcept(FE_DIVBYZERO)){
        printf("    FE_DIVBYZERO was raised\n");
        exit(EXIT_FAILURE);
      }
      #endif
      // Swap avvenuto
      exchAcc[m]+=1;
      // Swap delle temperature
      tempT=myChains[tempToreplica[m]].T;
      myChains[tempToreplica[m]].T=myChains[tempToreplica[m+1]].T;
      myChains[tempToreplica[m+1]].T=tempT;
      // Swap del parametro di accettanza
      tempAccparam=myChains[tempToreplica[m]].acc_Parameter;
      myChains[tempToreplica[m]].acc_Parameter=myChains[tempToreplica[m+1]].acc_Parameter;
      myChains[tempToreplica[m+1]].acc_Parameter = tempAccparam;
      // Tracking dello swap per tempToreplica
      tempindex=tempToreplica[m];
      tempToreplica[m]=tempToreplica[m+1];
      tempToreplica[m+1]=tempindex;
    }
    // Swap rifiutato
  }
}

// Overrelaxation
void update_overrelax(param_t *myParams, chain_t *single_chain, dsfmt_t *dsfmt, int rank){
  int iSite, j;
  double mean_f[NUM_COMPS], trial[NUM_COMPS];
  double norm;

  #if BC==0
  for(iSite=0; iSite<myParams->N; iSite++){
    #ifdef DEBUG
    double ene1;
    ene1=get_energy(single_chain, myParams);
    #endif
    // Calculate mean field
    mean_field(single_chain, myParams, iSite, mean_f);
    norm=spin_product(mean_f, mean_f);
    // If the norm is too small
    if(norm<1.0e-12){
      for(j=0;j<NUM_COMPS;j++){
        trial[j] = 2*dsfmt_genrand_close_open(&dsfmt[rank])-1.0;
      }
      for(j=0;j<NUM_COMPS;j++){
        trial[j] = trial[j]/sqrt(spin_product(trial, trial));
      }
      spin_copy(trial, single_chain->lattice[iSite]);
      return;
    }
    else{
      for(j=0;j<NUM_COMPS;j++){
        trial[j]=2.0*spin_product(single_chain->lattice[iSite], mean_f)*mean_f[j]/norm - single_chain->lattice[iSite][j];
      }
      spin_copy(trial, single_chain->lattice[iSite]);
    }
  }
  #elif BC==1
  for(iSite=0; iSite<myParams->N; iSite++){
    if((iSite%myParams->size!=0)&&((iSite/myParams->size)%myParams->size!=0)&&((iSite/((int)pow(myParams->size,2)))%myParams->size!=0)&&(iSite%myParams->size!=myParams->size-1)&&((iSite/myParams->size)%myParams->size!=myParams->size-1)&&((iSite/((int)pow(myParams->size,2)))%myParams->size!=myParams->size-1)){
      #ifdef DEBUG
      double ene1;
      ene1=get_energy(single_chain, myParams);
      #endif
      // Calculate mean field
      mean_field(single_chain, myParams, iSite, mean_f);
      norm=spin_product(mean_f, mean_f);
      // If the norm is too small
      if(norm<1.0e-12){
        for(j=0;j<NUM_COMPS;j++){
          trial[j] = 2*dsfmt_genrand_close_open(&dsfmt[rank])-1.0;
        }
        for(j=0;j<NUM_COMPS;j++){
          trial[j] = trial[j]/sqrt(spin_product(trial, trial));
        }
        spin_copy(trial, single_chain->lattice[iSite]);
        return;
      }
      else{
        for(j=0;j<NUM_COMPS;j++){
          trial[j]=2.0*spin_product(single_chain->lattice[iSite], mean_f)*mean_f[j]/norm - single_chain->lattice[iSite][j];
        }
        spin_copy(trial, single_chain->lattice[iSite]);
      }
    }
  }
  #endif
}

// Minimizzazione iterativa
void iterative_minimization(param_t *myParams, chain_t *single_chain, dsfmt_t *dsfmt, int rank){
  double mean_f[NUM_COMPS];
  double trial[NUM_COMPS];
  double diff[NUM_COMPS];
  double norm=0;
  double temp;
  double max_diff=0;
  int exit=0;
  unsigned long iter=0;
  // Se DEBUG è definito resetto gli errori
  #ifdef DEBUG
  resetErr();
  #endif
  while(exit==0){
    max_diff=0;
    iter+=1;
    #if BC==0
    // Spazzata sul reticolo
    for(int iSite=0;iSite<myParams->N;iSite++){
      norm=0;
      // Campo medio in iSite
      mean_field(single_chain, myParams, iSite, mean_f);
      // Calcolo la norma del campo medio e lo rinormalizzo
      for(int j=0;j<NUM_COMPS;j++){
        norm+=pow(mean_f[j],2);
      }
      norm=sqrt(norm);
      // If norm is too small random trial spin
      if(norm<1.0e-12){
        for(int j=0;j<NUM_COMPS;j++){
          trial[j] = 2*dsfmt_genrand_close_open(&dsfmt[rank])-1.0;
        }
        for(int j=0;j<NUM_COMPS;j++){
          trial[j] = trial[j]/sqrt(spin_product(trial, trial));
        }
      }
      else{
        // Creo lo spin di prova
        for(int j=0;j<NUM_COMPS;j++){
          trial[j]=-mean_f[j]/norm;
        }
      }
      // Calcolo la differenza tra spin di prova e spin originale
      spin_difference(single_chain->lattice[iSite], trial, diff);
      temp=sqrt(spin_product(diff,diff));
      if(temp>max_diff){
        max_diff=temp;
      }
      // Copio lo spin di prova nel sito iSite
      spin_copy(trial, single_chain->lattice[iSite]);
    }
    if(fabs(max_diff)<=1.0e-12){
      exit=1;
    }
    #elif BC==1
    // Spazzata sul reticolo con fixed bc
    for(int iSite=0;iSite<myParams->N;iSite++){
      if((iSite%myParams->size!=0)&&((iSite/myParams->size)%myParams->size!=0)&&((iSite/((int)pow(myParams->size,2)))%myParams->size!=0)&&(iSite%myParams->size!=myParams->size-1)&&((iSite/myParams->size)%myParams->size!=myParams->size-1)&&((iSite/((int)pow(myParams->size,2)))%myParams->size!=myParams->size-1)){
        norm=0;
        // Campo medio in iSite
        mean_field(single_chain, myParams, iSite, mean_f);
        // Calcolo la norma del campo medio e lo rinormalizzo
        for(int j=0;j<NUM_COMPS;j++){
          norm+=pow(mean_f[j],2);
        }
        norm=sqrt(norm);
        // If norm is too small random trial spin
        if(norm<1.0e-12){
          for(int j=0;j<NUM_COMPS;j++){
            trial[j] = 2*dsfmt_genrand_close_open(&dsfmt[rank])-1.0;
          }
          for(int j=0;j<NUM_COMPS;j++){
            trial[j] = trial[j]/sqrt(spin_product(trial, trial));
          }
        }
        else{
          // Creo lo spin di prova
          for(int j=0;j<NUM_COMPS;j++){
            trial[j]=-mean_f[j]/norm;
          }
        }
        // Calcolo la differenza tra spin di prova e spin originale
        spin_difference(single_chain->lattice[iSite], trial, diff);
        temp=sqrt(spin_product(diff,diff));
        if(temp>max_diff){
          max_diff=temp;
        }
        // Copio lo spin di prova nel sito iSite
        spin_copy(trial, single_chain->lattice[iSite]);
      }
    }
    if(fabs(max_diff)<=1.0e-12){
      exit=1;
    }
    #endif
  }
}

// Procedura di minimizzazione
void push_to_zero(param_t *myParams, chain_t *low_T_chain, chain_t *zero_T_replica, dsfmt_t *dsfmt, int rank){
  double min_energy=0;
  double temp;
  FILE *fptr;
  char *file_name[STD_STRING_LENGTH];
  // Creo il nome del file per le misurazioni
  for(int i=0;i<myParams->iTerm;i++){
    // Prima copio la replica a temperatura zero
    replica_copy(myParams, low_T_chain, zero_T_replica);
    // Minimizzazione sulla copia
    iterative_minimization(myParams, zero_T_replica, dsfmt, rank);
    // Rinormalizzao tutti gli spin
    renormalize_spin(zero_T_replica, myParams);
    // Misuro l'energia
    temp=get_energy(zero_T_replica, myParams);
    if(temp<min_energy){
      min_energy=temp;
      // Scrivo la configurazione di bassa temperatura
      // Comprese le misure
      save_zeroT_config(myParams, zero_T_replica);
    }
  }
}
