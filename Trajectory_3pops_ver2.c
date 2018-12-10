// nreps:   number of replicates.
// h: dominance effect.
// nsamK: number of samples in KhoeSan.
// nderK: number of derived alleles in KhoeSan.
// alpha: mean of prior on selection coefficient (exponential distribution).
// t: time when focal allele appears.
// mEP: fraction of European ancestry in Pastoralists.
// tEP: time of a single-wave migration from Europeans to Pastoralists (generatios ago).
// mEK: fraction of European ancestry in Khoesan.
// tEK: time of a single-wave migration from Europeans to Khoesan (generatios ago).
// mPK: fraction of Pastoralist ancestry in Khoesan.
// tPK: time of a single-wave migration from Pastoralists to Khoesan (generatios ago).
// nsamP: number of samples in KhoeSan.
// nderP: number of derived alleles in KhoeSan.
//// The function popsize(j), returns the diploid population size at time j generations in the past.
//// The user must supply this function. Examples are provided.
//// A Wright-Fisher model is assumed.  
//// The first line of the output is nreps.
//// The following line contains "npoints: n1 (where n1 is the number of time points until present)."
//// Following that are n1 lines where each line contains two numbers, the time point( = generation / (4 * N)), allele frequency.
//// Time is measured in units of 4N generations.
//// Compilation: gcc -o traj -lm Trajectory_3pops_ver2.c popsizeE.c popsizeP.c popsizeK.c binomial_mt.c mt.c
////              gcc -o stepftn stepftn_3pops.c
//// usage: "./traj nreps h nsamK nderK alpha t mEP tEP mEK tEK mPK tPK nsamP nderP seed | ./stepftn > XXX.out"

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

extern int npres    ;
extern double r ;

# define MAX 1000000
# define MAXLEN 1000

main(int argc, char **argv) {
    
    
    int nreps, rep, i, popmax, count, nsuccess, seed   ;
    int nsamK, nderK, nsamP, nderP  ;
    int t, tEP, tEK, tPK, genK, genP, genE  ;

    double s, h, alpha, pK, pE, pP, pEfwd, pPfwd, pKfwd, wbar, pprime, pKmax, pPmax, prob, ngen  ;
    double *freqE, *freqK, *freqP   ;
    double sE, mEP, mEK, mPK   ;

    int popsizeK(int)    ;      // Population size changes in Khoesan.
    int popsizeP(int)    ;      // Population size changes in Pastoralists.
    int popsizeE(int)    ;      // Population size changes in Europeans.
    double genrand_real3(void)  ;
    float bnldev(float,int) ;
    double expdist(double lam)  ;

    if(argc < 15) {
        fprintf(stderr, "./traj nreps h nsamK nderK alpha t mEP tEP mEK tEK mPK tPK nsamP nderP seed\n")    ;
        exit(1) ;
    }

    printf("// ")   ;
    for(i = 0; i < argc; i++)   printf(" %s", argv[i])    ;
    printf("\n")    ;
    
        nreps = strtol(argv[1], NULL, 10)   ;    // number of replicates.
    // Parameters in Khoesan.
        h = strtod(argv[2], NULL)   ;            // dominance effect.
        nsamK = strtol(argv[3], NULL, 10)    ;   // number of samples in KhoeSan.
        nderK = strtol(argv[4], NULL, 10)    ;   // number of derived alleles in KhoeSan.
        alpha = strtod(argv[5], NULL)  ;         // mean of prior on selection coefficient (exponential distribution).
    // Parameters in Europeans.
        t = strtol(argv[6], NULL, 10)  ;         // time when focal allele appears.
        mEP = strtod(argv[7], NULL)  ;           // fraction of European ancestry in Pastoralists.
        tEP = strtol(argv[8], NULL, 10)  ;       // time of a single-wave migration from Europeans to Pastoralists (generatios ago).
        mEK = strtod(argv[9], NULL)  ;           // fraction of European ancestry in Khoesan.
        tEK = strtol(argv[10], NULL, 10)  ;      // time of a single-wave migration from Europeans to Khoesan (generatios ago).
    // Parameters in Pastoralists.
        mPK = strtod(argv[11], NULL)  ;          // fraction of Pastoralist ancestry in Khoesan.
        tPK = strtol(argv[12], NULL, 10)  ;      // time of a single-wave migration from Pastoralists to Khoesan (generatios ago).
        nsamP = strtol(argv[13], NULL, 10)    ;  // number of samples in KhoeSan.
        nderP = strtol(argv[14], NULL, 10)    ;  // number of derived alleles in KhoeSan.
    // Others.
        sscanf(argv[15], "%x", &seed)   ;        // random seed.

    // Initialize Mersenne Twistor.
        void init_by_array(unsigned long init_key[], int key_length)    ;
        unsigned long init[4] = {seed, 0x265, 0x367, 0x569}, length=4 ;
        init_by_array(init, length) ;
    
    freqK = (double *)malloc((unsigned)(MAX+1) * sizeof(double)) ;
    freqP = (double *)malloc((unsigned)(MAX+1) * sizeof(double)) ;
    freqE = (double *)malloc((unsigned)(MAX+1) * sizeof(double)) ;
    pKmax = nderK / (double)nsamK  ;
    pPmax = nderP / (double)nsamP  ;

    popmax = 0  ;
    for(i = 0; i < MAX; i++)    if(popsizeK(i) > popmax) popmax = popsizeK(i) ;
    printf("%d N0: %d\n", nreps, popsizeK(0))    ;
    nsuccess = 0    ;
    count = 0   ;
                    
    while(nsuccess < nreps){
        count++ ;
                
        // Selection started in Europeans t generations ago.
            pE = 1.0 / (2.0 * popsizeE(t))  ;
            pEfwd = pE  ;
            freqE[0] = 0.0  ;
            genE = 1    ;
            freqE[genE] = pEfwd ;
            sE = expdist(alpha)  ;
            while((genE < t) && (pEfwd > 0))    {
                genE++  ;
                wbar = pEfwd * pEfwd * (1. + sE) + 2. * pEfwd * (1. - pEfwd) * (1. + sE * h) + (1. - pEfwd) * (1. - pEfwd)   ;
                pprime =( pEfwd * pEfwd * (1. + sE) + pEfwd * (1. - pEfwd) * (1. + sE * h)) / wbar  ;
                pEfwd = bnldev(pprime, 2 * popsizeE(t - genE)) / (2.0 * popsizeE(t - genE))  ;
                freqE[genE] = pEfwd ;
            }

            // Accept if the allele is fixed in Europeans.
            if(pEfwd == 1.0 && genE == t)    {
                           
                // Neutral in Pastoralists.
                    pP = mEP * freqE[t - tEP]   ;
                    pPfwd = pP  ;
                    genP = (t - tEP)    ;
                    for(i =0; i <  (t - tEP); i++)    freqP[i] = 0.0  ;    
                    freqP[genP] = pPfwd  ;
                    while((genP < t) && (pPfwd > 0.0) && (pPfwd < 1.0))   {
                        genP++   ;
                        pprime = pPfwd   ;
                        pPfwd = bnldev(pprime, 2 * popsizeP(t - genP)) / (2.0 * popsizeP(t - genP))    ;
                        freqP[genP] = pPfwd ;
                    }
                    prob = pow(pPfwd / pPmax, (double)nderP) * pow((1.0 - pPfwd) / (1.0 - pPmax), (double)(nsamP - nderP))   ;
                    if(genrand_real3() < prob)  {
                
                        // Sample selection coefficiet from a prior.
                            s = expdist(alpha)  ;

                        // Selection in KhoeSan.
                            pK = mPK * freqP[t - tPK]    ;
                            pKfwd = pK    ;
                            genK = (t - tPK) ;     
                            for(i =0; i <  (t - tPK); i++)    freqK[i] = 0.0  ;
                            freqK[genK] = pKfwd ;
                            while((genK < t) && (pKfwd > 0) && (pKfwd < 1.0)) {
                                genK++   ;
                                wbar = pKfwd * pKfwd * (1. + s) + 2. * pKfwd * (1. - pKfwd) * (1. + s * h) + (1. - pKfwd) * (1. - pKfwd)  ;
                                pprime = (pKfwd * pKfwd * (1. + s) + pKfwd * (1. - pKfwd) * (1. + s * h)) / wbar    ;
                            // A single-wave of migration from Europeans to Khoesan.
                                if(t - genK == tEK)  {
                                    pKfwd = bnldev(((1 - mEK) * pprime + mEK * freqE[genK]), 2 * popsizeK(t - genK)) / (2.0 * popsizeK(t - genK))    ;

                                } else {
                                    pKfwd = bnldev(pprime, 2 * popsizeK(t - genK)) / (2.0 * popsizeK(t - genK))    ;

                                }
                                freqK[genK] = pKfwd    ;
                            } 
                    // Accept trajectory.
                        prob = pow(pKfwd / pKmax, (double)nderK) * pow((1.0 - pKfwd) / (1.0 - pKmax), (double)(nsamK - nderK))   ;
                        if(genrand_real3() < prob)  {
                            nsuccess++  ;
                            printf("# s: %lf age: %lf freq: %lf\n", s, (double)tPK / (4.0 * popsizeK(0)), pK) ;
                            for(i = 0; i < t; i++)   printf("%lf\t%.10f\t%.10f\t%.10f\n", i / (4.0 * popsizeK(0)), freqK[i], freqP[i], freqE[i])   ;
                        }
                    }
            }
    }

    free(freqK)  ;    
        freqK = NULL ;    
    free(freqP)  ;    
        freqP = NULL ;    
    free(freqE)  ;    
        freqE = NULL ;    

    fprintf(stderr," !%.10f\n", nreps / (double)count) ;

    return 0;    

}

double expdist(double lam){

    double genrand_real3();  
    return -log(genrand_real3())/lam;
}

