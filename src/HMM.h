//
//  HMM.h
//  eventsExample
//
//  Created by Andre Baltazar on 10/12/13.
//
//

#ifndef __eventsExample__HMM__
#define __eventsExample__HMM__

#include "ofMain.h"


class testHMM : public ofBaseApp{
	
public:
    
    testHMM(vector<vector<int> >vect);
    
    void setup();
    void update();
    void draw();
    void Run(int key, vector<vector<int> >vect);
       
    int nstates;                /* number of states */
    int nobvs;                  /* number of observations */
    int nseq;                   /* number of data sequences  */
    int length;                 /* data sequencel length */
    double *prior;           /* initial state probabilities */
    double *trans;           /* state transition probabilities */
    double *obvs;            /* output probabilities */
    int *data;
    double *gmm;             /* gamma */
    double *xi;              /* xi */
    double *pi_;              /* pi_ */
    
    double logadd(double, double);
    double sum(double *, int);
    double forward_backward(int *, size_t, int);
    void viterbi(int *, size_t);
    void init_count();
    void update_prob();
    void usage();
    void freeall();

    int iterations;
    int mode;
    
    int c;
    double d;
    double *loglik;
    double p;
    int i, j, k;
    
    
    vector<vector<int> >vec_data;
    
    
    
    
    
    };














#endif /* defined(__eventsExample__HMM__) */
