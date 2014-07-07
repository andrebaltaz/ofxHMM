//
//  HMM.cpp
//  eventsExample
//
//  Created by Andre Baltazar on 10/12/13.
//
//
/*
 * Copyright (c) 2009, Chuan Liu <chuan@cs.jhu.edu>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#include "HMM.h"

#define handle_error(msg) \
do { perror(msg); exit(EXIT_FAILURE); } while (0)

#define IDX(i,j,d) (((i)*(d))+(j))

testHMM::testHMM(vector<vector<int> >vect){

    vec_data=vect;
    setup();
    //Run();

}

//--------------------------------------------------------------
void testHMM::setup(){
    
    nstates = 8;                /* number of states */
    nobvs = 12;                  /* number of observations */
    nseq = vec_data.size();                   /* number of data sequences  */
    length = vec_data[0].size();                 /* data sequencel length */
    //nseq=10;
    //length=10;
    iterations = 30;
    mode = 3;
    
    
    prior = (double *) malloc(sizeof(double) * nstates);
   // if (prior == NULL) handle_error();
    
    trans = (double *) malloc(sizeof(double) * nstates * nstates);
   // if (trans == NULL) handle_error("malloc");
    
    xi = (double *) malloc(sizeof(double) * nstates * nstates);
   // if (xi == NULL) handle_error("malloc");
    
    pi_ = (double *) malloc(sizeof(double) * nstates);
   // if (pi_ == NULL) handle_error("malloc");
    
    obvs = (double *) malloc(sizeof(double) * nstates * nobvs);
   // if (obvs == NULL) handle_error("malloc");
    
    gmm = (double *) malloc(sizeof(double) * nstates * nobvs);
   // if (gmm == NULL) handle_error("malloc");
    
    data = (int *) malloc (sizeof(int) * nseq * length);
   // if (data == NULL) handle_error("malloc");

  
    //initial state probabilities
       for (j = 0; j < nstates; j++) {
        prior[j] = log(1.0/nstates);
       // printf("prior %d=%f\n", j, prior[j]);
    }
    
    /* initial state transition  probabilities */
    int k=0;
    int ax=0;
    for (i=0; i<= nstates; i++) {
        
        k=i*nstates+(ax);
        ax=i+1;
        
        for (j=0; j<nstates; j++) {
           
            if (i==nstates-1 && j==(nstates-1)) {
                trans[IDX( i,j, nstates)] =log(1);
            }else
            
            if (IDX( i,j, nstates)==k) {
                trans[IDX( i,j, nstates)] = log(0.5);
                trans[IDX( i,j+1, nstates)] = log(0.5);
                j++;
            }
            else
            //IDX(i,j,d) (((i)*(d))+(j))
            trans[IDX( i,j, nstates)] =log(0.0001);
            
            
        }
    }
    
    /* state output probabilities */
    for (i=3+nstates; i<= 2+nstates*2 ; i++) {
        for (j=0; j<nobvs; j++) {
            obvs[IDX((i - 3 - nstates),j,nobvs)] = log(1.0/nobvs);
        }
    }

    //setup data
    
    for (i=0; i< nseq; i++) {
        for (j=0; j<length; j++) {
            data[i*length+j] = vec_data[i][j];
            //data[i*length+j] = aux_data[i*length+j];
        }
    }
    
    ////////////////////////////////////////////////////////////////////
    //PRINT INITIAL CONDITIONS
    printf("MODE %d\n\n", mode);
    
    printf("nstates= %d \n", nstates);
    printf("\n");
    printf("Nobvs=%d\n", nobvs);
    printf("\n");
    printf("# initial state probability\n");
    for (j = 0; j < nstates; j++) {
        printf(" %.4f", exp(prior[j]));
    }
    printf("\n");
    printf("\n");
    printf("# state transition probability\n");
    for (j = 0; j < nstates; j++) {
        for (k = 0; k < nstates; k++) {
            printf(" %.4f", exp(trans[IDX(j,k,nstates)]));
        }
        printf("\n");
    }
    printf("\n");
    printf("# initial state output probility\n");
    for (j = 0; j < nstates; j++) {
        for (k = 0; k < nobvs; k++) {
            printf(" %.4f", exp(obvs[IDX(j,k,nobvs)]));
        }
        printf("\n");
    }
    printf("\n");
    
    //printf the data vector
    printf("# data\n");
    
    for (i=0; i< nseq; i++) {
        printf("\n");
        for (j=0; j<length; j++) {
           printf ("data %d = %d\n",i*length+j, data[i*length+j]);
        }
    }

    printf("\n");

    
}

//--------------------------------------------------------------
void testHMM::update(){
    
}

//--------------------------------------------------------------
void testHMM::draw(){
    
}
    
    
void testHMM::Run(int key, vector<vector<int> >vect){
    
    vec_data=vect;
    mode=key;
    
    if (mode==1) {
        for (i=0; i< nseq; i++) {
            for (j=0; j<length; j++) {
                data[i*length+j] = vec_data[i][j];
            }
        }
        nseq=vec_data.size();
        length=vec_data[0].size();
    }
    
    ///////////////////////////////////////////////////////////////////
    
    
    if (mode == 3) {
        loglik = (double *) malloc(sizeof(double) * nseq);
      //  if (loglik == NULL) handle_error("malloc");
        for (i = 0; i < iterations; i++) {
            init_count();
            for (j = 0; j < nseq; j++) {
                loglik[j] = forward_backward(data + length * j, length, 1);
            }
            p = sum(loglik, nseq);
            
            update_prob();
            
            printf("iteration %d log-likelihood: %.4lf\n", i + 1, p);
            printf("updated parameters:\n");
            printf("# initial state probability\n");
            for (j = 0; j < nstates; j++) {
                printf(" %.4f", exp(prior[j]));
            }
            printf("\n");
            printf("# state transition probability\n");
            for (j = 0; j < nstates; j++) {
                for (k = 0; k < nstates; k++) {
                    printf(" %.4f", exp(trans[IDX(j,k,nstates)]));
                }
                printf("\n");
            }
            printf("# state output probility\n");
            for (j = 0; j < nstates; j++) {
                for (k = 0; k < nobvs; k++) {
                    printf(" %.4f", exp(obvs[IDX(j,k,nobvs)]));
                }
                printf("\n");
            }
            printf("\n");
        }
        free(loglik);
//////////////////////////////////////////////////////////////////////////////////////
    } else if (mode == 2) {
        for (i = 0; i < nseq; i++) {
            viterbi(data + length * i, length);
        }
        
        
////////////////////////////////////////////////////////////////////////////////////////
    } else if (mode == 1) {
        loglik = (double *) malloc(sizeof(double) * nseq);
     //   if (loglik == NULL) handle_error("malloc");
        for (i = 0; i < nseq; i++) {
            loglik[i] = forward_backward(data + length * i, length, 0);
        }
        p = sum(loglik, nseq);
        for (i = 0; i < nseq; i++)
            printf("%.4lf\n", loglik[i]);
        printf("total: %.4lf\n", p);
        free(loglik);
    }
    
    //freeall();
    return 0;
}

/* compute sum of the array using Kahan summation algorithm */
double testHMM::sum(double *data, int size)
{
    double sum = data[0];
    int i;
    double y, t;
    double c = 0.0;
    for (i = 1; i < size; i++) {
        y = data[i] - c;
        t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}

/* initilize counts */
void testHMM::init_count() {
    size_t i;
    for (i = 0; i < nstates * nobvs; i++)
        gmm[i] = - INFINITY;
    
    for (i = 0; i < nstates * nstates; i++)
        xi[i] = - INFINITY;
    
    for (i = 0; i < nstates; i++)
        pi_[i] = - INFINITY;
}

void testHMM::update_prob() {
    double pisum = - INFINITY;
    double gmmsum[nstates];
    double xisum[nstates];
    size_t i, j;
    
    for (i = 0; i < nstates; i++) {
        gmmsum[i] = - INFINITY;
        xisum[i] = - INFINITY;
        
        pisum = logadd(pi_[i], pisum);
    }
    
    for (i = 0; i < nstates; i++) {
        prior[i] = pi_[i] - pisum;
    }
    
    for (i = 0; i < nstates; i++) {
        for (j = 0; j < nstates; j++) {
            xisum[i] = logadd(xisum[i], xi[IDX(i,j,nstates)]);
        }
        for (j = 0; j < nobvs; j++) {
            gmmsum[i] = logadd(gmmsum[i], gmm[IDX(i,j,nobvs)]);
        }
    }
    
    for (i = 0; i < nstates; i++) {
        for (j = 0; j < nstates; j++) {
            trans[IDX(i,j,nstates)] = xi[IDX(i,j,nstates)] - xisum[i];
        }
        for (j = 0; j < nobvs; j++) {
            obvs[IDX(i,j,nobvs)] = gmm[IDX(i,j,nobvs)] - gmmsum[i];
        }
    }
    
}

/* forward backward algoritm: return observation likelihood */
double testHMM::forward_backward(int *data, size_t len, int backward)
{
    /* construct trellis */
    double alpha[len][nstates];
    double beta[len][nstates];
    
    size_t i, j, k;
    double p, e;
    double loglik;
    
    for (i = 0; i < len; i++) {
        for (j = 0; j < nstates; j++) {
            alpha[i][j] = - INFINITY;
            beta[i][j] = - INFINITY;
        }
    }
    
    /* forward pass */
    for (i = 0; i < nstates; i++) {
        alpha[0][i] = prior[i] + obvs[IDX(i,data[0],nobvs)];
    }
    for (i = 1; i < len; i++) {
        for (j = 0; j < nstates; j++) {
            for (k = 0; k < nstates; k++) {
                p = alpha[i-1][k] + trans[IDX(k,j,nstates)] + obvs[IDX(j,data[i],nobvs)];
                alpha[i][j] = logadd(alpha[i][j], p);
            }
        }
    }
    loglik = -INFINITY;
    for (i = 0; i < nstates; i++) {
        loglik = logadd(loglik, alpha[len-1][i]);
    }
    
    if (! backward)
        return loglik;
    
    /* backward pass & update counts */
    for (i = 0; i < nstates; i++) {
        beta[len-1][i] = 0;         /* 0 = log (1.0) */
    }
    for (i = 1; i < len; i++) {
        for (j = 0; j < nstates; j++) {
            
            e = alpha[len-i][j] + beta[len-i][j] - loglik;
            gmm[IDX(j,data[len-i],nobvs)] = logadd(gmm[IDX(j,data[len-i],nobvs)], e);
            
            for (k = 0; k < nstates; k++) {
                p = beta[len-i][k] + trans[IDX(j,k,nstates)] + obvs[IDX(k,data[len-i],nobvs)];
                beta[len-1-i][j] = logadd(beta[len-1-i][j], p);
                
                e = alpha[len-1-i][j] + beta[len-i][k]
                + trans[IDX(j,k,nstates)] + obvs[IDX(k,data[len-i],nobvs)] - loglik;
                xi[IDX(j,k,nstates)] = logadd(xi[IDX(j,k,nstates)], e);
            }
        }
    }
    p = -INFINITY;
    for (i = 0; i < nstates; i++) {
        p = logadd(p, prior[i] + beta[0][i] + obvs[IDX(i,data[0],nobvs)]);
        
        e = alpha[0][i] + beta[0][i] - loglik;
        gmm[IDX(i,data[0],nobvs)] = logadd(gmm[IDX(i,data[0],nobvs)], e);
        
        pi_[i] = logadd(pi_[i], e);
    }
    
#ifdef DEBUG
    /* verify if forward prob == backward prob */
    if (fabs(p - loglik) > 1e-5) {
        fprintf(stderr, "Error: forward and backward incompatible: %lf, %lf\n", loglik, p);
    }
#endif
    
    return loglik;
}

/* find the most probable sequence */
void testHMM::viterbi(int *data, size_t len)
{
    double lambda[len][nstates];
    int backtrace[len][nstates];
    int stack[len];
    
    size_t i, j, k;
    double p;
    
    for (i = 0; i < len; i++) {
        for (j = 0; j < nstates; j++) {
            lambda[i][j] = - INFINITY;
        }
    }
    
    for (i = 0; i < nstates; i++) {
        lambda[0][i] = prior[i] + obvs[IDX(i,data[0],nobvs)];
        backtrace[0][i] = -1;       /* -1 is starting point */
    }
    for (i = 1; i < len; i++) {
        for (j = 0; j < nstates; j++) {
            for (k = 0; k < nstates; k++) {
                p = lambda[i-1][k] + trans[IDX(k,j,nstates)] + obvs[IDX(j,data[i],nobvs)];
                if (p > lambda[i][j]) {
                    lambda[i][j] = p;
                    backtrace[i][j] = k;
                }
            }
        }
    }
    
    /* backtrace */
    for (i = 0; i < nstates; i++) {
        if (i == 0 || lambda[len-1][i] > p) {
            p = lambda[len-1][i];
            k = i;
        }
    }
    stack[len - 1] = k;
    for (i = 1; i < len; i++) {
        stack[len - 1 - i] = backtrace[len - i][stack[len - i]];
    }
    for (i = 0; i < len; i++) {
        printf("%d ", stack[i]);
    }
    printf("\n");
}

double testHMM::logadd(double x, double y) {
    if (y <= x)
        return x + log1p(exp(y - x));
    else
        return y + log1p(exp(x - y));
}

void testHMM::usage() {
    fprintf(stdout, "hmm [-hnt] [-c config] [-p(1|2|3)]\n");
    fprintf(stdout, "usage:\n");
    fprintf(stdout, "  -h   help\n");
    fprintf(stdout, "  -c   configuration file\n");
    fprintf(stdout, "  -t   output computation time\n");
    fprintf(stdout, "  -p1  compute the probability of the observation sequence\n");
    fprintf(stdout, "  -p2  compute the most probable sequence (Viterbi)\n");
    fprintf(stdout, "  -p3  train hidden Markov mode parameters (Baum-Welch)\n");
    fprintf(stdout, "  -n   number of iterations\n");
}

/* free all memory */
void testHMM::freeall() {
    if (trans) free(trans);
    if (obvs) free(obvs);
    if (prior) free(prior);
    if (data) free(data);
    if (gmm) free(gmm);
    if (xi) free(xi);
    if (pi_) free(pi_);
}


////////////Manual testing/////////////////

//initial state probabilities

//double prior[ ] = {log(0.04), log(0.02), log(0.06), log(0.04), log(0.11), log(0.11), log(0.01), log(0.09),
//                log(0.03), log(0.05), log(0.06), log(0.11), log(0.05), log(0.11), log(0.03), log(0.08) };




/* initial state transition  probabilities */

//    double aux []= {
//        0.08, 0.02, 0.10, 0.05, 0.07, 0.08, 0.07, 0.04, 0.08, 0.10, 0.07, 0.02, 0.01, 0.10, 0.09, 0.01,
//        0.06, 0.10, 0.11, 0.01, 0.04, 0.11, 0.04, 0.07, 0.08, 0.10, 0.08, 0.02, 0.09, 0.05, 0.02, 0.02,
//        0.08, 0.07, 0.08, 0.07, 0.01, 0.03, 0.10, 0.02, 0.07, 0.03, 0.06, 0.08, 0.03, 0.10, 0.10, 0.08,
//        0.08, 0.04, 0.04, 0.05, 0.07, 0.08, 0.01, 0.08, 0.10, 0.07, 0.11, 0.01, 0.05, 0.04, 0.11, 0.06,
//        0.03 ,0.03, 0.08, 0.10, 0.11, 0.04, 0.06, 0.03, 0.03, 0.08, 0.03, 0.07, 0.10, 0.11, 0.07, 0.03,
//        0.02 ,0.05, 0.01, 0.09, 0.05, 0.09, 0.05, 0.12, 0.09, 0.07, 0.01, 0.07, 0.05, 0.05, 0.11, 0.06,
//        0.11, 0.05, 0.10, 0.07, 0.01, 0.08, 0.05, 0.03, 0.03, 0.10, 0.01, 0.10, 0.08, 0.09, 0.07, 0.02,
//        0.03, 0.02, 0.16, 0.01, 0.05, 0.01, 0.14, 0.14, 0.02, 0.05, 0.01, 0.09, 0.07, 0.14, 0.03, 0.01,
//        0.01, 0.09, 0.13, 0.01, 0.02, 0.04, 0.05, 0.03, 0.10, 0.05, 0.06, 0.06, 0.11, 0.06, 0.03, 0.14,
//        0.09, 0.03, 0.04, 0.05, 0.04, 0.03, 0.12, 0.04, 0.07, 0.02, 0.07, 0.10, 0.11, 0.03, 0.06, 0.09,
//        0.09, 0.04, 0.06, 0.06, 0.05, 0.07, 0.05, 0.01, 0.05, 0.10, 0.04, 0.08, 0.05, 0.08, 0.08, 0.10,
//        0.07, 0.06, 0.01, 0.07, 0.06, 0.09, 0.01, 0.06, 0.07, 0.07, 0.08, 0.06, 0.01, 0.11, 0.09, 0.05,
//        0.03, 0.04, 0.06, 0.06, 0.06, 0.05, 0.02, 0.10, 0.11, 0.07, 0.09, 0.05, 0.05, 0.05, 0.11, 0.08,
//        0.04, 0.03, 0.04, 0.09, 0.10, 0.09, 0.08, 0.06, 0.04, 0.07, 0.09, 0.02, 0.05, 0.08, 0.04, 0.09,
//        0.05, 0.07, 0.02, 0.08, 0.06, 0.08, 0.05, 0.05, 0.07, 0.06, 0.10, 0.07, 0.03, 0.05, 0.06, 0.10,
//        0.11, 0.03, 0.02, 0.11, 0.11, 0.01, 0.02, 0.08, 0.05, 0.08, 0.11, 0.03, 0.02, 0.10, 0.01, 0.11
//    };

/* state output probabilities */


//    double aux_p [] = {
//        0.01, 0.99,
//        0.58, 0.42,
//        0.48, 0.52,
//        0.58, 0.42,
//        0.37, 0.63,
//        0.33, 0.67,
//        0.51, 0.49,
//        0.28, 0.72,
//        0.35, 0.65,
//        0.61, 0.39,
//        0.97, 0.03,
//        0.87, 0.13,
//        0.46, 0.54,
//        0.55, 0.45,
//        0.23, 0.77,
//        0.76, 0.24 };

//setup data
//
//    int aux_data[] ={
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//        1, 1, 0, 0, 1, 1, 1, 0, 0, 0,
//        1, 1, 0, 1, 0, 0, 0, 1, 0, 1,
//        1, 1, 1, 1, 1, 0, 1, 1, 1, 0,
//        0, 1, 0, 1, 1, 0, 1, 1, 1, 1,
//        1, 0, 1, 1, 0, 1, 0, 1, 1, 1,
//        1, 0, 1, 1, 1, 1, 0, 0, 1, 1,
//        0, 1, 0, 1, 1, 1, 0, 0, 0, 0,
//        0, 1, 1, 0, 0, 0, 1, 1, 1, 1,
//        0, 1, 1, 0, 0, 0, 0, 1, 1, 0};