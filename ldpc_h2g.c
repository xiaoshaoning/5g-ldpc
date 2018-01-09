/* Invert sparse binary H for  LDPC*/    
/* Author : Igor Kozintsev   igor@ifp.uiuc.edu   
   Please let me know if you find bugs in this code (I did test   
   it but I still have some doubts). All other comments are welcome    
   too :) !   
   I use a simple algorithm to invert H.   
   We convert H to [I | A]   
                   [junk ]   
   using column reodering and row operations (junk - a few rows of H   
   which are linearly dependent on the previous ones)   
   G is then found as G = [A'|I]   
   G is stored as array of doubles in Matlab which is very inefficient.   
   Internal representation in this programm is unsigned char. Please modify   
   the part which writes G if you wish.   
   */    
#include <math.h>    
#include "mex.h"    
    
/* Input Arguments: tentative H matrix*/    
#define H_IN    prhs[0]              
#define Q_IN    prhs[1] /* field base */    
    
/* Output Arguments: final matrices*/    
#define H_OUT   plhs[0]    
#define G_OUT   plhs[1]    
    
    
    
    
    
    
/************************************ GFq math *******************************/    
    
/* This file contains lookup tables and routines required   
 * to perform the field operations over GF(q), i.e. addition   
 * and multiplication.   
 *   
 * Addition is easy, we use exclusive-or operation.   
 * For multiplication we need two tables for each q.   
 * The first is the logarithm table, and the second   
 * is the exponential table.  We catch multiplication   
 * by zero and one separately.   
 *   
 * Source for tables: MacWilliams and Sloane, fig 4.5   
 *   
 * WARNING: for speed, no check is made that legal values   
 *          are supplied.   
 */    
const int log4[4]   = {0,0,1,2}; /* stores i at address \alpha^i */    
const int log8[8]   = {0,0,1,3,2,6,4,5};    
const int log16[16] = {0,0,1,4,2,8,5,10,3,14,9,7,6,13,11,12};    
const int log32[32] = {0,0,1,18,2,5,19,11,3,29,6,27,20,8,12,23,4,\    
            10,30,17,7,22,28,26,21,25,9,16,13,14,24,15};    
const int log64[64] = {0,0,1,6,2,12,7,26,3,32,13,35,8,48,27,18,4,24,\    
            33,16,14,52,36,54,9,45,49,38,28,41,19,56,5,62,\    
            25,11,34,31,17,47,15,23,53,51,37,44,55,40,10,\    
            61,46,30,50,22,39,43,29,60,42,21,20,59,57,58};    
const int log128[128] = {0,0,1,31,2,62,32,103,3,7,63,15,33,84,104,\    
            93, 4,124,8,121,64,79,16,115,34,11,85,38,105,46,94,51,\    
            5,82,125,60,9,44,122,77,65,67,80,42,17,69,116,23,35,118,\    
            12,28,86,25,39,57,106,19,47,89,95,71,52,110,6,14,83,92,126,\    
            30,61,102,10,37,45,50,123,120,78,114,66,41,68,22,81,59,43,76,\    
            18,88,70,109,117,27,24,56,36,49,119,113,13,91,29,101,87,108,\    
            26,55,40,21,58,75,107,54,20,74,48,112,90,100,96,97,72,98,53,73,111,99};    
const int log256[256] = {0,0,1,25,2,50,26,198,3,223,51,238,27,104,199,75,4,100,\    
            224,14,52,141,239,129,28,193,105,248,200,8,76,113,5,138,101,47,225,\    
            36,15,33,53,147,142,218,240,18,130,69,29,181,194,125,106,39,249,185,\    
            201,154,9,120,77,228,114,166,6,191,139,98,102,221,48,253,226,152,37,\    
            179,16,145,34,136,54,208,148,206,143,150,219,189,241,210,19,92,131,\    
            56,70,64,30,66,182,163,195,72,126,110,107,58,40,84,250,133,186,61,202,\    
            94,155,159,10,21,121,43,78,212,229,172,115,243,167,87,7,112,192,247,\    
            140,128,99,13,103,74,222,237,49,197,254,24,227,165,153,119,38,184,180,\    
            124,17,68,146,217,35,32,137,46,55,63,209,91,149,188,207,205,144,135,151,\    
            178,220,252,190,97,242,86,211,171,20,42,93,158,132,60,57,83,71,109,65,\    
            162,31,45,67,216,183,123,164,118,196,23,73,236,127,12,111,246,108,161,59,\    
            82,41,157,85,170,251,96,134,177,187,204,62,90,203,89,95,176,156,169,160,\    
            81,11,245,22,235,122,117,44,215,79,174,213,233,230,231,173,232,116,214,\    
            244,234,168,80,88,175};    
    
const int exp4[3]   = {1,2,3}; /* stores \alpha^i at address i */    
const int exp8[7]   = {1,2,4,3,6,7,5};    
const int exp16[15] = {1,2,4,8,3,6,12,11,5,10,7,14,15,13,9};    
const int exp32[31] = {1,2,4,8,16,5,10,20,13,26,17,7,14,28,29,31,\    
            27,19,3,6,12,24,21,15,30,25,23,11,22,9,18};    
const int exp64[63] = {1,2,4,8,16,32,3,6,12,24,48,35,5,10,20,40,19,\    
            38,15,30,60,59,53,41,17,34,7,14,28,56,51,37,\    
            9,18,36,11,22,44,27,54,47,29,58,55,45,25,50,\    
            39,13,26,52,43,21,42,23,46,31,62,63,61,57,49,33};    
const int exp128[127] = {1,2,4,8,16,32,64,9,18,36,72,25,50,100,65,11,\    
            22,44,88,57,114,109,83,47,94,53,106,93,51,102,69,3,6,12,24,\    
            48,96,73,27,54,108,81,43,86,37,74,29,58,116,97,75,31,62,124,\    
            113,107,95,55,110,85,35,70,5,10,20,40,80,41,82,45,90,61,122,\    
            125,115,111,87,39,78,21,42,84,33,66,13,26,52,104,89,59,118,101,\    
            67,15,30,60,120,121,123,127,119,103,71,7,14,28,56,112,105,91,63,\    
            126,117,99,79,23,46,92,49,98,77,19,38,76,17,34,68};    
    
const int exp256[255] = {1,2,4,8,16,32,64,128,29,58,116,232,205,135,19,38,76,\    
            152,45,90,180,117,234,201,143,3,6,12,24,48,96,192,157,39,78,156,\    
            37,74,148,53,106,212,181,119,238,193,159,35,70,140,5,10,20,40,80,\    
            160,93,186,105,210,185,111,222,161,95,190,97,194,153,47,94,188,101,\    
            202,137,15,30,60,120,240,253,231,211,187,107,214,177,127,254,\    
            225,223,163,91,182,113,226,217,175,67,134,17,34,68,136,13,26,52,104,\    
            208,189,103,206,129,31,62,124,248,237,199,147,59,118,236,197,151,51,\    
            102,204,133,23,46,92,184,109,218,169,79,158,33,66,132,21,42,84,168,\    
            77,154,41,82,164,85,170,73,146,57,114,228,213,183,115,230,209,191,99,\    
            198,145,63,126,252,229,215,179,123,246,241,255,227,219,171,75,150,49,\    
            98,196,149,55,110,220,165,87,174,65,130,25,50,100,200,141,7,14,28,56,\    
            112,224,221,167,83,166,81,162,89,178,121,242,249,239,195,155,43,86,172,\    
            69,138,9,18,36,72,144,61,122,244,245,247,243,251,235,203,139,11,22,44,\    
            88,176,125,250,233,207,131,27,54,108,216,173,71,142};    
    
    
/* For testing:   
   main(){   
      
   u_int i,j,q;   
      
   while(1){   
   printf("please enter a b q:");   
   scanf("%u %u %u",&i,&j,&q);   
   printf("\n a * b in GF(q) = %u\n",GFq_m(i,j,q));   
   }   
   }   
   */    
    
int GFq_m(int a, int b, int q)    
{    
      
  if ( a == 0 || b == 0 )  return 0 ;     
  if ( a == 1 )  return b ;    
  if ( b == 1 )  return a ;    
  switch (q){    
  case 256:    
    return exp256[(log256[a]+log256[b])%255];    
  case 128:    
    return exp128[(log128[a]+log128[b])%127];    
  case 64:    
    return exp64[(log64[a]+log64[b])%63];    
  case 32:    
    return exp32[(log32[a]+log32[b])%31];    
  case 16:    
    return exp16[(log16[a]+log16[b])%15];    
  case 8:    
    return exp8[(log8[a]+log8[b])%7];    
  case 4:    
    return exp4[(log4[a]+log4[b])%3];    
  }    
  mexErrMsgTxt(1,"GFq_m: I'm afraid I don't know how to multiply in GFq\n");    
  return 0 ;    
}    
    
int GFq_inv(int a, int q)    
{    
      
  if ( a == 0) mexErrMsgTxt(1,"GFq_inv: no inverse for 0!\n");;     
  if ( a == 1 )  return 1 ;    
  switch (q){    
  case 256:    
    return exp256[(255-log256[a])];    
  case 128:    
    return exp128[(127-log128[a])];    
  case 64:    
    return exp64[(63-log64[a])];    
  case 32:    
    return exp32[(31-log32[a])];    
  case 16:    
    return exp16[(15-log16[a])];    
  case 8:    
    return exp8[(7-log8[a])];    
  case 4:    
    return exp4[(3-log4[a])];    
  }    
  mexErrMsgTxt(1,"GFq_inv: not defined inverse for  GFq\n");    
  return 0 ;    
}    
    
int GFq_a(int a, int b)    
{    
    return a^b;    
}    
    
/************************************ end GFq math *******************************/    
    
    
    
    
    
    
    
void mexFunction(    
                 int nlhs,       mxArray *plhs[],    
                 int nrhs, const mxArray *prhs[]    
          )    
{    
  unsigned char **HH, **GG;    
  int ii, jj, *ir, *jc, rdep, tmp, d, q, scale;    
  double *sr1, *sr2, *g;    
  int N,M,K,i,j,k,kk,nz,*irs1,*jcs1, *irs2, *jcs2;    
    
  /* Check for proper number of arguments */    
  if (nrhs != 2) {    
    mexErrMsgTxt("h2g requires two input arguments.");    
  } else if (nlhs != 2) {    
    mexErrMsgTxt("h2g requires two output arguments.");    
  } else if (!mxIsSparse(H_IN)) {    
    mexErrMsgTxt("h2g requires sparse H matrix.");    
  }    
      
/* get the field base */    
  q = (int)mxGetScalar(Q_IN);    
    
/* read sparse matrix H */    
    sr1  = mxGetPr(H_IN);    
    irs1 = mxGetIr(H_IN);  /* row */    
    jcs1 = mxGetJc(H_IN);  /* column */    
    nz = mxGetNzmax(H_IN); /* number of nonzero elements (they are ones)*/    
    M = mxGetM(H_IN);       
    N = mxGetN(H_IN);    
    
        
/* create working array HH[row][column]*/    
    HH = (unsigned char **)mxMalloc(M*sizeof(unsigned char *));    
    for(i=0 ; i<M ; i++){    
      HH[i] = (unsigned char *)mxMalloc(N*sizeof(unsigned char));    
    }    
    for(i=0 ; i<M ; i++)    
      for(j=0 ; j<N ; j++)    
        HH[i][j] = 0; /* initialize all to zero */    
    
    k=0;    
    for(j=0 ; j<N ; j++) {    
      for(i=0 ; i<(jcs1[j+1]-jcs1[j]) ; i++) {    
        ii = irs1[k]; /* index in column j*/     
        HH[ii][j] = (unsigned char)sr1[k]; /* put  nonzeros */    
        k++;    
      }    
    }    
    
/* invert HH matrix here */    
    /* row and column indices */    
    ir = (int *)mxMalloc(M*sizeof(int));    
        jc = (int *)mxMalloc(N*sizeof(int));    
    for( i=0 ; i<M ; i++)    
        ir[i] = i;    
    for( j=0 ; j<N ; j++)    
        jc[j] = j;    
    
    
    
    /* perform Gaussian elimination on H, store reodering operations */    
    rdep = 0; /* number of dependent rows in H*/    
    d = 0;    /* current diagonal element */    
    
    while( (d+rdep) < M) { /* cycle through independent rows of H */    
        
        j = d; /* current column index along row ir[d] */    
        while( (HH[ir[d]][jc[j]] == 0) && (j<(N-1)) )    
            j++;            /* find first nonzero element in row i */    
        if( HH[ir[d]][jc[j]] ) { /* found nonzero element */    
    
    
            /* swap columns */    
            tmp = jc[d]; jc[d] = jc[j]; jc[j] = tmp;    
    
            if(q==2) { /* GF2 */    
              /* eliminate current column using row operations */    
                for(ii=0 ; ii<M ; ii++)     
                 if(HH[ir[ii]][jc[d]] && (ii != d)) /* nonzero and non-diagonal */    
                    for(jj=d ; jj<N ; jj++)     
                     HH[ir[ii]][jc[jj]] = (HH[ir[ii]][jc[jj]]+HH[ir[d]][jc[jj]])%2;    
            }    
            else { /* GFq */    
                    
                scale = GFq_inv(HH[ir[d]][jc[d]],q); /* inverse of the diag. element */    
                /* scale the current row to make the first element 1 */    
                for(jj=0 ; jj<N ; jj++)     
                    HH[ir[d]][jc[jj]] = GFq_m(HH[ir[d]][jc[jj]],scale,q);    
    
                /* eliminate current column using row operations */    
                for(ii=0 ; ii<M ; ii++) {    
                    if(HH[ir[ii]][jc[d]] && (ii != d)) {    
                        scale = HH[ir[ii]][jc[d]];    
                        for(jj=d ; jj<N ; jj++) {    
                            tmp = GFq_m(HH[ir[d]][jc[jj]],scale,q);    
                            HH[ir[ii]][jc[jj]] = GFq_a(HH[ir[ii]][jc[jj]],tmp);    
                        }    
                    }    
                }    
            }    
        }    
        else { /* all zeros -  need to delete this row and update indices */    
            rdep++; /* increase number of dependent rows */    
            tmp = ir[d];    
            ir[d] = ir[M-rdep];    
            ir[M-rdep] = tmp;    
            d--; /* no diagonal element is found */    
        }    
        d++; /* increase the number of diagonal elements */    
      
    }/*while i+rdep*/    
/* done inverting HH */    
    
    
    K = N-M+rdep; /* true K */    
    
/* create G matrix  G = [A'| I] if H = [I|A]*/    
    GG = (unsigned char **)mxMalloc(K*sizeof(unsigned char *));    
    for(i=0 ; i<K ; i++){    
        GG[i] = (unsigned char *)mxMalloc(N*sizeof(unsigned char));    
    }    
    for(i=0 ; i<K ; i++)    
        for(j=0 ; j<(N-K) ; j++) {    
            tmp = (N-K+i);    
            GG[i][j] = HH[ir[j]][jc[tmp]];    
        }    
    
    for(i=0 ; i<K ; i++)    
        for(j=(N-K); j<N ; j++)    
            if(i == (j-N+K) ) /* diagonal */    
                GG[i][j] = 1;    
            else    
                GG[i][j] = 0;    
    
/* NOTE, it is a very inefficient way to store G. Change to taste!*/    
    G_OUT = mxCreateDoubleMatrix(K, N, mxREAL);    
    /* Assign pointers to the output matrix */    
    g = mxGetPr(G_OUT);    
    for(i=0 ; i<K ; i++)    
        for(j=0 ; j<N; j++)    
            g[i+j*K] = GG[i][j];    
    
    
    H_OUT = mxCreateSparse(M,N,nz,mxREAL);    
    sr2  = mxGetPr(H_OUT);    
    irs2 = mxGetIr(H_OUT);  /* row */    
    jcs2 = mxGetJc(H_OUT);  /* column */    
    /* Write H_OUT swapping columns according to jc */    
    k = 0;     
    for (j=0; (j<N ); j++) {    
      jcs2[j] = k;    
      tmp = jcs1[jc[j]+1]-jcs1[jc[j]];    
      for (i=0; i<tmp ; i++) {    
        kk = jcs1[jc[j]]+i;    
        sr2[k] = sr1[kk];    
        irs2[k] = irs1[kk];    
        k++;            
      }    
    }    
    jcs2[N] = k;    
        
    
/* free the memory */    
  for( j=0 ; j<M ; j++) {    
    mxFree(HH[j]);    
  }    
  mxFree(HH);    
  mxFree(ir);    
  mxFree(jc);      
  for(i=0;i<K;i++){    
    mxFree(GG[i]);    
  }    
  mxFree(GG);    
  return;    
}           