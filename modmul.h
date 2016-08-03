#ifndef __MODMUL_H
#define __MODMUL_H

#include  <stdio.h>
#include <stdlib.h>

#include <string.h>
#include    <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#endif

#define limb_bits  GMP_LIMB_BITS
/*multiplication functions*/
void montMul(mpz_t r, const mpz_t x_copy, const mpz_t y_copy,  const mpz_t modulus, const mp_limb_t w);
void my_mul(mpz_t result, const mpz_t x_copy, const mpz_t y, const mpz_t n);

/*exponentiation functions*/
void mont_power(mpz_t result,  const mpz_t x_copy, const mpz_t e, const mpz_t n, const mp_limb_t w);
void mont_power_windowed(mpz_t result,  const mpz_t x_copy, const mpz_t e, const mpz_t n, const unsigned int window);
void slidingWindow(mpz_t result, const mpz_t x_copy, const mpz_t y, const mpz_t n, mp_size_t k);

/*montgomery reduction reduction */
void mont_red(mpz_t r, const mpz_t x, const mpz_t modular, const mp_limb_t w);
void mont_mod( mpz_t result, const mpz_t x_copy, const mpz_t n, const mp_limb_t w);

// functions to compute parameters
void calc_ro(mpz_t ro, const mpz_t n);
void calc_ro_squared(mpz_t ro_s, const mpz_t n);
mp_limb_t calc_omega(const mpz_t n);
void getBinaryRep(const mpz_t number, int** binary, int *binarysize);
void getKbitRep(const mpz_t number, int** kbitsRep, int *kbitSize, int windowSize);

// custom functions on limbs
mp_size_t my_mpn_add(mp_limb_t* rl, mp_limb_t* xl, mp_size_t xs, mp_limb_t* yl, mp_size_t ys);
mp_size_t my_mpn_mul(mp_limb_t* rl, const mp_limb_t* xl, mp_size_t xs, const mp_limb_t* yl, mp_size_t ys);



/*make a limb addition dealing with all the side preomputations
  rl = pointer to the resulted limbs
  xl and yl = pointer to addition terms limbs
  xs and ys = sizes of x and y
  returns rs = the size of r
*/
mp_size_t my_mpn_add(mp_limb_t* rl, mp_limb_t* xl, mp_size_t xs, mp_limb_t* yl, mp_size_t ys){
  mp_limb_t co;
  if(xs >= ys) /*check bigger size*/
    co = mpn_add(rl, xl, xs, yl, ys);
  else
    co = mpn_add(rl, yl, ys, xl, xs);
  mp_size_t rs = (ys > xs) ? ys : xs;
  if(co == 1){
    rs++;
    // rl = (mp_limb_t*)realloc(rl, sizeof(mp_limb_t)*rs);
    rl[rs-1] = 1;
  }
  return rs;
}

/*make a limb multiplication dealing with all the side preomputations
  rl = pointer to the resulted limbs
  xl and yl = pointer to multiplication terms limbs
  xs and ys = sizes of x and y
  returns rs = the size of r
*/
mp_size_t my_mpn_mul(mp_limb_t* rl, const mp_limb_t* xl, mp_size_t xs, const mp_limb_t* yl, mp_size_t ys){
  mp_limb_t co = 0;
  // compare sizes so multiplication hapens with the bigger one first
  if(xs> ys )
    co = mpn_mul(rl, xl, xs, yl, ys );
  else
     co = mpn_mul(rl, yl, ys, xl, xs);
  mp_size_t rs  = xs + ys -1;
  if(co != 0) {
    rs++;
    // rl = (mp_limb_t*)realloc(rl, sizeof(mp_limb_t)*rs);
    rl[rs-1] = co;
  }
  return rs;
}



void getBinaryRep(const mpz_t number, int** binary, int *binarysize){
  (*binarysize) = mpz_sizeinbase (number, 2);
  *binary = (int*) malloc(sizeof(int)*(*binarysize));
  if(binary == NULL) printf("Cannot allocate memory for binary\n");

  for(int i=0; i<(*binarysize);i++){
    (*binary)[i] = mpz_tstbit(number,i);
  }
}

void getKbitRep(const mpz_t number, int** kbitsRep, int *kbitSize, int windowSize){

  int* binaryRep = NULL;
  int binaryLength;
  int temp, kbit=0;

  getBinaryRep(number, &binaryRep, &binaryLength);

  if(windowSize != 0)
    (*kbitSize) = binaryLength/windowSize + 1;
  else{
    printf("windows size can't be 0\n");
    return;
  }

  *kbitsRep = (int*) malloc(sizeof(int)*(*kbitSize));
  if(kbitsRep == NULL) printf("Cannot allocate memory for binary\n");


  for(int i=0; i<binaryLength; i++){
      temp=0;
      for(int j=windowSize-1; j>=0;  j--){
        if( (i+j) < binaryLength){
          if(binaryRep[i+j]==0)
            temp <<= 1;
          else temp = (temp<<1) + 1 ;
        }
      }
      (*kbitsRep)[kbit] = temp % (1<<windowSize);
      i+=windowSize-1;
      kbit++;
  }
  (*kbitSize) = kbit;

  free(binaryRep);
}

// adds 2 integers using ad
void my_add(mpz_t result, const mpz_t x_copy, const mpz_t y_copy, const mpz_t n){
  mpz_t x;      mpz_init(x);    mpz_set(x, x_copy);
  mpz_t y;      mpz_init(y);    mpz_set(y, y_copy);

  while(mpz_cmp(x, n)>= 0) {mpz_sub(x, x, n); printf("da\n");}
  while(mpz_cmp(y, n)>= 0) {mpz_sub(y, y, n);printf("da\n");}

  mp_size_t xs = mpz_size(x);
  mp_size_t ys = mpz_size(y);
  mp_size_t rs = xs>ys ? xs : ys;

  mp_limb_t* rl = mpz_limbs_modify(result, rs+1);
    rs = my_mpn_add(rl, x->_mp_d, xs, y->_mp_d, ys);
  mpz_limbs_finish(result, rs);

  if(mpz_cmp(result, n) >= 0 ) mpz_sub(result, result, n);

  mpz_clears(x, y, NULL);
}

//  result = x*y %n, using multiplication in normal form
void my_mul(mpz_t result, const mpz_t x_copy, const mpz_t y, const mpz_t n){
  mpz_t x;   mpz_init_set(x, x_copy);

  if(mpz_cmp(x, n)>= 0){
    mp_limb_t w = calc_omega(n);
    mont_red(x, x, n, w);
  }


  mpz_set_ui(result, 0);

  // get binary rep of y
  int* binaryRep = NULL;
  int binaryLength;
  getBinaryRep(y, &binaryRep, &binaryLength);

  for(int i=binaryLength-1; i>=0; i--){
    my_add(result, result, result, n);

    if(binaryRep[i] == 1){
      my_add(result, result, x, n);
    }
  }

  free(binaryRep);
  mpz_clear(x);
}


/*computes a montgomery multiplication of 2 normal integers resulting in a montgomery form result
  r = x * y * ro^(-1)
  + : make copy of x and y so one can safely call montMul(x, x, x, modular, w);
  - : uses more memory as a result
 */
void montMul(mpz_t r, const mpz_t x_copy, const mpz_t y_copy, const mpz_t modulus, const mp_limb_t w){

  mpz_t x;      mpz_init(x);    mpz_set(x, x_copy);
  mpz_t y;      mpz_init(y);    mpz_set(y, y_copy);
  mpz_set_ui(r, 0);

  mp_size_t mods   = mpz_size(modulus);
  mp_size_t xs     = mpz_size(x);
  mp_size_t temps  = xs + 1;
  mp_size_t temp2s = mods + 1;
  mp_size_t rs     = mpz_size(r);

  mp_limb_t* mod_l  = mpz_limbs_read(modulus);
  mp_limb_t* xl     = mpz_limbs_read(x);
  mp_limb_t* templ  = (mp_limb_t*)malloc(sizeof(mp_limb_t)* temps);
  mp_limb_t* templ2 = (mp_limb_t*)malloc(sizeof(mp_limb_t)* temp2s);
  mp_limb_t* rl     = NULL;
  mp_limb_t yil     = 0;
  mp_limb_t u       = 0;

  for(int i=0; i< mods; i++){
    yil = mpz_getlimbn(y, i);

    // compute u
    u = (yil * x->_mp_d[0] + r->_mp_d[0]) * w;

    // compute x * yi
    temps   = my_mpn_mul(templ, xl, xs, &yil, 1);

    //compute  modulus * u in limbs
    temp2s  = my_mpn_mul(templ2, mod_l, mods, &u, 1);

    // add temp and temp2 to r and loose the first limb
    rl = mpz_limbs_modify(r, (rs) > (temps) ? (rs + temp2s +1) : (temps + temp2s + 1));
      rs = my_mpn_add(rl, templ, temps, rl, rs);
      rs = my_mpn_add(rl, templ2, temp2s, rl, rs);
      mpn_rshift(rl, rl, rs, limb_bits/2);
      mpn_rshift(rl, rl, rs, limb_bits/2);
      if(rs > 0) rs--;
    mpz_limbs_finish(r, rs);
  }

  if( 0 <= mpz_cmp(r, modulus)){
    mpz_sub(r, r, modulus);
  }

  mpz_clear(x);
  mpz_clear(y);
  free(templ); free(templ2);
}


/*computes r = x * ro^(-1) given x, a modulus, and omega
  + : makes copy of x so one can safely call mont_red(x, x, modular, w)
 */
void mont_red(mpz_t r, const mpz_t x_copy, const mpz_t modular, const mp_limb_t w){
  mpz_set(r, x_copy);
  mp_size_t mods  = mpz_size(modular);
  mp_size_t rs    = mpz_size(r);
  mp_size_t temps = mods + 1;

  mp_limb_t* rl    = NULL;
  mp_limb_t* modl  = mpz_limbs_read(modular);
  mp_limb_t* templ = (mp_limb_t*)malloc(sizeof(mp_limb_t)*temps);
  mp_limb_t u      = 0;
  mp_limb_t co     = 0;

  for(int i=0; i<mods; i++){
    // // compute u
    u = w * mpz_getlimbn(r, i);

    // /*temp = u * N*/
    temps = my_mpn_mul(templ, modl, mods, &u, 1);

    /*temp = u * N * b^i    */
    for(int j=0; j<i; j++){
      co = mpn_lshift(templ, templ, temps, limb_bits/2);
      if(co != 0) {
        temps++;
        templ = (mp_limb_t*)realloc(templ, sizeof(mp_limb_t)*temps);
        templ[temps-1] = co;
      }
      co = mpn_lshift(templ, templ, temps, limb_bits/2);
      if(co != 0) {
        temps++;
        templ = (mp_limb_t*)realloc(templ, sizeof(mp_limb_t)*temps);
        templ[temps-1] = co;
      }
    }

    /*r = r + temp */
    rl = mpz_limbs_modify(r, (rs > temps) ? rs+1 : temps+1);
      rs = my_mpn_add(rl, rl, rs, templ, temps);
    mpz_limbs_finish(r, rs);
  }

  /*r = r /  b^(size of N )    */
  rl = mpz_limbs_modify(r, rs);
    for(int i=0; i<mods; i++){
      mpn_rshift(rl, rl, rs, limb_bits/2);
      mpn_rshift(rl, rl, rs, limb_bits/2);
      if(rs>0) rs--;
    }
  mpz_limbs_finish(r, rs);

  free(templ);
}

/*function computes result = x * y % n using montgomery form multiplication
-  not recomended for one off multiplication, use my_mul() instead*/
void mont_mul(mpz_t result, const mpz_t x, const mpz_t y, const mpz_t n){
  mp_limb_t w = calc_omega( n);
  mpz_t ro2;      mpz_init(ro2);      calc_ro_squared(ro2, n);
  mpz_t x1;       mpz_init(x1);
  mpz_t y1;       mpz_init(y1);

  montMul(x1, x, ro2, n, w);
  montMul(y1, y, ro2, n, w);
  montMul(result, x1, y1, n,  w );
  mont_red(result, result, n, w);

  mpz_clears(ro2, x1, y1, NULL);
}



/*function computes result = x ^(e) % n using montgomery form
  + : more efficient than exponentiation in normal form
  + : parameter x is not affected by the function
  - : less eficitent than montgomery windowed exponentiation
  */
void mont_power(mpz_t result,  const mpz_t x_copy, const mpz_t e, const mpz_t n, const mp_limb_t w){
  int* binaryPow = NULL;
  int binaryLength;
  getBinaryRep(e, &binaryPow, &binaryLength);

  mpz_t x;       mpz_init(x);           mpz_set(x, x_copy);
  mpz_t t;       mpz_init(t);
  mpz_t t1;      mpz_init(t1);
  mpz_t x1;      mpz_init(x1);
  mpz_t ro_s;    mpz_init(ro_s);        calc_ro_squared(ro_s, n);

  // compute montgomery form result : result * ro % n with result being 1 in the beggining
  mont_red(t1, ro_s, n, w);

  // compute montgomery form of x : x*ro %n
  montMul(x1, x, ro_s, n, w);

  for(int i=binaryLength-1; i>=0; i--){
    // t1 = t*t*ro
    montMul(t1, t1, t1, n, w);

    if(binaryPow[i] == 1){
      // t1 = t*x*ro
      montMul(t1, t1, x1, n, w);
    }
  }
  // result  = t1 * ro^(-1) = t
  mont_red(result, t1, n, w);

  mpz_clear(x);
  mpz_clear(t);
  mpz_clear(t1);
  mpz_clear(x1);
  mpz_clear(ro_s);
  free(binaryPow);
}


/*function computes result = x ^(e) % n using montgomery form and a  k-windowed  algorithm
  + : parameter x is not affected by the function
  - : less eficient than montgomery sliding windowed exponentiation
  window must be bigger than 0
*/
void mont_power_windowed(mpz_t result,  const mpz_t x_copy, const mpz_t e, const mpz_t n, const unsigned int window){
  int* winRep = NULL;
  int winRepLength;
  int preCompSize = (1<<window)-1;
  getKbitRep(e, &winRep, &winRepLength, window);

  mpz_t ro_s;       mpz_init(ro_s);       calc_ro_squared(ro_s, n);
  mpz_t x;          mpz_init_set(x, x_copy);
  mp_limb_t w  = calc_omega(n);
  mpz_t* x_hat      = (mpz_t*)malloc(sizeof(mpz_t)*preCompSize);


  mpz_init(x_hat[0]);

  mont_red(result, ro_s, n, w);
  montMul(x_hat[0], x, ro_s, n, w);

  for(int i=1; i<preCompSize; i++){
    mpz_init(x_hat[i]); mpz_set_ui(x_hat[i], 1);
    montMul(x_hat[i], x_hat[i-1], x_hat[0], n, w);
  }

  for(int i=winRepLength-1; i>=0; i--){
    for(int j=0; j< window; j++)
      montMul(result, result, result, n, w) ;

    if(winRep[i] != 0)
      montMul(result, result, x_hat[winRep[i] -1 ], n, w);
  }
  mont_red(result, result, n, w);

  mpz_clear(ro_s);
  for(int i=0; i<preCompSize; i++)
    mpz_clear(x_hat[i]);
  free(x_hat); free(winRep);
  mpz_clear(x);
}

/*function computes result = x ^(e) % n using montgomery form and a sliding window ( of size max size k) algorithm
   + : parameter x is not affected by the function
   window must be bigger than 0
*/
void slidingWindow(mpz_t result, const mpz_t x_copy, const mpz_t y, const mpz_t n, mp_size_t k){
  mpz_t x;   mpz_init_set(x, x_copy);

  // compute binary representation of y
  int* binary = NULL;
  int binaryLength;
  getBinaryRep(y, &binary, &binaryLength);

  int l, count =0, i = binaryLength-1;
  int preCompSize = (1<<k) >> 1;

  // compute omega
  mp_limb_t w = calc_omega(n);
  mp_limb_t u;

  // compute ro squared
  mpz_t ro_squared;     mpz_init(ro_squared);  calc_ro_squared(ro_squared, n);

  // allocate memory for elements to be procomputed
  mpz_t* x_hat = (mpz_t*) malloc(sizeof(mpz_t)*preCompSize);

  // compute montgomery form of result ( = ro at the beginning)
  mont_red(result, ro_squared, n, w);

  // compute the first precomputation ( happens even if window k  = 0)
  mpz_init(x_hat[0]);         montMul(x_hat[0], x, ro_squared, n, w);

  // compute montgomery form of x squared
  mpz_t x_sq1;         mpz_init(x_sq1);            montMul(x_sq1, x_hat[0], x_hat[0], n, w);

  //make the precomputations array with all elements in montgomery form
  for(int i=1; i<preCompSize; i++){
    mpz_init(x_hat[i]);
    montMul(x_hat[i], x_hat[i-1], x_sq1, n, w);
  }
  mpz_clears(ro_squared, x_sq1, NULL);

  /*while digits still left*/
  while(i >=0){
    /*dont use any window */
    if(binary[i] == 0){
      l = i;
      u = 0;
    }
    /*decide how big to make the window */
    else{
      l = (i - k + 1) > 0 ? i-k+1  : 0;
      count = 0;
      u = 0;
      while(binary[l]==0 && l < i ) l++;
      for(count=i; count>=l; count--){
        u = u*2;
        if(binary[count] == 1)
          u++;
      }
    }
    /*use the resulted window to compute the result*/
    for(count = 0; count< i - l + 1; count++){
      montMul(result, result, result, n, w);
    }
    if(u != 0)
      montMul(result, result, x_hat[(u-1)>>1], n, w);
    i =  l-1;
  }
  /*return the result into normal form */
  mont_red(result, result, n, w);

  free(binary);
  for(int i=0; i<preCompSize; i++)
    mpz_clear(x_hat[i]);
  free(x_hat);
  mpz_clear(x);
}

/*computes x % n using montgomery form */
void mont_mod( mpz_t result, const mpz_t x_copy, const mpz_t n, const mp_limb_t w){
  if(mpz_cmp(x_copy, n) < 0 ){
    mpz_set(result, x_copy);
    return;
  }
  if(mpz_cmp(x_copy, n) == 0){
    mpz_set(result, 0);
  }

  mpz_t x;            mpz_init(x);            mpz_set(x, x_copy);
  mpz_t ro_squared;   mpz_init(ro_squared);   calc_ro_squared(ro_squared, n);

  mont_red(x, x, n, w);
  montMul(result, x, ro_squared, n, w);

  mpz_clear(ro_squared);
  mpz_clear(x);
}


// computes omega for a modulus n
mp_limb_t calc_omega(const mpz_t n){
  mp_limb_t size = 1;
  mpn_lshift(&size, &size, 1, limb_bits-1);
  size -= 1;

  mp_limb_t w1 = - mpz_getlimbn(n, 0);
  mp_limb_t w  = 1;

  for(int i=0; i<limb_bits; i++){
    w = w * w * w1;
  }
  return w;
}

// compute ro for a modulus n
void calc_ro(mpz_t ro, const mpz_t n){
  mp_size_t ros = mpz_size(n) + 1;
  mp_limb_t* rol = mpz_limbs_modify(ro, ros);
    for(int i=0; i<ros-1; i++)
      rol[i] = 0;
    rol[ros-1] = 1;
  mpz_limbs_finish(ro, ros);
}

// compute ro squared for a modulus independently of ro
void calc_ro_squared(mpz_t ro_s, const mpz_t n){
  mp_size_t bits = 2 * mpz_size(n) * limb_bits;
  mpz_set_ui(ro_s, 1);
  for(int i=0; i< bits; i++){
    mpz_mul_2exp(ro_s, ro_s, 1);
    if(mpz_cmp(ro_s, n)>= 0)
      mpz_sub(ro_s, ro_s, n);
  }
}

// prints limbs of a mpz integer
void print_limbs(mpz_t x){
  mp_size_t xs = mpz_size(x);
  for(int i=0; i<xs; i++)
    gmp_printf("limb %d is %llu\n", i,  x->_mp_d[i]);
}
