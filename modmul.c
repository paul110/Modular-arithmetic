#include "modmul.h"
#include <check.h>

#define limb_bits  GMP_LIMB_BITS
/*
Perform stage 1:

- read each 3-tuple of N, e and m from stdin,
- compute the RSA encryption c, then
- write the ciphertext c to stdout.
*/






void stage1() {
  // struct timeval timstr;        /* structure to hold elapsed time */
  // struct rusage ru;             /* structure to hold CPU time--system and user */
  // double tic,toc;               /* floating point numbers to calculate elapsed wallclock time */
  // double usrtim;                /* floating point number to record elapsed user CPU time */

  // double systim;                /* floating point number to record elapsed system CPU time */
  // gettimeofday(&timstr,NULL);
  // tic=timstr.tv_sec+(timstr.tv_usec/1000000.0);


  mpz_t n, e, m, c;

  mpz_init(c);
  mpz_init(n);
  mpz_init(e);
  mpz_init(m);

  while(3 == gmp_scanf("%ZX\n %ZX\n %ZX\n", n, e, m )){
    /*compute c = m ^ e  %n   with a window of 6*/
    slidingWindow(c, m, e, n, 6);
    gmp_printf("%ZX\n", c);
  }

  mpz_clear(c);
  mpz_clear(m);
  mpz_clear(e);
  mpz_clear(n);

  // gettimeofday(&timstr,NULL);
  // toc=timstr.tv_sec+(timstr.tv_usec/1000000.0);
  // getrusage(RUSAGE_SELF, &ru);
  // timstr=ru.ru_utime;
  // usrtim=timstr.tv_sec+(timstr.tv_usec/1000000.0);
  // timstr=ru.ru_stime;
  // systim=timstr.tv_sec+(timstr.tv_usec/1000000.0);

  // printf("Elapsed time:\t\t\t%.6f (s)\n", toc-tic);
  // printf("Elapsed user CPU time:\t\t%.6f (s)\n", usrtim);
  // printf("Elapsed system CPU time:\t%.6f (s)\n", systim);
}

/*
Perform stage 2:

- read each 9-tuple of N, d, p, q, d_p, d_q, i_p, i_q and c from stdin,
- compute the RSA decryption m, then
- write the plaintext m to stdout.
*/

void stage2() {
  // struct timeval timstr;        /* structure to hold elapsed time */
  // struct rusage ru;             /* structure to hold CPU time--system and user */
  // double tic,toc;               /* floating point numbers to calculate elapsed wallclock time */
  // double usrtim;                /* floating point number to record elapsed user CPU time */

  // double systim;                /* floating point number to record elapsed system CPU time */
  // gettimeofday(&timstr,NULL);
  // tic=timstr.tv_sec+(timstr.tv_usec/1000000.0);
  mpz_t n;         mpz_init(n);
  mpz_t d;         mpz_init(d);
  mpz_t q;         mpz_init(q);
  mpz_t p;         mpz_init(p);
  mpz_t c;         mpz_init(c);
  mpz_t m;         mpz_init(m);
  mpz_t dp;        mpz_init(dp);
  mpz_t dq;        mpz_init(dq);
  mpz_t ip;        mpz_init(ip);
  mpz_t iq;        mpz_init(iq);
  mpz_t phi_p;     mpz_init(phi_p);
  mpz_t phi_q;     mpz_init(phi_q);

  mpz_t temp;      mpz_init(temp);
  mpz_t tempc;     mpz_init(tempc);
  mp_limb_t w;
  mpz_t ro;        mpz_init(ro);

  while(9 == gmp_scanf("%ZX\n %ZX\n %ZX\n %ZX\n %ZX\n %ZX\n %ZX\n %ZX\n %ZX\n", n, d, p, q, dp, dq, ip, iq, c)){
      mpz_set_ui(m, 0);

      /*compute omega for p*/
      w = calc_omega( p);

      /*compute c % p so that we can feed it to montgomery exponentiation*/
      mont_mod(tempc, c, p, w);

      /*compute temp = c ^ dp  % p */
      slidingWindow(temp, tempc, dp, p, 6);

      /*compute (temp * q * iq) % n  and then add it to m*/
      my_mul(temp, temp, q, n);
      my_mul(temp, temp, iq, n);
      my_add(m, m, temp, n);

      /*compute omega for q*/
      w = calc_omega(q);

      /*compute c % q so that we can feed it to montgomery exponentiation*/
      mont_mod(tempc, c, q, w);

      /*compute temp = c ^ dq  % q */
      slidingWindow(temp, tempc, dq, q, 6);

      /*compute (temp * p * ip) % n and then add it to m */
      my_mul(temp, temp, p, n);
      my_mul(temp, temp, ip, n);
      my_add(m, m, temp, n);

      gmp_printf("%ZX\n", m);
  }

  mpz_clear(n);
  mpz_clear(d);
  mpz_clear(p);
  mpz_clear(q);
  mpz_clear(dp);
  mpz_clear(dq);
  mpz_clear(ip);
  mpz_clear(iq);
  mpz_clear(phi_q);
  mpz_clear(phi_p);
  mpz_clear(c);
  mpz_clear(m);
  mpz_clear(tempc);
  mpz_clear(temp);
  mpz_clear(ro);
  // gettimeofday(&timstr,NULL);
  // toc=timstr.tv_sec+(timstr.tv_usec/1000000.0);
  // getrusage(RUSAGE_SELF, &ru);
  // timstr=ru.ru_utime;
  // usrtim=timstr.tv_sec+(timstr.tv_usec/1000000.0);
  // timstr=ru.ru_stime;
  // systim=timstr.tv_sec+(timstr.tv_usec/1000000.0);

  // printf("Elapsed time:\t\t\t%.6f (s)\n", toc-tic);
  // printf("Elapsed user CPU time:\t\t%.6f (s)\n", usrtim);
  // printf("Elapsed system CPU time:\t%.6f (s)\n", systim);
}

/*
Perform stage 3:

- read each 5-tuple of p, q, g, h and m from stdin,
- compute the ElGamal encryption c = (c_1,c_2), then
- write the ciphertext c to stdout.
*/

void stage3() {
  // struct timeval timstr;        /* structure to hold elapsed time */
  // struct rusage ru;             /* structure to hold CPU time--system and user */
  // double tic,toc;               /* floating point numbers to calculate elapsed wallclock time */
  // double usrtim;                /* floating point number to record elapsed user CPU time */

  // double systim;                /* floating point number to record elapsed system CPU time */
  // gettimeofday(&timstr,NULL);
  // tic=timstr.tv_sec+(timstr.tv_usec/1000000.0);

  FILE *fp = fopen("/dev/urandom", "r");
  mpz_t random;     mpz_init(random);   mpz_set_ui(random, 0);
  mpz_t randKey;    mpz_init(randKey);
  gmp_randstate_t state;
  gmp_randinit_mt(state);
  int c;

  mpz_t p;    mpz_init(p);
  mpz_t q;    mpz_init(q);
  mpz_t g;    mpz_init(g);
  mpz_t h;    mpz_init(h);
  mpz_t m;    mpz_init(m);
  mpz_t r;    mpz_init(r);
  mpz_t c1;   mpz_init(c1);
  mpz_t c2;   mpz_init(c2);
  mpz_t h1;   mpz_init(h1);

  mpz_t ro;  mpz_init(ro);

  mp_limb_t w;
  while( 5 == gmp_scanf("%ZX\n %ZX\n %ZX\n %ZX\n %ZX\n", p, q, g, h, m)){

    /*  generate key by reading 20 characters from dev/urandom
        resulting in a 160 bits random number
    */
    for(int i=0; i<20; i++){
      c = fgetc(fp);
      mpz_mul_2exp(random, random, 8);
      mpz_add_ui(random, random, c);
    }
    /*    generat a random key by feeding the previously obtained random number into the
          gmp seeding function */
    gmp_randseed(state, random);
    mpz_urandomb(randKey, state, 160);

    mpz_set(r, randKey);

    // compute omega for q
    w = calc_omega(q);

    mont_mod(r, r, q, w);

    /*compute c1 = (g ^ r) % p with a window of 6*/
    slidingWindow(c1, g, r, p, 6);
    gmp_printf("%ZX\n", c1);

    /*compute h = (h ^ r) % p with a window of 6*/
    slidingWindow(h, h, r, p, 6);

    /*compute c2 = (m * h) % p */
    my_mul(c2, m, h, p);

    gmp_printf("%ZX\n", c2);
  }

  fclose(fp);
  gmp_randclear(state);
  mpz_clear(randKey);
  mpz_clear(random);
  mpz_clear(p);
  mpz_clear(q);
  mpz_clear(g);
  mpz_clear(h);
  mpz_clear(m);
  mpz_clear(r);
  mpz_clear(c1);
  mpz_clear(c2);
  mpz_clear(h1);
  mpz_clear(ro);

  // gettimeofday(&timstr,NULL);
  // toc=timstr.tv_sec+(timstr.tv_usec/1000000.0);
  // getrusage(RUSAGE_SELF, &ru);
  // timstr=ru.ru_utime;
  // usrtim=timstr.tv_sec+(timstr.tv_usec/1000000.0);
  // timstr=ru.ru_stime;
  // systim=timstr.tv_sec+(timstr.tv_usec/1000000.0);

  // printf("Elapsed time:\t\t\t%.6f (s)\n", toc-tic);
  // printf("Elapsed user CPU time:\t\t%.6f (s)\n", usrtim);
  // printf("Elapsed system CPU time:\t%.6f (s)\n", systim);
}

/*
Perform stage 4:

- read each 5-tuple of p, q, g, x and c = (c_1,c_2) from stdin,
- compute the ElGamal decryption m, then
- write the plaintext m to stdout.
*/

void stage4() {
  // struct timeval timstr;        /* structure to hold elapsed time */
  // struct rusage ru;             /* structure to hold CPU time--system and user */
  // double tic,toc;               /* floating point numbers to calculate elapsed wallclock time */
  // double usrtim;                /* floating point number to record elapsed user CPU time */

  // double systim;                /* floating point number to record elapsed system CPU time */
  // gettimeofday(&timstr,NULL);
  // tic=timstr.tv_sec+(timstr.tv_usec/1000000.0);

  mpz_t p;  mpz_init(p);
  mpz_t q;  mpz_init(q);
  mpz_t g;  mpz_init(g);
  mpz_t c1; mpz_init(c1);
  mpz_t c2; mpz_init(c2);
  mpz_t x;  mpz_init(x);
  mpz_t m;  mpz_init(m);

  mpz_t ro;    mpz_init(ro);
  mp_limb_t w;

  while(6 == gmp_scanf("%ZX\n %ZX\n %ZX\n %ZX\n %ZX\n %ZX\n", p, q, g, x, c1, c2)){

    /*compute omega for q*/
    w = calc_omega(q);

    /*x = x % q*/
    mont_mod(x, x, q, w);

    /*x = -x % q*/
    mpz_sub(x, q, x);

    /*compute m = (c1 ^ x) % p with a window of 6 */
    slidingWindow(m, c1, x, p, 6);

    /*m = (m * c2) % p*/
    my_mul(m, m, c2, p);

    gmp_printf("%ZX\n", m);

  }


  mpz_clear(p);
  mpz_clear(q);
  mpz_clear(g);
  mpz_clear(c1);
  mpz_clear(c2);
  mpz_clear(x);
  mpz_clear(ro);
  mpz_clear(m);

  // gettimeofday(&timstr,NULL);
  // toc=timstr.tv_sec+(timstr.tv_usec/1000000.0);
  // getrusage(RUSAGE_SELF, &ru);
  // timstr=ru.ru_utime;
  // usrtim=timstr.tv_sec+(timstr.tv_usec/1000000.0);
  // timstr=ru.ru_stime;
  // systim=timstr.tv_sec+(timstr.tv_usec/1000000.0);

  // printf("Elapsed time:\t\t\t%.6f (s)\n", toc-tic);
  // printf("Elapsed user CPU time:\t\t%.6f (s)\n", usrtim);
  // printf("Elapsed system CPU time:\t%.6f (s)\n", systim);
}

/*
The main function acts as a driver for the assignment by simply invoking
the correct function for the requested stage.
*/

int main( int argc, char* argv[] ) {
  if( 2 != argc ) {
    abort();
  }

  if     ( !strcmp( argv[ 1 ], "stage1" ) ) {
    stage1();
  }
  else if( !strcmp( argv[ 1 ], "stage2" ) ) {
    stage2();
  }
  else if( !strcmp( argv[ 1 ], "stage3" ) ) {
    stage3();
  }
  else if( !strcmp( argv[ 1 ], "stage4" ) ) {
    stage4();
  }
  else {
    abort();
  }

  return 0;
}




/*tests if mpn addition is correct
  used during development
*/
void test_add(mp_limb_t* x, mp_size_t xn, mp_limb_t* y, mp_size_t yn, int nr){
  mpz_t result, xz, yz, resultLimbs;
  mpz_init(result); mpz_init(resultLimbs);
  mpz_init(xz); mpz_init(yz);

  mpz_import(xz, xn, -1, sizeof(mp_limb_t), 0, 0, x);
  mpz_import(yz, yn, -1, sizeof(mp_limb_t), 0, 0, y);
  mpz_add(result, xz, yz);

  mp_size_t rs;

  mp_limb_t* rl = mpz_limbs_modify(resultLimbs, xn+yn+1);
  mp_limb_t co = 0;
  if(xn > yn)
    co = mpn_add(rl, x, xn, y, yn);
  else
    co = mpn_add(rl, y, yn, x, xn);
  rs = (xn > yn) ? xn : yn;
  if(co == 1){
    rl[rs] = 1;
    rs =rs +1;
  }
    mpz_limbs_finish(resultLimbs, rs);

  if(mpz_cmp(result, resultLimbs) != 0){
    gmp_printf("\nWRONG AAAAAAADDDDDD %d\n", nr);
    gmp_printf("size of correct is %llu, size of wrong is %llu\n", mpz_size(result), mpz_size(resultLimbs));
    gmp_printf("res is %Zd\n", result);
    gmp_printf("res is %Zd\n", resultLimbs);
    print_limbs(result);
    print_limbs(resultLimbs);
  }
}
