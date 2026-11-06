/* Opaque, overallocated scalar storage; all arithmetic is through FLINT's
   exported public scalar functions. No FLINT fields are accessed. */
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
typedef struct { long align; unsigned char storage[248]; } scalar;
extern void arb_init(void*), arb_clear(void*), arb_set_ui(void*,unsigned long), arb_log(void*,const void*,long), arb_mul(void*,const void*,const void*,long);
extern int arb_set_str(void*,const char*,long);
extern char* arb_get_str(const void*,long,unsigned long);
extern void flint_free(void*);
extern void acb_init(void*),acb_clear(void*),acb_one(void*),acb_set_arb_arb(void*,const void*,const void*),acb_get_real(void*,const void*),acb_get_imag(void*,const void*);
extern void acb_mul_arb(void*,const void*,const void*,long),acb_exp(void*,const void*,long),acb_sub_ui(void*,const void*,unsigned long,long),acb_mul(void*,const void*,const void*,long),acb_div_arb(void*,const void*,const void*,long),acb_neg(void*,const void*),acb_add_ui(void*,const void*,unsigned long,long),acb_conj(void*,const void*),acb_mul_ui(void*,const void*,unsigned long,long),acb_div(void*,const void*,const void*,long),acb_add(void*,const void*,const void*,long);
static scalar *logs_=NULL,*squares_=NULL,*denoms_=NULL;static uint64_t *primes_=NULL;static size_t count_=0;static long bits_=0;
int init_euler(const uint64_t *p,size_t n,long bits){
 if(logs_)return -1;count_=n;bits_=bits;
 logs_=calloc(n,sizeof(scalar));squares_=calloc(n,sizeof(scalar));denoms_=calloc(n,sizeof(scalar));primes_=malloc(n*sizeof(uint64_t));
 if(!logs_||!squares_||!denoms_||!primes_)return -2;
 memcpy(primes_,p,n*sizeof(uint64_t));scalar v;arb_init(&v);
 for(size_t i=0;i<n;i++){
  arb_init(logs_+i);arb_init(squares_+i);arb_init(denoms_+i);
  arb_set_ui(&v,p[i]);arb_log(logs_+i,&v,bits);arb_mul(squares_+i,logs_+i,logs_+i,bits);
  arb_set_ui(&v,p[i]-1);arb_mul(denoms_+i,&v,&v,bits);
 }arb_clear(&v);return 0;
}
int get_ab(const char*r,char**out){
 scalar rr,zero,z,q,t,u,A,B,tmp,aa,bb;arb_init(&rr);arb_init(&zero);
 if(arb_set_str(&rr,r,bits_)){arb_clear(&rr);arb_clear(&zero);return -1;}
 acb_init(&z);acb_init(&q);acb_init(&t);acb_init(&u);acb_init(&A);acb_init(&B);acb_init(&tmp);arb_init(&aa);arb_init(&bb);
 acb_set_arb_arb(&z,&zero,&rr);acb_one(&A);
 for(size_t i=0;i<count_;i++){
  acb_mul_arb(&t,&z,logs_+i,bits_);acb_neg(&t,&t);acb_exp(&q,&t,bits_);
  acb_sub_ui(&t,&q,1,bits_);acb_mul(&t,&t,&t,bits_);acb_div_arb(&t,&t,denoms_+i,bits_);acb_neg(&t,&t);acb_add_ui(&t,&t,1,bits_);acb_mul(&A,&A,&t,bits_);
  acb_conj(&u,&q);acb_mul_ui(&u,&u,primes_[i],bits_);acb_sub_ui(&u,&u,1,bits_);acb_mul(&u,&u,&u,bits_);
  acb_set_arb_arb(&tmp,squares_+i,&zero);acb_div(&t,&tmp,&u,bits_);acb_add(&B,&B,&t,bits_);
 }
 acb_get_real(&aa,&A);acb_get_imag(&bb,&A);out[0]=arb_get_str(&aa,40,0);out[1]=arb_get_str(&bb,40,0);
 acb_get_real(&aa,&B);acb_get_imag(&bb,&B);out[2]=arb_get_str(&aa,40,0);out[3]=arb_get_str(&bb,40,0);
 arb_clear(&rr);arb_clear(&zero);arb_clear(&aa);arb_clear(&bb);acb_clear(&z);acb_clear(&q);acb_clear(&t);acb_clear(&u);acb_clear(&A);acb_clear(&B);acb_clear(&tmp);return 0;
}
void release_string(void*p){flint_free(p);}
int tail_sums(const uint64_t*p,size_t n,char**out){
 scalar s0,s2,v,lp,t,d;arb_init(&s0);arb_init(&s2);arb_init(&v);arb_init(&lp);arb_init(&t);arb_init(&d);
 for(size_t i=0;i<n;i++){
  arb_set_ui(&v,p[i]);arb_log(&lp,&v,128);arb_mul(&lp,&lp,&lp,128);
  arb_set_ui(&v,p[i]-1);arb_mul(&d,&v,&v,128);
  extern void arb_inv(void*,const void*,long),arb_add(void*,const void*,const void*,long),arb_div(void*,const void*,const void*,long);
  arb_inv(&t,&d,128);arb_add(&s0,&s0,&t,128);arb_div(&t,&lp,&d,128);arb_add(&s2,&s2,&t,128);
 }
 out[0]=arb_get_str(&s0,40,0);out[1]=arb_get_str(&s2,40,0);
 arb_clear(&s0);arb_clear(&s2);arb_clear(&v);arb_clear(&lp);arb_clear(&t);arb_clear(&d);return 0;
}
