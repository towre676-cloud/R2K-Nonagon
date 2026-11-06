#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stddef.h>
/* Only exported FLINT functions touch opaque storage; no FLINT layout is assumed.
 * 256 bytes is deliberately oversized for scalar/matrix handles on this runtime.
 * The optional Python implementation supplies an independent fallback/control. */
typedef union { max_align_t align; unsigned char data[256]; } O;
typedef void *V;
extern void arb_init(V),arb_clear(V),arb_set(V,V),arb_zero(V),arb_one(V);
extern int arb_set_str(V,const char*,long);
extern void arb_set_si(V,long),arb_add(V,V,V,long),arb_sub(V,V,V,long),arb_mul(V,V,V,long),arb_div(V,V,V,long),arb_mul_si(V,V,long,long),arb_div_ui(V,V,unsigned long,long),arb_sqrt(V,V,long),arb_sinc_pi(V,V,long);
extern int arb_is_positive(V),arb_is_finite(V);
extern char *arb_get_str(V,long,unsigned long);
extern void flint_free(V);
extern void arb_mat_init(V,long,long),arb_mat_clear(V),arb_mat_mul(V,V,V,long),arb_mat_transpose(V,V);
extern V arb_mat_entry_ptr(V,long,long);
static long na,nb,nc,bits;static int ready=0;
static O la,lb,KA,KB,AB,AC,BC,CC,Pm,Qm,Tm,ABt,ABC;
static O *ax,*aw,*bx,*bw,*cx,*cw,*cr,*cqw;
static V ent(O *a,long i,long j){return arb_mat_entry_ptr(a,i,j);}
static O *vec(long n){O *v=calloc(n,sizeof(O));for(long i=0;i<n;i++)arb_init(v+i);return v;}
static void unvec(O *v,long n){for(long i=0;i<n;i++)arb_clear(v+i);free(v);}
static void norm(O *out,O *a,long n,long m){O t;arb_init(&t);arb_zero(out);for(long i=0;i<n;i++)for(long j=0;j<m;j++){arb_mul(&t,ent(a,i,j),ent(a,i,j),bits);arb_add(out,out,&t,bits);}arb_clear(&t);}
static void inner(O *out,O *a,O *b,long n,long m){O t;arb_init(&t);arb_zero(out);for(long i=0;i<n;i++)for(long j=0;j<m;j++){arb_mul(&t,ent(a,i,j),ent(b,i,j),bits);arb_add(out,out,&t,bits);}arb_clear(&t);}
static void kernel(O *out,O *x,O *w,long n,O *y,O *v,long m,O *shift){O t,u;arb_init(&t);arb_init(&u);for(long i=0;i<n;i++)for(long j=0;j<m;j++){arb_add(&t,y+j,shift,bits);arb_sub(&t,&t,x+i,bits);arb_sinc_pi(&t,&t,bits);arb_mul(&u,w+i,v+j,bits);arb_mul(ent(out,i,j),&t,&u,bits);}arb_clear(&t);arb_clear(&u);}
static void points(O *x,O *w,long n,const char **roots,const char **weights,O *length){O t;arb_init(&t);for(long i=0;i<n;i++){arb_set_str(&t,roots[i],bits);arb_mul(x+i,&t,length,bits);arb_div_ui(x+i,x+i,2,bits);arb_set_str(&t,weights[i],bits);arb_mul(w+i,&t,length,bits);arb_div_ui(w+i,w+i,2,bits);arb_sqrt(w+i,w+i,bits);}arb_clear(&t);}
static void cleanup(void){if(!ready)return;O *m[]={&KA,&KB,&AB,&AC,&BC,&CC,&Pm,&Qm,&Tm,&ABt,&ABC};for(int i=0;i<11;i++)arb_mat_clear(m[i]);unvec(ax,na);unvec(aw,na);unvec(bx,nb);unvec(bw,nb);unvec(cx,nc);unvec(cw,nc);unvec(cr,nc);unvec(cqw,nc);arb_clear(&la);arb_clear(&lb);ready=0;}
int counts_config(long a,long b,long c,long prec,const char *l,const char *m,const char **ra,const char **wa,const char **rb,const char **wb,const char **rc,const char **wc){cleanup();na=a;nb=b;nc=c;bits=prec;arb_init(&la);arb_init(&lb);arb_set_str(&la,l,bits);arb_set_str(&lb,m,bits);ax=vec(na);aw=vec(na);bx=vec(nb);bw=vec(nb);cx=vec(nc);cw=vec(nc);cr=vec(nc);cqw=vec(nc);points(ax,aw,na,ra,wa,&la);points(bx,bw,nb,rb,wb,&lb);for(long i=0;i<nc;i++){arb_set_str(cr+i,rc[i],bits);arb_set_str(cqw+i,wc[i],bits);}arb_mat_init(&KA,na,na);arb_mat_init(&KB,nb,nb);arb_mat_init(&AB,na,nb);arb_mat_init(&AC,na,nc);arb_mat_init(&BC,nb,nc);arb_mat_init(&CC,nc,nc);arb_mat_init(&Pm,na,nb);arb_mat_init(&Qm,na,nb);arb_mat_init(&Tm,na,na);arb_mat_init(&ABt,nb,na);arb_mat_init(&ABC,na,nc);O zero;arb_init(&zero);arb_zero(&zero);kernel(&KA,ax,aw,na,ax,aw,na,&zero);kernel(&KB,bx,bw,nb,bx,bw,nb,&zero);arb_clear(&zero);ready=1;return 0;}
char *counts_F(const char *hs,int overlap){if(!ready)return NULL;O h,lo,hi,c,center,half,tmp,zero;O ab,aab,abb,aabb,abab,ac,bc,cc,abc,out,v;
O *z[]={&h,&lo,&hi,&c,&center,&half,&tmp,&zero,&ab,&aab,&abb,&aabb,&abab,&ac,&bc,&cc,&abc,&out,&v};for(int i=0;i<19;i++)arb_init(z[i]);arb_zero(&zero);arb_set_str(&h,hs,bits);
kernel(&AB,ax,aw,na,bx,bw,nb,&h);arb_mat_mul(&Pm,&KA,&AB,bits);arb_mat_mul(&Qm,&AB,&KB,bits);arb_mat_transpose(&ABt,&AB);arb_mat_mul(&Tm,&AB,&ABt,bits);norm(&ab,&AB,na,nb);inner(&aab,&Pm,&AB,na,nb);inner(&abb,&Qm,&AB,na,nb);inner(&aabb,&Pm,&Qm,na,nb);norm(&abab,&Tm,na,na);
/* -AB +2 AAB+2 ABB -4 AABB -2 ABAB +2 AB^2 */
arb_zero(&out);arb_sub(&out,&out,&ab,bits);arb_mul_si(&tmp,&aab,2,bits);arb_add(&out,&out,&tmp,bits);arb_mul_si(&tmp,&abb,2,bits);arb_add(&out,&out,&tmp,bits);arb_mul_si(&tmp,&aabb,4,bits);arb_sub(&out,&out,&tmp,bits);arb_mul_si(&tmp,&abab,2,bits);arb_sub(&out,&out,&tmp,bits);
if(overlap){/* Python passes the analytically known overlap regime. */
/* lo=max(-L/2,h-M/2), hi=min(L/2,h+M/2), selected by positive differences. */
arb_div_ui(&lo,&la,2,bits);arb_mul_si(&lo,&lo,-1,bits);arb_div_ui(&tmp,&lb,2,bits);arb_sub(&v,&h,&tmp,bits);arb_sub(&tmp,&v,&lo,bits);if(arb_is_positive(&tmp))arb_set(&lo,&v);
arb_div_ui(&hi,&la,2,bits);arb_div_ui(&tmp,&lb,2,bits);arb_add(&v,&h,&tmp,bits);arb_sub(&tmp,&hi,&v,bits);if(arb_is_positive(&tmp))arb_set(&hi,&v);
arb_sub(&c,&hi,&lo,bits);if(!arb_is_positive(&c)){for(int i=0;i<19;i++)arb_clear(z[i]);return NULL;}arb_add(&center,&lo,&hi,bits);arb_div_ui(&center,&center,2,bits);arb_div_ui(&half,&c,2,bits);
for(long i=0;i<nc;i++){arb_mul(cx+i,cr+i,&half,bits);arb_add(cx+i,cx+i,&center,bits);arb_mul(cw+i,cqw+i,&half,bits);arb_sqrt(cw+i,cw+i,bits);}kernel(&AC,ax,aw,na,cx,cw,nc,&zero);/* B nodes must be translated, instead shift C by -h. */arb_mul_si(&tmp,&h,-1,bits);kernel(&BC,bx,bw,nb,cx,cw,nc,&tmp);kernel(&CC,cx,cw,nc,cx,cw,nc,&zero);norm(&ac,&AC,na,nc);norm(&bc,&BC,nb,nc);norm(&cc,&CC,nc,nc);arb_mat_mul(&ABC,&AB,&BC,bits);inner(&abc,&ABC,&AC,na,nc);
arb_add(&out,&out,&c,bits);arb_mul_si(&tmp,&ac,2,bits);arb_sub(&out,&out,&tmp,bits);arb_mul_si(&tmp,&bc,2,bits);arb_sub(&out,&out,&tmp,bits);arb_mul_si(&tmp,&cc,2,bits);arb_sub(&out,&out,&tmp,bits);arb_mul_si(&tmp,&abc,8,bits);arb_add(&out,&out,&tmp,bits);arb_sub(&v,&c,&ab,bits);
}else{arb_set(&v,&ab);}arb_mul(&v,&v,&v,bits);arb_mul_si(&v,&v,2,bits);arb_add(&out,&out,&v,bits);
char *result=arb_is_finite(&out)?arb_get_str(&out,55,0):NULL;for(int i=0;i<19;i++)arb_clear(z[i]);return result;}
void counts_release(char *s){flint_free(s);}
void counts_cleanup(void){cleanup();}
