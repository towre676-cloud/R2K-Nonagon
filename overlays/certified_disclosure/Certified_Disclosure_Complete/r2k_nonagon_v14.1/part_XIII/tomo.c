#include <time.h>
/* tomo.c -- one-port known-support network tomography, exact modulo a word prime.
 *
 * Input (text, stdin):
 *   p n port m ng
 *   perm[0..n-1]                 fill-reducing elimination order (old indices)
 *   m lines:  a c w              branch (old indices), weight residue mod p
 *   ng lines: i g                ground (old index), residue mod p
 *   nr                           number of extra R=1 sources requested
 *   nr source indices (old)
 * Pipeline: sparse LDL^T of K mod p in the given order; selected inverse Z = K^{-1} on the filled pattern
 * (Takahashi recurrences); port column by triangular solves; full columns for requested extra sources;
 * row system (K g_s)_i = delta_is in the unknowns (branch weights, grounds); sparse incremental echelon form
 * mod p with a heap-ordered dense accumulator; back substitution.
 * Output (stdout): rank U nrows readings, then the solution residues (one per unknown) if full rank, then
 * the list of unknown indices without pivot (if deficient).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef uint64_t u64; typedef int64_t i64;
static u64 P;
static u64 mulm(u64 a, u64 b) { return (a * b) % P; }
static u64 addm(u64 a, u64 b) { u64 c = a + b; return c >= P ? c - P : c; }
static u64 subm(u64 a, u64 b) { return a >= b ? a - b : a + P - b; }
static u64 powm(u64 a, u64 e) { u64 r = 1; a %= P; while (e) { if (e & 1) r = mulm(r, a); a = mulm(a, a); e >>= 1; } return r; }
static u64 invm(u64 a) { return powm(a, P - 2); }
static void *xm(size_t s) { void *q = malloc(s ? s : 1); if (!q) { fprintf(stderr, "oom\n"); exit(2); } return q; }
static void *xr(void *q, size_t s) { q = realloc(q, s ? s : 1); if (!q) { fprintf(stderr, "oom\n"); exit(2); } return q; }
static int cmpint(const void *a, const void *b) { int x = *(const int *)a, y = *(const int *)b; return (x > y) - (x < y); }

int n, port, m, ng;
int *perm, *iperm;             /* perm[new] = old ; iperm[old] = new */
int *ea, *ec; u64 *ew;         /* branches (old) */
int *gi; u64 *gv;              /* grounds (old) */
/* K in new indexing: adjacency lists lower+upper */
int *kdeg, **kadj; u64 **kval; u64 *kdiag;
/* L structure: S[j] sorted rows > j, Lv[j] values; D[j] */
int *Slen, **S; u64 **Lv; u64 *D, *invD; int *parent;
/* Z selected inverse on S */
u64 **Zv; u64 *Zd;

static void add_k(int a, int b, u64 v) { /* add v to K[a][b] (new idx, a != b) */
    for (int t = 0; t < kdeg[a]; t++) if (kadj[a][t] == b) { kval[a][t] = addm(kval[a][t], v); return; }
    kadj[a] = xr(kadj[a], sizeof(int) * (kdeg[a] + 1)); kval[a] = xr(kval[a], sizeof(u64) * (kdeg[a] + 1));
    kadj[a][kdeg[a]] = b; kval[a][kdeg[a]] = v; kdeg[a]++;
}
static int find_pos(int j, int i) { /* position of row i in S[j], or -1 */
    int lo = 0, hi = Slen[j] - 1;
    while (lo <= hi) { int mid = (lo + hi) >> 1; if (S[j][mid] == i) return mid; if (S[j][mid] < i) lo = mid + 1; else hi = mid - 1; }
    return -1;
}
static int zget(int i, int k, u64 *out) { /* Z[i][k], new indices */
    if (i == k) { *out = Zd[i]; return 1; }
    int a = i < k ? i : k, b = i < k ? k : i; int pos = find_pos(a, b);
    if (pos < 0) return 0; *out = Zv[a][pos]; return 1;
}
/* solve K x = e_r (new index r) -> x (new indexing) */
static void solve_col(int r, u64 *x) {
    memset(x, 0, sizeof(u64) * n); x[r] = 1;
    for (int j = 0; j < n; j++) if (x[j]) for (int t = 0; t < Slen[j]; t++) x[S[j][t]] = subm(x[S[j][t]], mulm(Lv[j][t], x[j]));
    for (int j = 0; j < n; j++) x[j] = mulm(x[j], invD[j]);
    for (int j = n - 1; j >= 0; j--) { u64 acc = x[j]; for (int t = 0; t < Slen[j]; t++) acc = subm(acc, mulm(Lv[j][t], x[S[j][t]])); x[j] = acc; }
}
/* ---- sparse echelon ---- */
typedef struct { int len; int *idx; u64 *val; u64 rhs; } Row;
Row *piv; int U; u64 *acc; char *mark; int *heap, hsz;
static void hpush(int v) { int i = hsz++; heap[i] = v; while (i && heap[(i - 1) / 2] > heap[i]) { int t = heap[i]; heap[i] = heap[(i - 1) / 2]; heap[(i - 1) / 2] = t; i = (i - 1) / 2; } }
static int hpop(void) { int r = heap[0]; heap[0] = heap[--hsz]; int i = 0; for (;;) { int l = 2 * i + 1, q = l + 1, s = i; if (l < hsz && heap[l] < heap[s]) s = l; if (q < hsz && heap[q] < heap[s]) s = q; if (s == i) break; int t = heap[i]; heap[i] = heap[s]; heap[s] = t; i = s; } return r; }
static int insert_row(int len, int *idx, u64 *val, u64 rhs) { /* returns 1 if new pivot */
    hsz = 0; u64 r = rhs;
    for (int t = 0; t < len; t++) if (val[t]) { int c = idx[t]; if (!mark[c]) { mark[c] = 1; acc[c] = 0; hpush(c); } acc[c] = addm(acc[c], val[t]); }
    while (hsz) {
        int c = hpop(); mark[c] = 0; u64 a = acc[c]; acc[c] = 0;
        if (!a) continue;
        if (piv[c].len) { /* eliminate */
            Row *R = &piv[c]; /* R normalized: val[0] = 1 at idx[0] = c */
            for (int t = 1; t < R->len; t++) { int cc = R->idx[t]; if (!mark[cc]) { mark[cc] = 1; acc[cc] = 0; hpush(cc); } acc[cc] = subm(acc[cc], mulm(a, R->val[t])); }
            r = subm(r, mulm(a, R->rhs));
        } else { /* new pivot at c: collect remaining */
            int cap = hsz + 1; Row *R = &piv[c]; R->idx = xm(sizeof(int) * cap); R->val = xm(sizeof(u64) * cap); u64 ia = invm(a);
            R->idx[0] = c; R->val[0] = 1; int L = 1;
            /* drain heap in order */
            while (hsz) { int cc = hpop(); mark[cc] = 0; u64 v = acc[cc]; acc[cc] = 0; if (v) { if (L >= cap) { cap *= 2; R->idx = xr(R->idx, sizeof(int) * cap); R->val = xr(R->val, sizeof(u64) * cap); } R->idx[L] = cc; R->val[L] = mulm(v, ia); L++; } }
            R->len = L; R->rhs = mulm(r, ia); return 1;
        }
    }
    if (r) { fprintf(stderr, "inconsistent row\n"); exit(3); }
    return 0;
}
int main(void) {
    if (scanf("%lu %d %d %d %d", &P, &n, &port, &m, &ng) != 5) return 1;
    perm = xm(sizeof(int) * n); iperm = xm(sizeof(int) * n);
    for (int i = 0; i < n; i++) { if (scanf("%d", &perm[i]) != 1) return 1; iperm[perm[i]] = i; }
    ea = xm(sizeof(int) * m); ec = xm(sizeof(int) * m); ew = xm(sizeof(u64) * m);
    for (int e = 0; e < m; e++) if (scanf("%d %d %lu", &ea[e], &ec[e], &ew[e]) != 3) return 1;
    gi = xm(sizeof(int) * (ng + 1)); gv = xm(sizeof(u64) * (ng + 1));
    for (int g = 0; g < ng; g++) if (scanf("%d %lu", &gi[g], &gv[g]) != 2) return 1;
    int nr; if (scanf("%d", &nr) != 1) nr = 0; int *rs = xm(sizeof(int) * (nr + 1)); for (int t = 0; t < nr; t++) if (scanf("%d", &rs[t]) != 1) return 1;
    /* optional: one radius per extra source (default 1). Rows i with 1 <= dist(s,i) <= radius are used. */
    int *rad = xm(sizeof(int) * (nr + 1)); for (int t = 0; t < nr; t++) rad[t] = 1;
    { int nrad; if (scanf("%d", &nrad) == 1 && nrad == nr) for (int t = 0; t < nr; t++) if (scanf("%d", &rad[t]) != 1) return 1; }
    /* K new-indexed */
    kdeg = calloc(n, sizeof(int)); kadj = calloc(n, sizeof(int *)); kval = calloc(n, sizeof(u64 *)); kdiag = calloc(n, sizeof(u64));
    for (int e = 0; e < m; e++) { int a = iperm[ea[e]], c = iperm[ec[e]]; u64 w = ew[e];
        kdiag[a] = addm(kdiag[a], w); kdiag[c] = addm(kdiag[c], w); add_k(a, c, subm(0, w)); add_k(c, a, subm(0, w)); }
    for (int g = 0; g < ng; g++) { int a = iperm[gi[g]]; kdiag[a] = addm(kdiag[a], gv[g]); }
    /* symbolic: S_j = {i>j: K_ij != 0} U  U_{children c} S_c \ {j} */
    Slen = calloc(n, sizeof(int)); S = calloc(n, sizeof(int *)); parent = xm(sizeof(int) * n);
    int **kids = calloc(n, sizeof(int *)); int *nk = calloc(n, sizeof(int)), *ck = calloc(n, sizeof(int));
    int *tmp = xm(sizeof(int) * n); char *seen = calloc(n, 1); long long nnzL = 0;
    for (int j = 0; j < n; j++) {
        int L = 0;
        for (int t = 0; t < kdeg[j]; t++) if (kadj[j][t] > j && !seen[kadj[j][t]]) { seen[kadj[j][t]] = 1; tmp[L++] = kadj[j][t]; }
        for (int q = 0; q < nk[j]; q++) { int c = kids[j][q]; for (int t = 0; t < Slen[c]; t++) { int i = S[c][t]; if (i > j && !seen[i]) { seen[i] = 1; tmp[L++] = i; } } }
        for (int t = 0; t < L; t++) seen[tmp[t]] = 0;
        qsort(tmp, L, sizeof(int), cmpint); S[j] = xm(sizeof(int) * L); memcpy(S[j], tmp, sizeof(int) * L); Slen[j] = L; nnzL += L;
        parent[j] = L ? S[j][0] : -1;
        if (L) { int pj = S[j][0]; if (nk[pj] == ck[pj]) { ck[pj] = ck[pj] ? 2 * ck[pj] : 4; kids[pj] = xr(kids[pj], sizeof(int) * ck[pj]); } kids[pj][nk[pj]++] = j; }
    }
    fprintf(stderr, "nnz(L) = %lld t=%.1fs\n", nnzL, (double)clock()/CLOCKS_PER_SEC);
    /* row lists: for each row j, the columns k<j with j in S_k, and position */
    int *rl = calloc(n, sizeof(int)); for (int k = 0; k < n; k++) for (int t = 0; t < Slen[k]; t++) rl[S[k][t]]++;
    int **rk = xm(sizeof(int *) * n), **rp = xm(sizeof(int *) * n); int *rf = calloc(n, sizeof(int));
    for (int j = 0; j < n; j++) { rk[j] = xm(sizeof(int) * (rl[j] + 1)); rp[j] = xm(sizeof(int) * (rl[j] + 1)); }
    for (int k = 0; k < n; k++) for (int t = 0; t < Slen[k]; t++) { int j = S[k][t]; rk[j][rf[j]] = k; rp[j][rf[j]] = t; rf[j]++; }
    /* numeric left-looking LDL^T */
    Lv = xm(sizeof(u64 *) * n); D = xm(sizeof(u64) * n); u64 *w = calloc(n, sizeof(u64));
    for (int j = 0; j < n; j++) {
        w[j] = kdiag[j]; for (int t = 0; t < Slen[j]; t++) w[S[j][t]] = 0;
        for (int t = 0; t < kdeg[j]; t++) if (kadj[j][t] > j) w[kadj[j][t]] = kval[j][t];
        for (int q = 0; q < rl[j]; q++) { int k = rk[j][q], pos = rp[j][q]; u64 ljk = Lv[k][pos], f = mulm(ljk, D[k]);
            w[j] = subm(w[j], mulm(f, ljk));
            for (int t = pos + 1; t < Slen[k]; t++) { int i = S[k][t]; w[i] = subm(w[i], mulm(Lv[k][t], f)); } }
        D[j] = w[j]; if (!D[j]) { fprintf(stderr, "zero pivot mod p at %d\n", j); return 4; }
        u64 id = invm(D[j]); Lv[j] = xm(sizeof(u64) * (Slen[j] + 1));
        for (int t = 0; t < Slen[j]; t++) { Lv[j][t] = mulm(w[S[j][t]], id); w[S[j][t]] = 0; }
        w[j] = 0;
    }
    invD = xm(sizeof(u64) * n); for (int j = 0; j < n; j++) invD[j] = invm(D[j]);
    fprintf(stderr, "factorized t=%.1fs\n", (double)clock()/CLOCKS_PER_SEC);
    /* selected inverse */
    Zv = xm(sizeof(u64 *) * n); Zd = xm(sizeof(u64) * n);
    { /* Takahashi on the filled pattern. For k = S_j[b], every i in S_j with i > k lies in S_k (elimination-tree closure),
         so Z(i,k) is read by a forward merge over S_k. Products < 2^62 are accumulated in 128 bits and reduced once. */
      int maxS = 1; for (int j = 0; j < n; j++) if (Slen[j] > maxS) maxS = Slen[j];
      unsigned __int128 *y = xm(sizeof(unsigned __int128) * maxS);
      for (int j = n - 1; j >= 0; j--) {
        int Lj = Slen[j]; int *Sj = S[j]; u64 *lj = Lv[j]; Zv[j] = xm(sizeof(u64) * (Lj + 1));
        for (int a = 0; a < Lj; a++) y[a] = 0;
        for (int b = 0; b < Lj; b++) {
            int k = Sj[b]; u64 lb = lj[b]; y[b] += (unsigned __int128)Zd[k] * lb;
            int *Sk = S[k]; u64 *zk = Zv[k]; int Lk = Slen[k], t = 0;
            for (int a = b + 1; a < Lj; a++) { int i = Sj[a]; while (t < Lk && Sk[t] < i) t++;
                if (t >= Lk || Sk[t] != i) { fprintf(stderr, "pattern miss\n"); return 5; }
                u64 z = zk[t]; y[a] += (unsigned __int128)z * lb; y[b] += (unsigned __int128)z * lj[a]; }
        }
        for (int a = 0; a < Lj; a++) { u64 r = (u64)(y[a] % P); Zv[j][a] = r ? P - r : 0; }
        unsigned __int128 acc2 = 0; for (int b = 0; b < Lj; b++) acc2 += (unsigned __int128)lj[b] * Zv[j][b];
        Zd[j] = subm(invD[j], (u64)(acc2 % P));
      }
      free(y);
    }
    fprintf(stderr, "selected inverse done t=%.1fs\n", (double)clock()/CLOCKS_PER_SEC);
    /* port column and requested columns, old indexing accessors */
    u64 *xp = xm(sizeof(u64) * n); solve_col(iperm[port], xp);
    u64 **xc = xm(sizeof(u64 *) * (nr + 1)); for (int t = 0; t < nr; t++) { xc[t] = xm(sizeof(u64) * n); solve_col(iperm[rs[t]], xc[t]); }
    /* adjacency (old) and unknown index: edges 0..m-1, grounds m..m+ng-1 */
    int *deg = calloc(n, sizeof(int)); for (int e = 0; e < m; e++) { deg[ea[e]]++; deg[ec[e]]++; }
    int **inc = xm(sizeof(int *) * n); int *fi = calloc(n, sizeof(int)); for (int i = 0; i < n; i++) inc[i] = xm(sizeof(int) * (deg[i] + 1));
    for (int e = 0; e < m; e++) { inc[ea[e]][fi[ea[e]]++] = e; inc[ec[e]][fi[ec[e]]++] = e; }
    int *gidx = xm(sizeof(int) * n); for (int i = 0; i < n; i++) gidx[i] = -1; for (int g = 0; g < ng; g++) gidx[gi[g]] = m + g;
    U = m + ng; piv = calloc(U, sizeof(Row)); acc = calloc(U, sizeof(u64)); mark = calloc(U, 1); heap = xm(sizeof(int) * (U + 1));
    int *ri = xm(sizeof(int) * 64); u64 *rv = xm(sizeof(u64) * 64); int rank = 0; long nrows = 0;
    /* G accessor: column given by vector (old row index) or via Z */
    #define GCOL(vec, i) (vec[iperm[i]])
    int *ball = xm(sizeof(int) * n), *dist = xm(sizeof(int) * n); for (int i = 0; i < n; i++) dist[i] = -1;
    for (int pass = 0; pass < 3; pass++) {
        int ns = pass == 0 ? 1 : (pass == 1 ? n : nr);
        for (int q = 0; q < ns; q++) {
            int s; u64 *vec = NULL;
            if (pass == 0) { s = port; vec = xp; }
            else if (pass == 1) { s = q; if (s == port) continue; }
            else { s = rs[q]; vec = xc[q]; }
            int rows_from = (pass == 1) ? s : 0, rows_to = (pass == 1) ? s : n - 1;
            int nball = 0;
            if (pass == 2) { /* BFS ball of radius rad[q] around s; rows at distance >= 1 (s itself already in pass 1) */
                nball = 0; ball[nball++] = s; dist[s] = 0;
                for (int h = 0; h < nball; h++) { int u = ball[h]; if (dist[u] >= rad[q]) continue;
                    for (int t = 0; t < deg[u]; t++) { int e = inc[u][t]; int o = ea[e] == u ? ec[e] : ea[e]; if (dist[o] < 0) { dist[o] = dist[u] + 1; ball[nball++] = o; } } }
                rows_from = 1; rows_to = nball - 1;
            }
            for (int ii = rows_from; ii <= rows_to; ii++) {
                int i = pass == 2 ? ball[ii] : ii;
                int L = 0; int need = deg[i] + 1; if (need > 64) { ri = xr(ri, sizeof(int) * need); rv = xr(rv, sizeof(u64) * need); }
                u64 gis; if (vec) gis = GCOL(vec, i); else if (!zget(iperm[i], iperm[s], &gis)) { fprintf(stderr, "miss\n"); return 6; }
                for (int t = 0; t < deg[i]; t++) { int e = inc[i][t]; int o = ea[e] == i ? ec[e] : ea[e]; u64 gos;
                    if (vec) gos = GCOL(vec, o); else if (!zget(iperm[o], iperm[s], &gos)) { fprintf(stderr, "miss2\n"); return 6; }
                    ri[L] = e; rv[L] = subm(gis, gos); L++; }
                if (gidx[i] >= 0) { ri[L] = gidx[i]; rv[L] = gis; L++; }
                rank += insert_row(L, ri, rv, i == s ? 1 : 0); nrows++;
            }
            if (pass == 2) for (int h = 0; h < nball; h++) dist[ball[h]] = -1;
        }
        fprintf(stderr, "pass %d rank %d / %d t=%.1fs\n", pass, rank, U, (double)clock()/CLOCKS_PER_SEC);
    }
    printf("%d %d %ld\n", rank, U, nrows);
    if (rank == U) { /* back substitution */
        u64 *x = calloc(U, sizeof(u64));
        for (int c = U - 1; c >= 0; c--) { Row *R = &piv[c]; u64 v = R->rhs; for (int t = 1; t < R->len; t++) v = subm(v, mulm(R->val[t], x[R->idx[t]])); x[c] = v; }
        for (int c = 0; c < U; c++) printf("%lu\n", x[c]);
    } else { for (int c = 0; c < U; c++) if (!piv[c].len) printf("free %d\n", c);
        /* kernel basis vectors: support of each (free column set to 1, other free columns 0) */
        u64 *x = calloc(U, sizeof(u64));
        for (int f = 0; f < U; f++) if (!piv[f].len) {
            memset(x, 0, sizeof(u64) * U); x[f] = 1;
            for (int c = U - 1; c >= 0; c--) if (piv[c].len) { Row *R = &piv[c]; u64 v = 0; for (int t = 1; t < R->len; t++) v = subm(v, mulm(R->val[t], x[R->idx[t]])); x[c] = v; }
            printf("kern %d", f); int cnt = 0; for (int c = 0; c < U; c++) if (x[c]) { if (cnt < 400) printf(" %d", c); cnt++; } printf(" | %d\n", cnt);
        }
    }
    long long fill = 0; for (int c = 0; c < U; c++) fill += piv[c].len; fprintf(stderr, "echelon nnz %lld\n", fill);
    return 0;
}
