#include <bits/stdc++.h>
using namespace std;

using ll = long long;

const int MOD = 998244353; // 119 * 2^23 + 1
const int G = 3;           // Raíz primitiva módulo MOD

int binpow(int a, ll b) {
    int res = 1;
    while (b > 0) {
        if (b & 1) res = (int)(1LL * res * a % MOD);
        a = (int)(1LL * a * a % MOD);
        b >>= 1;
    }
    return res;
}

void ntt(vector<int>& a, bool invert) {
    int n = (int)a.size();
    if (n == 1) return;

    // Bit-reversal permutation
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) swap(a[i], a[j]);
    }

    // Precalculamos potencias de wlen por nivel (mitad de % MOD en el loop interno)
    vector<int> w(n / 2);
    w[0] = 1;
    for (int len = 2; len <= n; len <<= 1) {
        int half = len >> 1;
        int wlen = binpow(G, (MOD - 1) / len);
        if (invert) wlen = binpow(wlen, MOD - 2);
        for (int j = 1; j < half; j++) {
            w[j] = (int)(1LL * w[j - 1] * wlen % MOD);
        }
        for (int i = 0; i < n; i += len) {
            for (int j = 0; j < half; j++) {
                int u = a[i + j];
                int v = (int)(1LL * a[i + j + half] * w[j] % MOD);
                int s = u + v;
                int d = u - v;
                a[i + j] = (s < MOD ? s : s - MOD);
                a[i + j + half] = (d >= 0 ? d : d + MOD);
            }
        }
    }

    if (invert) {
        int n_inv = binpow(n, MOD - 2);
        for (int& x : a) x = (int)(1LL * x * n_inv % MOD);
    }
}

vector<int> multiply(const vector<int>& a, const vector<int>& b) {
    if (a.empty() || b.empty()) return {};

    int need = (int)a.size() + (int)b.size() - 1;
    int n = 1;
    while (n < need) n <<= 1;

    vector<int> fa(n, 0), fb(n, 0);
    for (int i = 0; i < (int)a.size(); i++) {
        int x = a[i] % MOD;
        fa[i] = (x < 0 ? x + MOD : x);
    }
    for (int i = 0; i < (int)b.size(); i++) {
        int x = b[i] % MOD;
        fb[i] = (x < 0 ? x + MOD : x);
    }

    ntt(fa, false);
    ntt(fb, false);
    for (int i = 0; i < n; i++) {
        fa[i] = (int)(1LL * fa[i] * fb[i] % MOD);
    }
    ntt(fa, true);

    fa.resize(need);
    // Descomentar SOLO si necesitas forma canónica de polinomio (sin ceros finales):
    // while (fa.size() > 1 && fa.back() == 0) fa.pop_back();

    return fa;
}