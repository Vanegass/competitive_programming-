#include <bits/stdc++.h>
using namespace std;

using ll = long long;

//sin cambios en el arreglo, arreglo estatico
struct SparseTable {
    vector<vector<int>> st;

    SparseTable(const vector<int>& a) {
        int n = (int)a.size();
        if (n == 0) return;
        int lg = __lg(n) + 1;
        st.assign(lg, vector<int>(n));
        st[0] = a;
        for (int j = 1; j < lg; j++)
            for (int i = 0; i + (1 << j) <= n; i++)
                st[j][i] = min(st[j - 1][i], st[j - 1][i + (1 << (j - 1))]);
    }

    int query(int l, int r) const { // [l, r] 0-indexed
        int k = __lg(r - l + 1);
        return min(st[k][l], st[k][r - (1 << k) + 1]);
    }
};

//segment tree iterativo
struct SegTree {
    int n;
    vector<ll> t;
    const ll ID = 0; // Neutro: 0 (suma), LLONG_MAX (min), LLONG_MIN (max)

    ll merge(ll a, ll b) const { return a + b; }

    SegTree(int n) : n(n), t(2 * n, ID) {}

    SegTree(const vector<ll>& a) : n((int)a.size()), t(2 * n, ID) {
        for (int i = 0; i < n; i++) t[i + n] = a[i];
        for (int i = n - 1; i > 0; i--) t[i] = merge(t[i << 1], t[i << 1 | 1]);
    }

    // a[p] = val (0-indexed). Para sumar en un punto, cambia "=" por "+="
    void update(int p, ll val) {
        for (t[p += n] = val; p >>= 1; )
            t[p] = merge(t[p << 1], t[p << 1 | 1]);
    }

    // Consulta en rango inclusivo [l, r] (0-indexed)
    ll query(int l, int r) const {
        ll resL = ID, resR = ID;
        for (l += n, r += n + 1; l < r; l >>= 1, r >>= 1) {
            if (l & 1) resL = merge(resL, t[l++]);
            if (r & 1) resR = merge(t[--r], resR);
        }
        return merge(resL, resR);
    }
};

//segment tree recursivo + lazy propagation
struct LazySegTree {
    int n;
    vector<ll> t, lazy;
    const ll ID = 0; // Neutro de la query: 0 (suma), LLONG_MAX (min), LLONG_MIN (max)

    ll merge(ll a, ll b) const { return a + b; }

    // ÚNICO lugar donde ocurre la matemática de la actualización
    void apply(int v, int l, int r, ll val) {
        t[v] += val * (r - l + 1); // Si la query es min/max, quita "* (r - l + 1)"
        lazy[v] += val;
    }

    void push(int v, int l, int r) {
        if (lazy[v] == 0) return;
        int m = (l + r) >> 1;
        apply(v << 1, l, m, lazy[v]);
        apply(v << 1 | 1, m + 1, r, lazy[v]);
        lazy[v] = 0;
    }

    void build(const vector<ll>& a, int v, int l, int r) {
        if (l == r) {
            t[v] = a[l];
            return;
        }
        int m = (l + r) >> 1;
        build(a, v << 1, l, m);
        build(a, v << 1 | 1, m + 1, r);
        t[v] = merge(t[v << 1], t[v << 1 | 1]);
    }

    void update(int v, int l, int r, int ql, int qr, ll val) {
        if (ql > r || qr < l) return;
        if (ql <= l && r <= qr) {
            apply(v, l, r, val);
            return;
        }
        push(v, l, r);
        int m = (l + r) >> 1;
        update(v << 1, l, m, ql, qr, val);
        update(v << 1 | 1, m + 1, r, ql, qr, val);
        t[v] = merge(t[v << 1], t[v << 1 | 1]);
    }

    ll query(int v, int l, int r, int ql, int qr) {
        if (ql > r || qr < l) return ID;
        if (ql <= l && r <= qr) return t[v];
        push(v, l, r);
        int m = (l + r) >> 1;
        return merge(query(v << 1, l, m, ql, qr),
                     query(v << 1 | 1, m + 1, r, ql, qr));
    }

    LazySegTree(int n) : n(n), t(4 * n, ID), lazy(4 * n, 0) {}

    LazySegTree(const vector<ll>& a) : n((int)a.size()), t(4 * n, ID), lazy(4 * n, 0) {
        if (n > 0) build(a, 1, 0, n - 1);
    }

    void update(int l, int r, ll val) {
        if (l <= r) update(1, 0, n - 1, l, r, val);
    }

    ll query(int l, int r) {
        return l <= r ? query(1, 0, n - 1, l, r) : ID;
    }
};