#include <bits/stdc++.h>
using namespace std;

using ll = long long;
const double PI = acos(-1.0);

struct cd {
    double x, y;
    cd(double x = 0, double y = 0) : x(x), y(y) {}
    cd operator+(const cd& o) const { return {x + o.x, y + o.y}; }
    cd operator-(const cd& o) const { return {x - o.x, y - o.y}; }
    cd operator*(const cd& o) const { return {x * o.x - y * o.y, x * o.y + y * o.x}; }
    void operator/=(double d) { x /= d; y /= d; }
};

void fft(vector<cd>& a, bool invert) {
    int n = (int)a.size();

    // Bit-reversal permutation
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) swap(a[i], a[j]);
    }

    // Raíces precalculadas por nivel (evita error acumulativo de w *= wlen)
    vector<cd> w(n / 2);
    for (int len = 2; len <= n; len <<= 1) {
        int half = len >> 1;
        double ang = 2 * PI / len * (invert ? -1 : 1);
        for (int j = 0; j < half; j++) {
            w[j] = cd(cos(ang * j), sin(ang * j));
        }
        for (int i = 0; i < n; i += len) {
            for (int j = 0; j < half; j++) {
                cd u = a[i + j];
                cd v = a[i + j + half] * w[j];
                a[i + j] = u + v;
                a[i + j + half] = u - v;
            }
        }
    }

    if (invert) {
        for (cd& c : a) c /= n;
    }
}

vector<ll> multiply(const vector<ll>& a, const vector<ll>& b) {
    if (a.empty() || b.empty()) return {};

    int need = (int)a.size() + (int)b.size() - 1;
    int n = 1;
    while (n < need) n <<= 1;

    vector<cd> fa(n), fb(n);
    for (int i = 0; i < (int)a.size(); i++) fa[i] = cd((double)a[i]);
    for (int i = 0; i < (int)b.size(); i++) fb[i] = cd((double)b[i]);

    fft(fa, false);
    fft(fb, false);
    for (int i = 0; i < n; i++) fa[i] = fa[i] * fb[i];
    fft(fa, true);

    vector<ll> res(need);
    for (int i = 0; i < need; i++) {
        res[i] = llround(fa[i].x);
    }

    // Descomentar SOLO si necesitas forma canónica de polinomio (sin ceros finales):
    // while (res.size() > 1 && res.back() == 0) res.pop_back();

    return res;
}