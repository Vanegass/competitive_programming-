//link problema: https://www.spoj.com/problems/POLYMUL/
//algoritmo: FFT O(nlogn)

#include <bits/stdc++.h>

//mi plantilla - jaja :D
using namespace std;
using ll = long long;
using vi = vector<int>;
using vo = vector<bool>;
using pii = pair<int,int>;
using pll = pair<ll,ll>;
using vpii = vector<pii>;
const int INF = 1e9;

#define pb push_back
#define all(x) (x).begin(),(x).end()
#define allr(x) (x).rbegin(),(x).rend()
#define rep(i,a,b) for(int i = (a); i < (b); ++i)

using cd = complex<double>;

const double PI = acos(-1);

void fft(vector<cd> &a, bool invert) {
    int n = a.size();
    
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * PI / len * (invert ? -1 : 1);
        cd wlen(cos(ang), sin(ang));
        
        for (int i = 0; i < n; i += len) {
            cd w(1);
            for (int j = 0; j < len / 2; j++) {
                cd u = a[i + j];
                cd v = a[i + j + len / 2] * w;
                
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                
                w *= wlen;
            }
        }
    }

    if (invert) {
        for (cd &x : a) {
            x /= n;
        }
    }
}

vector<ll> multiply(vector<ll> const& a, vector<ll> const& b) {
    if (a.empty() || b.empty()) return {};
    
    vector<cd> fa(a.begin(), a.end());
    vector<cd> fb(b.begin(), b.end());
    
    int n = 1;
    while (n < a.size() + b.size()) n <<= 1;
    
    fa.resize(n);
    fb.resize(n);

    fft(fa, false);
    fft(fb, false);
    
    for (int i = 0; i < n; i++) {
        fa[i] *= fb[i];
    }
    
    fft(fa, true);

    vector<ll> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = round(fa[i].real());
    }
    
    while (result.size() > 1 && result.back() == 0) {
        result.pop_back();
    }
    
    return result;
}

int main(){
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int t;
    cin >> t;

    while(t--){
        int n;
        cin >> n;

        vector<ll> a(n + 1), b(n + 1);

        for(int i = n; i >= 0; i--){
            cin >> a[i];
        }

        for(int i = n; i >= 0; i--){
            cin >> b[i];
        }

        vector<ll> ans = multiply(a, b);

        for(int i = ans.size() - 1; i >= 0; i--){
            cout << ans[i] << (i == 0 ? '\n' : ' ');
        }
    }

    return 0;
}