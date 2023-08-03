#ifndef _FFT_HPP
#define _FFT_HPP

#include <iostream>
#include <vector>
#include <complex>
#include <math.h>
#include <string.h>

using namespace std;

using fft_function = vector<complex<double>> (*)(const vector<complex<double>>&, bool);

void PrintCompVec(const vector<complex<double>> &v) {
    for (auto x : v)
        cout << x << " ";
    cout << endl;
}



template<class T>
size_t Smallest2PowGreaterThan(T n) {
    size_t N = 1;
    while (N < n)
        N *= 2;
    return N;
}

vector<complex<double>> PlainDft(const vector<complex<double>> &a, bool inverse = false) {
    size_t n = a.size();
    vector<complex<double>> b(n, complex<double>(0, 0));
    double circle = (inverse ? -2 * M_PI: 2 * M_PI);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            b[i] += a[j] * exp(complex<double>(0, circle * i * j / n));
        }
        if (inverse)
            b[i] /= complex<double>(n, 0);
    }
    return b;
}

vector<complex<double>> RecursiveFft(const vector<complex<double>> &a, bool inverse = false) {
    size_t n = a.size();
    if (n == 1) {
        return a;
    }
    else {
        vector<complex<double>> even(n/2), odd(n/2);
        for (int i = 0; i < n/2; ++i) {
            even[i] = a[2*i];
            odd[i] = a[2*i+1];
        }
        even = RecursiveFft(even, inverse);
        odd = RecursiveFft(odd, inverse);
        vector<complex<double>> ans(n);
        for (int i = 0; i < n/2; ++i) {
            auto omega = (inverse ? exp(complex<double>(0, -2.0*M_PI*i/n)) : exp(complex<double>(0, 2.0*M_PI*i/n)));
            if (inverse) {
                auto comp_2 = complex<double>(2, 0);
                ans[i] = (even[i] + omega * odd[i]) / comp_2;
                ans[i+n/2] = (even[i] - omega * odd[i]) / comp_2;
            }
            else {
                ans[i] = even[i] + omega * odd[i];
                ans[i+n/2] = even[i] - omega * odd[i];
            }
        }
        return ans;
    }
}

void ReplaceForBatterfly(vector<complex<double>> &a) {
    size_t n = a.size();
    int h = 0;
    while (1<<h < n)
        ++h;

    for (int i = 0; i < n; ++i) {
        int j = 0;
        for (int k = 0; k < h; ++k)
            j |= (i>>k & 1) << (h-1-k);
        if (i < j)
            swap(a[i], a[j]);
    }
}

vector<complex<double>> CooleyTukeyFft(const vector<complex<double>> &a, bool inverse = false) {
    size_t n = a.size();
    auto ans = a;
    ReplaceForBatterfly(ans);
    // cout << "a" << endl;
    // PrintCompVec(a);
    // cout << "ans" << endl;
    // PrintCompVec(ans);

    for (int b = 1; b < n; b *= 2) {
        for (int j = 0; j < b; j++) {
            double circle = (inverse ? -2*M_PI: 2*M_PI);
            complex<double> w = exp(complex<double>(0, circle / (2.0 * b) * j));
            for (int k = 0; k < n; k += b * 2) {
                complex<double> s = ans[j + k];        
                complex<double> t = ans[j + k + b] * w;
                ans[j + k] = s + t;                
                ans[j + k + b] = s - t;                
            }
        }
    }
    if (inverse) {
        for (int i = 0; i < n; i++)
            ans[i] /= n;
    }
    return ans;
}

// template<class T>
vector<int> Convolution(const vector<int> &a, const vector<int> &b, fft_function Fft) {
    if (a.size() != b.size()) {
        cout << "convolution: don't match length\n";
        exit(1);
    }
    int n = a.size();
    if (n != Smallest2PowGreaterThan(n)) {
        cout << "convolution: not pow of 2\n";
        exit(1);
    }

    vector<complex<double>> a_comp(n), b_comp(n); 
    for (int i = 0; i < n; ++i) {
        a_comp[i] = complex<double>(a[i], 0);
        b_comp[i] = complex<double>(b[i], 0);
    }
    // cout << "a_comp:";
    // PrintCompVec(a_comp);
    // cout << "b_comp:";
    // PrintCompVec(b_comp);
    auto A = Fft(a_comp, false), B = Fft(b_comp, false);
    // cout << "A:";
    // PrintCompVec(A);
    // cout << "B:";
    // PrintCompVec(B);
    vector<complex<double>> C(n);
    for (int i = 0; i < n; ++i)
        C[i] = A[i] * B[i];
    
    auto c = Fft(C, true);
    // cout << "C:";
    // PrintCompVec(c);

    vector<int> ans(n);
    for (size_t i = 0; i < n; ++i) {
        // cout << i << " " << c[i].real() << endl;
        ans[i] = round(c[i].real());
    }
    return ans;
}


template<class T>
vector<T> GetHammingWeight(const vector<T> &s, const vector<T> &t, fft_function Fft) {
    size_t n = s.size(), m = t.size(), N = Smallest2PowGreaterThan(n+m-1);
    vector<T> s_input(N, 0), u_input(N, 0);
    for (int i = 0; i < n; ++i)
        s_input[i] = s[i];
    for (int i = m-1; i >= 0; --i)
        u_input[i] = t[m-1-i];
    
    auto c = Convolution(s_input, u_input, Fft);
    vector<T> ans(n);
    for (int i = 0; i < n; ++i)
        ans[i] = c[m-1+i];
    return ans; 
}

// template<class T>
vector<size_t> GetMatchingLocations(size_t alphabet_size, const vector<int> &s, const vector<int> &t, fft_function Fft) {
    vector<int> hamming_weights(s.size());
    for (int c = 0; c < alphabet_size; ++c) {
        vector<int> s_(s.size()), t_(t.size());
        // cout << "s_:\n";
        for (size_t i = 0; i < s.size(); ++i) {
            s_[i] = static_cast<int>(s[i] == c);
            // cout << s_[i] << " ";
        }
        // cout << endl << "t_:\n";
        for (size_t i = 0; i < t.size(); ++i) {
            t_[i] = static_cast<int>(t[i] == c);
            // cout << t_[i] << " ";
        }
        // cout << endl;

        auto distance = GetHammingWeight(s_, t_, Fft);
        // cout << (char)(c+'a') << endl;
        for (size_t i = 0; i < s.size(); ++i) {
            hamming_weights[i] += distance[i];
            // cout << distance[i] << " ";
        }
        // cout << endl<< endl;
    }

    vector<size_t> ans;
    for (size_t i = 0; i < s.size(); ++i) {
        if (hamming_weights[i] == t.size())
            ans.emplace_back(i);
    }
    return ans;
}

// template<class T>
vector<size_t> GetMatchingLocationsFast(size_t alphabet_size, const vector<int> &s, const vector<int> &t, fft_function Fft) {
    size_t n = s.size(), m = t.size();
    vector<size_t> ans;
    for (size_t k = 0; k < n / (m+1); ++k) {
        vector<int> s_tmp;
        for (int l = k*m + k; l < min(n, (k+2)*m + k); ++l)
            s_tmp.push_back(s[l]);
        auto res_tmp = GetMatchingLocations(alphabet_size, s_tmp, t, Fft);
        for (auto x: res_tmp)
            ans.push_back(x + k*m+k);
    }
    return ans;
}



#endif