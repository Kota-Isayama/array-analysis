#include <iostream>
#include <vector>
#include <complex>
#include <math.h>
#include <string.h>
#include <random>

#include "fft.hpp"

using namespace std;

template<class T>
bool Check(vector<size_t> loc, const vector<T> &a, const vector<T> &b) {
    int cnt = 0, n = a.size(), m = b.size();
    bool flag = true;
    for (int i = 0; i < n-m+1; ++i) {
        bool match = true;
        for (int j = 0; j < m; ++j)
            match = match && (a[i+j] == b[j]);
        if (match) {
            // cout << i << " ";
            if (cnt == loc.size()) {
                flag = false;
            }
            else if (loc[cnt] != i) {
                flag = false;
            }
            else {
                ++cnt;
            }
        }
        else if (cnt < loc.size() && loc[cnt] == i) {
            return flag = false;
        }
    }
    // cout << endl;
    return flag;
}

int main(int argc, char *argv[]) {
    std::mt19937 rand_src(12345);

    int n = atoi(argv[1]);
    int m = atoi(argv[2]);
    if (n < m) {
        cout << "n < m\n";
        exit(1);
    }

    bool flag_plain = true, flag_recursive = true, flag_recursive_fast = true, flag_cooley_tukey = true;
    for (int i = 0; i < 100; ++i) {
        cout << "test " << i << endl;
        string s, t;
        vector<int> a(n), b(m);
#if 0
        for (int j = 0; j < n; ++j) {
            a[j] = j%26;
            s.push_back((char)(a[j] + 'a'));
        }
        for (int j = 0; j < m; ++j) {
            b[j] = j % 26;
            t.push_back((char)(b[j] + 'a'));
        }
#else
        for (int j = 0; j < n; ++j) {
            a[j] = rand_src() % 26;
            s.push_back((char)(a[j] + 'a'));
        }
        for (int j = 0; j < m; ++j) {
            b[j] = rand_src() % 26;
            t.push_back((char)(b[j] + 'a'));
        }
#endif

        cout << "s: " << s << endl;
        cout << "t: " << t << endl;

        // auto loc_plain = GetMatchingLocations(26, a, b, PlainFft);
        // cout << "plain:     ";
        // for (int x: loc_plain)
        //     cout << x << " ";
        // cout << endl;
        // auto loc_recursive = GetMatchingLocations(26, a, b, RecursiveFft);
        // cout << "recursive: ";
        // for (int x: loc_recursive)
        //     cout << x << " ";
        // cout << endl;
        // auto loc_recursive_fast = GetMatchingLocationsFast(26, a, b, RecursiveFft);
        // bool flag_plain_ = Check(loc_plain, a, b);
        // bool flag_recursive_ = Check(loc_recursive, a, b);
        // bool flag_recursive_fast_ = Check(loc_recursive_fast, a, b);
        auto loc = GetMatchingLocations(26, a, b, CooleyTukeyFft);
        bool flag_cooley_tukey_ = Check(loc, a, b);

        // cout << "plain:         " << ((flag_plain_) ? "OK" : "NG") << endl;
        // cout << "recursive:     " << ((flag_recursive_) ? "OK" : "NG") << endl;
        // cout << "recursive_fast:" << ((flag_recursive_fast_) ? "OK" : "NG") << endl;
        cout << "cooley:    " << ((flag_cooley_tukey_) ? "OK" : "NG") << endl;
        // flag_plain &= flag_plain_;
        // flag_recursive &= flag_recursive_;
        // flag_recursive_fast &= flag_recursive_fast_;

        cout << endl;
    }

    // cout << (flag_plain ? "plain All OK" : "NG") << endl;
    // cout << (flag_recursive ? "recursive All OK" : "NG") << endl;
    // cout << (flag_recursive_fast ? "recursive_fast All OK" : "NG") << endl;
    cout << "All cooley:    " << ((flag_cooley_tukey) ? "All OK" : "NG") << endl;
    return 0;
}