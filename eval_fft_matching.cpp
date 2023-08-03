#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <math.h>
#include <string.h>
#include <random>
#include <chrono>
#include <map>

#include "fft.hpp"
// #define TO_STRING(VariableName) # VariableName

using namespace std;

string GetFilename(string s, int n, int m) {
    string filename = s;
    filename += to_string(n);
    filename += "_";
    filename += to_string(m);

    return filename;
}

vector<int> ToIntVector(string s) {
    int n = s.length();
    vector<int> ans(n);
    for (int i = 0; i < n; ++i)
        ans[i] = (int)(s[i] - 'a');
    return ans;
}

// template<class T>
using matching_function = vector<size_t> (*)(size_t, const vector<int> &, const vector<int> &, fft_function);

// template<class T>
void EvalMathing(string input_filename, string output_filename, size_t alphabet_size, int N_, int N, fft_function Fft, matching_function Matching) {
    for (int n = N_; n <= N; n *= 10) {
        for (int m = 1; m <= n; m *= 10) {
            string input_filename2 = GetFilename(input_filename, n, m);
            string output_filename2 = GetFilename(output_filename, n, m);
            ifstream ifs(input_filename2);
            ofstream ofs(output_filename2);
            if (!ifs) {
                cout << "cannot open " << input_filename2 << endl;
                exit(1);
            }
            if (!ofs) {
                cout << "cannot open" << output_filename2 << endl;
                exit(1);
            }

            ofs << "n: " << n << "m: " << m << endl << flush;

            string s, t;
            ifs >> n >> m;
            ifs >> s >> t;
            auto a = ToIntVector(s);
            auto b = ToIntVector(t);

            auto start = std::chrono::system_clock::now(); 
            auto locs = GetMatchingLocations(26, a, b, RecursiveFft);
            auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
            auto dur = end - start;        // 要した時間を計算
            auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
            // 要した時間をミリ秒（1/1000秒）に変換して表示
            ofs << msec << " milli sec \n" << flush;

            ifs.close();
            ofs.close();
        } 
    }

    return;
}

int main(int argc, char *argv[]) {
    int N = atoi(argv[1]);

    string output_filename = "fft_matching_time/";
    string input_filename = "samples/sample_";
    
    EvalMathing(input_filename, output_filename+"recursive/not_fast/", 26, 1000000, N, RecursiveFft, GetMatchingLocations);
    // EvalMathing(input_filename, output_filename+"recursive/fast/", 26, 10, N, RecursiveFft, GetMatchingLocationsFast);
    // EvalMathing(input_filename, output_filename+"cooley_tukey/not_fast/", 26, 10, N, CooleyTukeyFft, GetMatchingLocations);
    // EvalMathing(input_filename, output_filename+"cooley_tukey/fast/", 26, 10, N, CooleyTukeyFft, GetMatchingLocationsFast);
    // EvalMathing(input_filename, output_filename+"plain/not_fast/", 26, 10, N, PlainFft, GetMatchingLocations);
    // EvalMathing(input_filename, output_filename+"plain/fast/", 26, 10, N, PlainFft, GetMatchingLocationsFast);
    return 0;
}