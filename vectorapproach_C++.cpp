#include <iostream>
#include <vector>
#include <math.h> 
#include <chrono>
#include <sys/time.h>
using namespace std;

int main()
{
   // Eigen::MatrixXf A = Eigen::MatrixXf::Random(3, 2);
   // std::cout << "A: " << A << std::endl;
    // Start measuring time
    // Start measuring time
    struct timeval begin, end;
    gettimeofday(&begin, 0);
    // auto begin = std::chrono::high_resolution_clock::now();
   // Eigen::JacobiSVD<Eigen::MatrixXf> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
   // std::cout << "U: " << svd.matrixU() << std::endl;;
   // std::cout << "V: " << svd.matrixV() << std::endl;;
   // std::cout << "Sigma: " << svd.singularValues() << std::endl;;
   int M=1600;
   int N=40001;
   int K= 2755;


   vector<int> matched_filter(420, 10);

   vector<vector<int> > time_delay_arr_P(N, vector<int>(M, 0));

   vector<vector<int> > slowtime(M, vector<int>(N, 0));

   vector<vector<int> > rangecell(N, vector<int>(M, 0));

   vector<int> Pnu(N*K, 1);

   vector<int> nu(M, 1);

   vector<int> ph_cor_coeffc(N, 1);

   vector<int> result_matrix(420, 1);


   // Eigen::MatrixXd idx_values::Zero(420,N);
   // // idx_values.setRandom(420,N);
   // cout<< idx_values;
   vector<vector<int> > idx_values(420, vector<int>(N, 0));

//    std::complex<double> c;
//    Eigen::MatrixXd R1; R1.setRandom(2,2);
//    Eigen::MatrixXcd C1 = c*R1; // multiply complex*real
//    Eigen::MatrixXcd C2 = c*C1; // complex scalar times complex matrix
//    C1(0,0) = c; // assign complex value.

   for(int m=0;m<M;m++){
         // std::complex<double> ph_cor_coeffc = exp (c * 2 * 3 * 1 * m);
         // cout<<"k"<<time_delay_arr_P(Eigen::all,m);
         for(int r=0;r<N;r++){
            // cout<<Pnu(idx_values(Eigen::all,r));
            // Pnu(idx_values(Eigen::all,r)) += ph_cor_coeffc(r) * matched_filter * nu(m);
            
            // result_matrix+=ph_cor_coeffc(r) * matched_filter * nu(m);
            for(int l=0;l<420;l++){
                int s=idx_values[l][r];
                Pnu[s] += ph_cor_coeffc[r] * matched_filter[l] * nu[m];
            }
            
            // cout<<"Hello";
         }
   }
    gettimeofday(&end, 0);
    long seconds = end.tv_sec - begin.tv_sec;
    long microseconds = end.tv_usec - begin.tv_usec;
    double elapsed = seconds + microseconds*1e-6;

    
    cout<<"Time measured: %.3f seconds.\n"<<elapsed;
//    cout<<"hello";
//        // Stop measuring time and calculate the elapsed time
//     auto end = std::chrono::high_resolution_clock::now();
//     auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

//     cout<<"Time measured: %.3f seconds.\n"<< (elapsed.count() * 1e-9);
}
