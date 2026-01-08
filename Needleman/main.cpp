#include "needleman_wunsch.hpp"
#include <iostream>
#include <chrono>


using namespace std;


int main () {
    int LENA = 100;
    string SEQ_A = "CGGTATTCTTCGAAAAGCGTAGGCGAAGCCCGACGGCATGCGGACATGGGCGAGAGCCTACTTGAGGGTAGCACTACATAGGTTTGATCGTACTGGGGGA";
    
    int LENB = 100;
    string SEQ_B = "ACCCGGAAGATGCGGCGGAGAAAACTCGGGTATGCCCCAGACAGTTGAGTGGGTTTCGTATATTCATGCGCCGTGCCAACCGGGGAGATGCTAGGCATGA";

    cout << "Partie Needleman-wunsch:" << endl;
    cout << "Sequence A: " << SEQ_A << endl;
    cout << "Sequence B: " << SEQ_B << endl;

	auto t_start = std::chrono::high_resolution_clock::now();
    cout << "Score en utilisant la version sequentielle:\t" << nw_distance_seq(SEQ_A.c_str(), LENA, SEQ_B.c_str(), LENB) << endl; 
	auto t_middle = std::chrono::high_resolution_clock::now();
    cout << "Score en utilisant la version omp:\t\t" << nw_distance_omp(SEQ_A.c_str(), LENA, SEQ_B.c_str(), LENB) << endl;
	auto t_end = std::chrono::high_resolution_clock::now();

    cout << endl;
    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;

    duration<double, std::milli> ms_double_seq = t_middle - t_start;
    std::cout<<"temps SEQ : " << ms_double_seq.count() << "ms" << endl;

    duration<double, std::milli> ms_double_omp = t_end - t_middle;
    std::cout<<"temps OMP : " << ms_double_omp.count() << "ms" << endl;
    
    return 0;
}