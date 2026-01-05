#include "needleman_wunsch.hpp"
#include <iostream>

using namespace std;


int main () {
    int LENA = 100;
    string SEQ_A = "CCCCATGTCCTGAGGCTGCTCCCATCCCTTGGTTTCAGGGCTCCCTGTTTTATGGAGAAGCCTGATTTTATTCATTCTCATTCAATTCAGCCAACACTTG";
    
    int LENB = 100;
    string SEQ_B = "TCCATCCACCTTGGCCTCCCAAAGTGCTGGGATTTTAAGTGTGAGCCAATGCACCCGGCCCTATACCCAGGTTTTAATTTTTATTATTATTTTCTTTTGA";

    cout << "Partie Needleman-wunsch:" << endl;
    cout << "Sequence A: " << SEQ_A << endl;
    cout << "Sequence B: " << SEQ_B << endl;
    cout << "Score en utilisant la version sequentielle: " << nw_distance_seq(SEQ_A.c_str(), LENA, SEQ_B.c_str(), LENB) << endl; 
    return 0;
}