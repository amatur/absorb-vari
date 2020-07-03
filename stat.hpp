//
//  stat.hpp
//  UST
//
//  Created by Amatur Rahman on 23/5/20.
//  Copyright © 2020 medvedevgroup. All rights reserved.
//

#ifndef stat_h
#define stat_h

#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
using namespace std;

class Stat{
public:
    double time_input  = -1;
    double time_output  = -1;
    
    int nUnitigs = -1;
    uint_fast64_t nKmers = -1;
    
    uint_fast64_t C_bcalm = 0;
    uint_fast64_t E_bcalm = 0;
    uint_fast64_t V_bcalm = 0;
    
    int charLowerbound = 0;
    
    int isolated_node_count = 0;
    int sink_count = 0;
    int source_count = 0;
    int sharedparent_count = 0;
    int sharparentCntRefined = 0;
    int onecount = 0;
    
    int walkstarting_node_count = 0;
    float upperbound = 0;
    
    template <class T>
    void statPrinter(ofstream& FOUT, string label, T a, bool consoleOutput = false){
        FOUT<<label<<"="<<a<<endl;
        if(consoleOutput){
            cout<<label<<"="<<a<<endl;
        }
    }
};

class USTStat: public Stat{
    public:
    int V_ust = 0;
    int C_ust = 0;
};

class ESSTipStat: public Stat{
    public:
    int V_esstip = 0;
    int C_esstip = 0;

    int V_nontip = 0;
    int C_nondna_esstip = 0;
};


class ESSStat: public Stat{
    public:
    int V_ess = 0;
    int C_ess = 0;
    
    int V_ust = 0;
    
    int V_non_absorbed = 0;
    int C_nondna_ess = 0;
    
    
    int absorbGraphNumCC_endAbosrb = 0;
};

vector<string> split (const string &s, char delim) {
    vector<string> result;
    stringstream ss (s);
    string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }

    return result;
}

//void substat(string globalStatFileName){
//    ifstream gl(globalStatFileName);
//    string fileline;
//    while (getline (gl, fileline)) {
//        vector<string> v = split (fileline, '=');
//        v[0]
//        for (auto i : v) cout << i << endl;
//        result.push_back (item);
//    }
//
//    getline(unitigFile, line);
//
//
//}



#endif /* stat_h */
