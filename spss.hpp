//
//  spss.hpp
//  Created by Amatur Rahman on 21/5/20.
//  Copyright © 2020 medvedevgroup. All rights reserved.
//

#ifndef spss_h
#define spss_h
using namespace std;
#include "misc.hpp"

typedef int32_t unitig_t;

enum DEBUGFLAG_T { NONE = 0,  UKDEBUG = 6, VERIFYINPUT = 1, INDEGREEPRINT = 2, DFSDEBUGG = 3, PARTICULAR = 4, NODENUMBER_DBG = 5, OLDNEWMAP = 9, PRINTER = 10, SINKSOURCE = 12};

enum ALGOMODE_T { TWOWAYEXT = 10, PROFILE_ONLY = 11, BRACKETCOMP = 15};

DEBUGFLAG_T DBGFLAG = NONE; //NODENUMBER_DBG
ALGOMODE_T ALGOMODE = TWOWAYEXT;

string tToMode[] = {"ess", "esstip", "ust"};

string mapmode[] = {"basic", "indegree_dfs", "indegree_dfs_initial_sort_only", "outdegree_dfs", "outdegree_dfs_initial_sort_only", "inverted_indegree_dfs", "plus_indegree_dfs", "random_dfs", "node_assign", "source_first", "twoway", "profile_only", "endpoint_priority", "graph_print", "tight_ub", "tip"
};

string modefilename[] = {"Fwd", "indegree_dfs", "indegree_dfs_initial_sort_only", "outdegree_dfs", "outdegree_dfs_initial_sort_only", "inverted_indegree_dfs", "plus_indegree_dfs", "random_dfs", "node_assign", "source_first", "", "profile_only", "endpoint_priority", "graph_print", "tight_ub", "Tip"
};

class SPSS{
public:
    typedef struct {
        bool left;  //1 means +, 0 means -
        bool right; //1 means +, 0 means -
        int toNode;
    } edge_t;
    vector<vector<edge_t> > adjList;
    
    typedef struct {
        edge_t edge;
        int fromNode;
    } edge_both_t;
    
    typedef struct {
        int serial;
        int ln;
    } unitig_struct_t;
    vector<unitig_struct_t> unitigs;
    
    
    struct node_sorter {
        int node;
        int sortkey;
        //bool operator() (struct node_sorter  i, struct node_sorter  j) { return (i.sortkey<j.sortkey);}
    };
    bool sort_by_key (struct node_sorter i, struct node_sorter j) { return (i.sortkey<j.sortkey); }
    bool sort_by_key_inverted (struct node_sorter i, struct node_sorter j) { return (i.sortkey>j.sortkey); }
    struct node_sorter * sortStruct;
    

    class GroupMerger {
        public:
            map<int, bool> fwdVisited;
            map<int, bool> bwdVisited;
            map<int, int> bwdWalkId;
            map<int, int> fwdWalkId;
            GroupMerger() {
            }
            void connectGroups(int from, int to){
                fwdVisited[from] = false;
                bwdVisited[to] = false;
                fwdWalkId[from] = to;
                bwdWalkId[to] = from;
            }
            ~GroupMerger() {
            }
        };
        GroupMerger gmerge;
    
        class DisjointSet {
            unordered_map<int, int> parent;
            
        public:
            DisjointSet() {
            }
            void make_set(int id) {
                this->parent[id] = -1;
            }
            
            void Union(int xId, int yId) {
                int xset = find_set(xId);
                int yset = find_set(yId);
                
                if(xset != yset)
                {
                    parent[xset] = yset;
                }
            }
            
            int find_set(int id) {
                if (parent[id] == -1)
                    return id;
                return find_set(parent[id]);
            }
            ~DisjointSet(){
            }
        };
        DisjointSet disSet;
    
        bool charToBool(char c) {
            if (c == '+') {
                return true;
            } else {
                if (c != '-') cout << "Erroneus character in BCALM2 output." << endl;
                return false;
            }
        }
        inline char boolToCharSign(bool sign) {
            return (sign == true) ? '+' : '-';
        }
        inline bool sRight(edge_t plusminusedge){
            return !(plusminusedge.right == true);
        }
        inline bool sLeft(edge_t plusminusedge){
            return (plusminusedge.left == true);
        }
        
        inline int maximumUnitigLength(){
            int m = 0;
            for(unitig_struct_t u: unitigs){
                if(u.ln > m){
                    m = u.ln;
                }
            }
            return m;
        }
    
//        void printBCALMGraph(vector<vector<edge_t> >& adjList) {
//            for (int i = 0; i < adjList.size(); i++) {
//                cout << i << "# ";
//                for (edge_t edge : adjList.at(i)) {
//                    cout << boolToCharSign(edge.left) << ":" << edge.toNode << ":" << boolToCharSign(edge.right) << ", ";
//                }
//                cout << endl;
//            }
//        }
    
    
    /**/
    //FILENAMES, IO
    //filenames
    ofstream globalStatFile;
        
    //POINTERS
    bool* nodeSign;
    
    
    unitig_t V_bcalm = 0;
    int K = 0;
    bool FLG_NEWUB = true;
    bool FLG_ABUNDANCE = false;
};


#endif /* spss_h */
