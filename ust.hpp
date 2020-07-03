//
//  ust.h
//  UST
//
//  Created by Amatur Rahman on 21/5/20.
//  Copyright © 2020 medvedevgroup. All rights reserved.
//
#ifndef ust_h
#define ust_h

#include<tuple>
#include<vector>
#include<string>
#include<map>
#include<unordered_map>


#include "param.hpp"
#include "stat.hpp"
#include "spss.hpp"

using namespace std;

class UST : public SPSS{
    
//different in UST and ESS
private:
    typedef struct {
        bool isWalkEnd = false;
        int pos_in_walk = -100;
        int finalWalkId = -1; // renders some walkId as invalid
    } new_node_info_t;
    UST::new_node_info_t* oldToNew;
    
    Param param;
    
    typedef tuple<int,int,int> threetuple; // uid, walkid, pos
    vector<threetuple> sorter;
    bool static sort_by_walkId (const threetuple &lhs, const threetuple &rhs){
        return get<1>(lhs) < get<1>(rhs);
    }
    bool static sort_by_pos (const threetuple &lhs, const threetuple &rhs){
        return get<2>(lhs) < get<2>(rhs);
    }
    
    USTStat stat;
    int countNewNode = 0; // number of paths in initial phase

    
    void collectInput(int argc, char** argv, string & graph_file_name, int & K, bool & FLG_ABUNDANCE);
    void readUnitigFile(const string& unitigFileName, vector<unitig_struct_t>& unitigs, vector<vector<edge_t> >& adjList);
    void indegreePopulateAndLowerBound();
    vector<threetuple> sorterMaker();
    void ustOutputToDisk(vector<threetuple>& sorter);

public:
    
    //void run(int argc, char** argv, bool runFlag);
    void run(string graph_file_name, int K, bool runFlag);
    ~UST();
};




void UST::run(string graph_file_name, int K, bool runFlag){
    
    cout<<"UST v2.0"<<endl;
    cout<<"--------------------"<<endl;
    
    this->FLG_ABUNDANCE = false;
    this->K= K;
    param.UNITIG_FILE = graph_file_name;
    //collectInput(argc, argv, graph_file_name, K, FLG_ABUNDANCE); //input: argc, argv
    
    cout<<"--------------------"<<endl;
    cout << "## Please wait while we read input file.... " << graph_file_name << ": k = "<<K<<endl;
    cout<<".........................."<<endl;
    double startTime = readTimer();
    
    this->readUnitigFile(param.UNITIG_FILE, unitigs, adjList); //input: graph_filename, output: last 2
    
    double TIME_READ_SEC = readTimer() - startTime;
    cout<<".........................."<<endl;
    cout<<"Done. TIME to read file "<<TIME_READ_SEC<<" sec."<<endl;
    cout<<"--------------------"<<endl;
    
    //initialization phase
    param.OUTPUT_FILENAME = getFileName(param.UNITIG_FILE)+".ust.fa";
    
    //system(("rm -rf "+getFileName(graph_file_name)+".stats.txt").c_str());
    globalStatFile.open((getFileName(param.UNITIG_FILE)+".stats.txt").c_str(), std::fstream::out);
    
    
    V_bcalm = adjList.size();
    nodeSign = new bool[V_bcalm];
    oldToNew = new new_node_info_t[V_bcalm];
    sortStruct = new struct node_sorter[V_bcalm];
    
    for (int i = 0; i < V_bcalm; i++) {
        
        if(ALGOMODE == TWOWAYEXT || ALGOMODE == BRACKETCOMP ){
            disSet.make_set(i);
        }
        
        oldToNew[i].finalWalkId = -1;
        sortStruct[i].sortkey = 0;
        sortStruct[i].node = i;
    }
    
    
    //stat collector
    stat.E_bcalm = 0;
    for (int i = 0; i < V_bcalm; i++) { //count total number of edges
        stat.E_bcalm += adjList[i].size();
    }
    
    stat.nKmers = 0;
    stat.C_bcalm = 0;
    for (unitig_struct_t unitig : unitigs) {    //count total number of kmers and characters
        stat.C_bcalm += unitig.ln;
        stat.nKmers +=  unitig.ln - K + 1;
    }
    
    stat.V_bcalm = V_bcalm;
    
    cout<<"--------------------"<<endl;
    cout<<"## Please wait while we gather info about lower bound.... "<<endl;
    cout<<".........................."<<endl;
    double time_a = readTimer();
    
    indegreePopulateAndLowerBound();
    
    cout<<".........................."<<endl;
    cout<<"Done. TIME to gather information about unitig graph: "<<readTimer() - time_a<<" sec."<<endl;
    cout<<"--------------------"<<endl;
    
    stat.walkstarting_node_count = ceil((stat.sharedparent_count + stat.sink_count + stat.source_count)/2.0) + stat.isolated_node_count;
    stat.charLowerbound = stat.C_bcalm-(K-1)*(V_bcalm - stat.walkstarting_node_count*1);
    stat.upperbound = (1-((stat.C_bcalm-(K-1)*(V_bcalm - stat.walkstarting_node_count*1.0))/stat.C_bcalm))*100.0;
    stat.statPrinter(globalStatFile, "K", K);
    stat.statPrinter(globalStatFile, "N_KMER", stat.nKmers);
    stat.statPrinter(globalStatFile, "V_BCALM", stat.V_bcalm);
    stat.statPrinter(globalStatFile, "E_BCALM", stat.E_bcalm);
    stat.statPrinter(globalStatFile, "C_BCALM", stat.C_bcalm);
    stat.statPrinter(globalStatFile, "C_LB", stat.charLowerbound);
    stat.statPrinter(globalStatFile, "V_LB", stat.walkstarting_node_count);
    stat.statPrinter(globalStatFile, "N_ISOLATED",  stat.isolated_node_count);
    stat.statPrinter(globalStatFile, "N_SINK",  stat.sink_count);
    stat.statPrinter(globalStatFile, "N_SOURCE", stat.source_count);
    stat.statPrinter(globalStatFile, "N_SPECIAL_NEW", stat.sharedparent_count);
    
    //OPTIONAL;DERIVABLE
    stat.statPrinter(globalStatFile, "BITSKMER_LB", (stat.charLowerbound*2.0)/stat.nKmers);
    stat.statPrinter(globalStatFile, "PERCENT_UB", stat.upperbound/100.0);
    stat.statPrinter(globalStatFile, "PERCENT_N_SPECIAL", stat.sharedparent_count*1.0/stat.V_bcalm);
    stat.statPrinter(globalStatFile, "PERCENT_N_ISOLATED", stat.isolated_node_count*1.0/stat.V_bcalm);
    stat.statPrinter(globalStatFile, "PERCENT_N_DEADEND", (stat.sink_count+stat.source_count)*1.0/stat.V_bcalm);
    
    // Iterating the map and printing ordered values
    //    for (auto i = inOutCombo.begin(); i != inOutCombo.end(); i++) {
    //        cout << "(" << i->first.first<< ", "<< i->first.second << ")" << " := " << i->second << '\n';
    //    }
    
    //if(ALGOMODE == PROFILE_ONLY){
    //printf("\n");
    //exit(1);
    //}
    
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
    //##################################################//
    
    //DFS
    cout<<"## Please wait while DFS step is going on... "<<endl;
    sorter = this->sorterMaker();
    
    time_a = readTimer();
    
    if(runFlag==1){
        return;
    }
    
    
    cout<<"## Writing strings to file.... "<<endl;
    this->ustOutputToDisk(sorter);
    
    
    cout<<"Done. TIME to output: "<<readTimer() - time_a<<" sec."<<endl;
    double TIME_TOTAL_SEC = readTimer() - startTime;
    
    
    if(param.VALIDATE){
        cout<<"## [4] Validating UST...\n";
        system((param.DSK_PATH +"dsk -file "+param.UNITIG_FILE+" -kmer-size "+to_string(K)+" -abundance-min 1  -verbose 0").c_str());
        system((param.DSK_PATH + "dsk -file absorbDecompressed.fa -kmer-size "+to_string(K)+" -abundance-min 1  -verbose 0").c_str());
        system((param.DSK_PATH+"dsk2ascii -file list_reads.unitigs.h5 -out output-bcalm.txt  -verbose 0").c_str());
        system((param.DSK_PATH + "dsk2ascii -file absorbDecompressed.h5 -out output-my.txt   -verbose 0").c_str());
        //cout<<"doing highly  accurate validation................"<<endl;
        system("sort -k 1 -n output-bcalm.txt -o a.txt; sort -k 1 -n output-my.txt -o b.txt");
        system("cmp a.txt b.txt && echo '### SUCCESS: Files Are Identical! ###' || echo '### WARNING: Files Are Different! ###'");
        system("rm -rf a.txt b.txt output-bcalm.txt output-my.txt list_reads.unitigs.h5 absorbDecompressed.h5");
    }
    
    
    
    float percent_saved_c = (1-(stat.C_ust*1.0/stat.C_bcalm))*1.0;
    float ustitchBitsPerKmer = stat.C_ust*2.0/stat.nKmers;
    
    
    stat.statPrinter(globalStatFile, "V_UST", stat.V_ust);
    stat.statPrinter(globalStatFile, "PERCENT_C_SAVED", percent_saved_c);
    stat.statPrinter(globalStatFile, "PERCENT_C_GAP_WITH_UB", stat.upperbound - percent_saved_c);
    stat.statPrinter(globalStatFile, "C_UST", stat.C_ust);
    stat.statPrinter(globalStatFile, "BITSKMER_UST", ustitchBitsPerKmer);
    
    
    globalStatFile << "TIME_READINPUT_SEC_" << mapmode[ALGOMODE].c_str() <<  "=" <<TIME_READ_SEC << endl;
    globalStatFile << "TIME_TOTAL_SEC_" << mapmode[ALGOMODE].c_str() <<  "=" <<TIME_TOTAL_SEC << endl;
    
    int maxulen =maximumUnitigLength();
    //stat.statPrinter(globalStatFile, "MAX_UNITIG_LEN MB", maxulen*1.0/1024.0/1024.0, true);
    //stat.statPrinter(globalStatFile, "E MB", stat.E_bcalm*8.0/1024.0/1024.0, true);
    //stat.statPrinter(globalStatFile, "V MB", stat.V_bcalm*8.0/1024.0/1024.0, true);
    
    globalStatFile.close();
    
    cout << "------------ UST completed successfully. Output is in file "<<param.OUTPUT_FILENAME << " ------------"<<endl;
    cout << "Total number of unique "<<K<<"-mers " <<  "= " << stat.nKmers << endl;
    cout << "Lower bound on any SPSS representation " <<  "= " << (stat.charLowerbound)*1.0/stat.nKmers << " neucleotide/k-mer"<< endl;
    
    cout << "Size of UST output" <<  "= " <<ustitchBitsPerKmer/2.0 << " neucleotide/k-mer"<< endl;
    //cout << "Gap with lower bound" << "= " << upperbound - percent_saved_c << "%" << endl;
    
    cout << "------------------------------------------------------"<<endl;
    
    
    
}


void UST::readUnitigFile(const string& unitigFileName, vector<unitig_struct_t>& unitigs, vector<vector<edge_t> >& adjList)
{
    //bcalm_file_type = 0
    
    ifstream unitigFile;
    if(!unitigFile.good()){
        fprintf(stderr, "Error: File named \"%s\" cannot be opened.\n", unitigFileName.c_str());
        exit(EXIT_FAILURE);
    }
    
    unitigFile.open(unitigFileName);
    string line;
    
    int nodeNum;
    char lnline[20];
    char kcline[20];
    char kmline[20];
    char edgesline[100000];
    bool doCont = false;
    
    int smallestK = 9999999;
    
    {
        getline(unitigFile, line);
        do {
            unitig_struct_t unitig_struct;
            
            if(FLG_ABUNDANCE){
                //>3 LN:i:24 ab:Z:10 10 10 10 10 7 7 7 7 7 3 3 3 3   L:-:0:+ L:-:1:+  L:+:0:-
                edgesline[0] = '\0';
                sscanf(line.c_str(), "%*c %d %s", &unitig_struct.serial, lnline);
                
                if(    line.find("ab:Z") == string::npos){
                    cout<<"Incorrect input format. Try using flag -a 0."<<endl;
                    exit(3);
                }
                
                sscanf(lnline, "%*5c %d", &unitig_struct.ln);
                
                int abpos = line.find("ab") + 5;
                int Lpos = line.find("L:");
                
                if(Lpos < 0){
                    Lpos = line.length() ;
                }
                // initialize string stream
                //cout<<line.substr(abpos, Lpos - abpos);
                stringstream ss(line.substr(abpos, Lpos - abpos));
                string abun;
                
                sscanf(line.substr(Lpos, line.length() - Lpos).c_str(), "%[^\n]s", edgesline);
                
                if(unitig_struct.ln < smallestK){
                    smallestK = unitig_struct.ln ;
                }
                if(unitig_struct.ln < K){
                    printf("Wrong k! Try again with correct k value. \n");
                    exit(2);
                }
                
            }else{
                edgesline[0] = '\0';
                sscanf(line.c_str(), "%*c %d %s  %s  %s %[^\n]s", &unitig_struct.serial, lnline, kcline, kmline, edgesline);
                
                if(    line.find("KC") == string::npos){
                    cout<<"Incorrect input format. Try using flag -a 1."<<endl;
                    exit(3);
                }
                
                //>0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:-
                sscanf(lnline, "%*5c %d", &unitig_struct.ln);
                
                
                if(unitig_struct.ln < smallestK){
                    smallestK = unitig_struct.ln ;
                }
                if(unitig_struct.ln < K){
                    printf("Wrong k! Try again with correct k value. \n");
                    exit(2);
                }
            }
            
            char c1, c2;
            stringstream ss(edgesline);
            
            vector<edge_t> edges;
            while (getline(ss, line, ' ')) {
                if (delSpaces(line).length() != 0) {
                    if(DBGFLAG==VERIFYINPUT){
                        cout<<line<<endl;
                    }
                    
                    sscanf(line.c_str(), "%*2c %c %*c %d  %*c  %c", &c1, &nodeNum, &c2); //L:-:0:-
                    edge_t newEdge;
                    
                    bool DELSELFLOOP=true;
                    if(DELSELFLOOP){
                        if((unitig_struct.serial)!= nodeNum){
                            newEdge.left = charToBool(c1);
                            newEdge.right = charToBool(c2);
                            newEdge.toNode = nodeNum;
                            edges.push_back(newEdge);
                        }
                    }else{
                        newEdge.left = charToBool(c1);
                        newEdge.right = charToBool(c2);
                        newEdge.toNode = nodeNum;
                        edges.push_back(newEdge);
                    }
                }
            }
            adjList.push_back(edges);
            
            doCont = false;
            while (getline(unitigFile, line)) {
                if (line.substr(0, 1).compare(">")) {
                    //unitig_struct.sequence = unitig_struct.sequence + line;
                    unitigs.push_back(unitig_struct);
                } else {
                    doCont = true;
                    break;
                }
            }
        } while (doCont);
        unitigFile.close();
    }
    
    
    if(smallestK > K ){
        cout<<"\n :::: :::: :::: :::: !!!!!!!!! WARNING !!!!!!!!!!! :::: :::: :::: ::::"<<endl;
        cout<<"The length of the smallest string we found was " << smallestK << ". Please make sure you are using the correct value of 'k' to ensure correctness of output."<<endl;
        cout << "------------------------------------------------------"<<endl;
    }
    //cout << "Complete reading input unitig file (bcalm2 file)." << endl;
}



void UST::collectInput(int argc, char** argv, string & graph_file_name, int & K, bool & FLG_ABUNDANCE){
    char * nvalue;
    char c;
    while( ( c = getopt (argc, argv, "i:k:a:") ) != -1 )
    {
        switch(c)
        {
            case 'a':
                if(optarg) {
                    if(strcmp(optarg, "0")==0 || strcmp(optarg, "1")==0){
                        FLG_ABUNDANCE = static_cast<bool>(std::atoi(optarg));
                    }else{
                        fprintf(stderr, "Error: use either -a 0 or -a 1 \n");
                        exit(EXIT_FAILURE);
                    }
                }
                break;
            case 'i':
                if(optarg) nvalue = optarg;
                break;
            case 'k':
                if(optarg) {
                    K = std::atoi(optarg) ;
                    if(K<=0){
                        fprintf(stderr, "Error: Specify a positive k value.\n");
                        exit(EXIT_FAILURE);
                    }
                }else{
                    fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name> [-a <0: for no counts, 1: to output separate count, default = 0>]\n",
                            argv[0]);
                    exit(EXIT_FAILURE);
                }
                break;
            default:
                fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name> [-a <0: for no counts, 1: to output separate count, default = 0>]\n",
                        argv[0]);
                exit(EXIT_FAILURE);
        }
    }
    
    if(K==0 || strcmp(nvalue, "")==0){
        fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name> [-a <0: for no counts, 1: to output separate count, default = 0>]\n",
                argv[0]);
        exit(EXIT_FAILURE);
    }
    graph_file_name = string(nvalue);
}


vector<UST::threetuple> UST::sorterMaker() {
    
    bool* saturated = new bool[V_bcalm];
    
    
    char* color = new char[V_bcalm];
    unitig_t * p_dfs = new unitig_t[V_bcalm];
    vector<list<unitig_t> > newToOld;
    
    double time_a = readTimer();
    
    for (unitig_t i = 0; i < V_bcalm; i++) {
        color[i] = 'w';
        p_dfs[i] = -1;
        saturated[i] = false;
    }
    
    cout<<"Basic V loop time: "<<readTimer() - time_a<<" sec"<<endl;
    time_a = readTimer();
    
    for (unitig_t j = 0; j < V_bcalm; j++) {
        unitig_t u;
        
        if(0==1){
            u = sortStruct[j].node;
        }else{
            u = j;
        }
        
        if (color[u] == 'w') {  //DFS_visit(u)
            
            unordered_map<int, vector<edge_t> > sinkSrcEdges; //int is the unitig id (old id)
            
            //            if(ALGOMODE == BRACKETCOMP){
            //                if(global_issinksource[u]==1){
            //                    vector<edge_t> adju = adjList.at(u);
            //                    vector<edge_t> myvector;
            //                    for (edge_t e : adju) {
            //                        myvector.push_back(e);
            //                    }
            //                    sinkSrcEdges[u] = myvector;
            //                    continue;
            //                }
            //            }
            
            stack<edge_t> s;
            edge_t uEdge;
            uEdge.toNode = u;
            s.push(uEdge);
            
            while (!s.empty()) {
                edge_t xEdge = s.top();
                
                int x = xEdge.toNode;
                s.pop();
                
                if (color[x] == 'w') {
                    color[x] = 'g';
                    s.push(xEdge);
                    vector<edge_t> adjx = adjList.at(x);
                    
                    // Now our branching code ::
                    // For a white x
                    // Consider 2 case:
                    // Case 1. p[x] = -1, it can happen in two way, x is the first one ever in this connected component, or no one wanted to take x
                    // either way, if p[x] = -1, i can be representative of a new node in new graph
                    // Case 2. p[x] != -1, so x won't be the representative/head of a newHome. x just gets added to its parent's newHome.
                    int u = unitigs.at(x).ln; //unitig length
                    assert(u >= K);
                    
                    if (p_dfs[x] == -1) {
                        list<unitig_t> xxx;
                        xxx.push_back(x);
                        newToOld.push_back(xxx);
                        oldToNew[x].finalWalkId = countNewNode++; // countNewNode starts at 0, then keeps increasing
                        oldToNew[x].pos_in_walk = 1;
                    } else {
                        newToOld[oldToNew[p_dfs[x]].finalWalkId].push_back(x);
                        oldToNew[x].finalWalkId = oldToNew[p_dfs[x]].finalWalkId;
                        
                        { //ALGOMODE==TWOWAYEXT || ALGOMODE==BRACKETCOMP
                            disSet.Union(x, p_dfs[x]);
                        }
                        
                        oldToNew[x].pos_in_walk = oldToNew[p_dfs[x]].pos_in_walk + 1;
                    }
                    
                    // x->y is the edge, x is the parent we are extending
                    for (edge_t yEdge : adjx) { //edge_t yEdge = adjx.at(i);
                        int y = yEdge.toNode;
                        
                        //                        if(ALGOMODE == BRACKETCOMP){
                        //                            if(global_issinksource[y] == true){
                        //                                continue;
                        //                            }
                        //                        }
                        
                        if (color[y] == 'w') { //Normal DFS
                            s.push(yEdge);
                        }
                        
                        if (y == x) { // self-loop
                            cout<<"FAIL: should not have self-loop."<<endl;
                            assert(false);
                            edge_both_t e;
                            e.edge = yEdge;
                            e.fromNode = x;
                        } else if (saturated[x]) {
                            // Since x is saturated, we only add resolveLater edges
                            // no need to check for consistency
                            if (y != p_dfs[x]) {
                                edge_both_t e;
                                e.edge = yEdge;
                                e.fromNode = x;
                            }
                        } else {
                            // If x has space to take a child, meaning x is not saturated
                            // hunting for potential child
                            
                            if (color[y] == 'w' && p_dfs[y] == -1) {
                                // y has white color & got no parent => means it's homeless, so let's see if we can take it as a child of x
                                //But just see if it is eligible to be a child, i.e. is it consistent (sign check)?
                                
                                //2 case, Does x's child have grandparent?
                                // If No:
                                if (p_dfs[x] == -1) {
                                    // case 1: child has no grandparent
                                    // so extend path without checking any sign
                                    nodeSign[x] = yEdge.left;
                                    nodeSign[y] = yEdge.right;
                                    p_dfs[y] = x;
                                    saturated[x] = true; //found a child
                                    
                                } else if (nodeSign[x] == yEdge.left) {
                                    // case 2: child (=y) has grandparent, i.e. x's parent exists
                                    nodeSign[y] = yEdge.right;
                                    p_dfs[y] = x;
                                    saturated[x] = true; //found a child
                                    
                                } else {
                                    // do we reach this case?
                                    edge_both_t e;
                                    e.edge = yEdge;
                                    e.fromNode = x;
                                }
                                
                            } else {
                                //merger
                                {
                                    // y is not white
                                    bool consistentEdge = (nodeSign[y] == yEdge.right && (p_dfs[x]==-1 || (p_dfs[x]!=-1&& nodeSign[x] == yEdge.left)) );
                                    if(p_dfs[y]==-1 && consistentEdge && oldToNew[x].finalWalkId != oldToNew[y].finalWalkId){
                                        
                                        //cout<<"x: "<<x<<":" <<disSet.find_set(x)<<" ";
                                        //cout<<"y: "<<y<<":" <<disSet.find_set(y) <<endl;
                                        
                                        //not in same group already, prevent cycle
                                        if(disSet.find_set(x)!=disSet.find_set(y)){
                                            nodeSign[x] = yEdge.left;
                                            nodeSign[y] = yEdge.right;
                                            p_dfs[y] = x;
                                            saturated[x] = true; //found a child
                                            // oldToNew[y].serial
                                            
                                            disSet.Union(x, y);
                                            gmerge.connectGroups(oldToNew[x].finalWalkId,oldToNew[y].finalWalkId );
                                            
                                        }
                                    }
                                }
                                
                                if (y != p_dfs[x]) {
                                    edge_both_t e;
                                    e.edge = yEdge;
                                    e.fromNode = x;
                                }
                            }
                        }
                    }
                } else if (color[x] == 'g') {
                    color[x] = 'b';
                }
            }
        }
    }
    delete [] color;
    delete [] p_dfs;
    delete [] saturated;
    cout<<"DFS time: "<<readTimer() - time_a<<" sec"<<endl;
    
    
    /***MERGE START***/
    bool* merged = new bool[countNewNode];  // now we will do union-find with path compresison for both way merge
    for (int i = 0; i<countNewNode; i++) {
        merged[i] = false;
    }
    
    { //both way merging
        for ( const auto& p: gmerge.fwdWalkId)
        {
            if(gmerge.fwdVisited[p.first] == false){
                int fromnode =p.first;
                int tonode = p.second;
                deque<int> lst;
                
                lst.push_back(fromnode);
                lst.push_back(tonode);
                
                gmerge.fwdVisited[fromnode] = true;
                gmerge.bwdVisited[tonode] = true;
                
                if(gmerge.fwdVisited.count(tonode)>0){
                    while(gmerge.fwdVisited[tonode] == false){
                        gmerge.fwdVisited[tonode] = true;
                        tonode = gmerge.fwdWalkId[tonode];
                        gmerge.bwdVisited[tonode] = true;
                        
                        lst.push_back(tonode);
                        if(gmerge.fwdVisited.count(tonode)==0)
                            break;
                    }
                }
                if(gmerge.bwdVisited.count(fromnode)>0){
                    while(gmerge.bwdVisited[fromnode] == false){
                        gmerge.bwdVisited[fromnode] = true;
                        fromnode = gmerge.bwdWalkId[fromnode];
                        gmerge.fwdVisited[fromnode] = true;
                        
                        lst.push_front(fromnode);
                        if(gmerge.bwdVisited.count(fromnode)==0)
                            break;
                    }
                }
                
                assert(!lst.empty());
                int commonWalkId = lst.at(0);
                int posOffset = 1;
                
                for(auto i: lst){
                    // i is new walk id before merging
                    merged[i] = true;
                    // travesing the walk list of walk ID i
                    for(int uid: newToOld[i]){
                        oldToNew[uid].finalWalkId = commonWalkId;
                        oldToNew[uid].pos_in_walk = posOffset++;
                    }
                }
                oldToNew[newToOld[lst.back()].back()].isWalkEnd = true;
                stat.V_ust++;
            }
        }
        
        for (int newNodeNum = 0; newNodeNum<countNewNode; newNodeNum++){
            if(merged[newNodeNum] == false){
                oldToNew[newToOld[newNodeNum].back()].isWalkEnd = true;
                stat.V_ust++;
            }
        }
        delete [] merged;
        vector<list<unitig_t> >().swap(newToOld);
    }
    
    //sorter of all walks and printing them
    vector<threetuple> sorter;
    for(int uid = 0 ; uid< V_bcalm; uid++){
        new_node_info_t nd = oldToNew[uid];
        //sorter.push_back(make_tuple(uid, nd.finalWalkId, nd.pos_in_walk, nd.isTip));
        sorter.push_back(make_tuple(uid, nd.finalWalkId, nd.pos_in_walk));
    }
    
    stable_sort(sorter.begin(),sorter.end(),sort_by_pos);
    stable_sort(sorter.begin(),sorter.end(),sort_by_walkId);
    return sorter;
}

void UST::ustOutputToDisk(vector<threetuple>& sorter){
    string uidSeqFilename = "uidSeq.usttemp"; //"uidSeq"+ mapmode[ALGOMODE] +".txt"
    ofstream uidSequence;
    uidSequence.open(uidSeqFilename);
    
    int finalUnitigSerial = 0;
    for(threetuple n : sorter){
        int uid = get<0>(n);
        int bcalmid = unitigs.at(uid).serial;
        uidSequence << finalUnitigSerial <<" "<< bcalmid << endl;
        finalUnitigSerial++;
    }
    uidSequence.close();
    
    //keep the sequences only
    system(("awk '!(NR%2)' "+param.UNITIG_FILE+" > seq.usttemp").c_str());
    system("sort -n -k 2 -o uidSeq.usttemp uidSeq.usttemp");
    if(FLG_ABUNDANCE){
        system(("awk '(NR%2)' "+param.UNITIG_FILE+" | cut -f 5 -d ':' | cut -f 1 -d 'L' > count.usttemp").c_str()); // get a separate count file
        system("paste -d' ' uidSeq.usttemp seq.usttemp count.usttemp > merged.usttemp ");
        system("sort -n -k 1 -o merged.usttemp merged.usttemp");
        system(("cat  merged.usttemp  | awk '{for (i=4;i<=NF;i+=1) print $i}' > "+getFileName(param.UNITIG_FILE)+".ust.counts").c_str());
    }else{
        system("paste -d' ' uidSeq.usttemp seq.usttemp > merged.usttemp ");
        system("sort -n -k 1 -o merged.usttemp merged.usttemp");
    }
    system("cat  merged.usttemp  | cut -d' ' -f3 >  seq.usttemp");
    
    
    ifstream sequenceStringFile ("seq.usttemp");
    ofstream ustOutputFile (param.OUTPUT_FILENAME);
    
    
    //both string and abundance sort
    //keep string only and output
    //open the string file
    
    int lastWalk = -1;
    string walkString = "";
    string unitigString = "";
    for(threetuple n : sorter){
        int uid = get<0>(n);
        int finalWalkId = get<1>(n);
        
        //for each line in file
        string sequenceFromFile = "";//getline
        getline (sequenceStringFile,sequenceFromFile);
        if(nodeSign[uid] == false){
            unitigString =  reverseComplement(sequenceFromFile);
        }else{
            unitigString =  sequenceFromFile;
        }
        
        if(finalWalkId!=lastWalk){
            if(lastWalk != -1){
                //print previous walk
                if(walkString.length()>=K){
                    ustOutputFile<<">"<<endl;
                    stat.C_ust+=walkString.length();
                    ustOutputFile<< walkString<<endl;
                }
            }
            
            //start a new walk
            walkString = "";
            lastWalk = finalWalkId;
        }
        walkString = plus_strings(walkString, unitigString, K);
        //ustOutputFile<<">"<<uid<<" " <<finalWalkId<<" "<<pos_in_walk<<endl;
    }
    if(walkString.length()>=K){
        ustOutputFile<<">"<<endl;
        stat.C_ust+=walkString.length();
        ustOutputFile<< walkString<<endl;
    }else{
    }
    
    sequenceStringFile.close();
    //smallKFile.close();
    system("rm -rf *.usttemp");
    system("rm -rf uidSeq.txt");
    
    
    // UST (TWOWAYEXT) DONE!
    //delete []  global_issinksource;
}


void UST::indegreePopulateAndLowerBound(){
    map<pair <int, int>, int> inOutCombo;
    
    int* global_indegree;
    int* global_outdegree;
    int* global_plusindegree;
    int* global_plusoutdegree;
    bool* global_issinksource;
    
    global_indegree = new int[V_bcalm];
    global_outdegree = new int[V_bcalm];
    global_plusindegree = new int[V_bcalm];
    global_plusoutdegree = new int[V_bcalm];
    global_issinksource = new bool[V_bcalm];
    bool* countedForLowerBound = new bool[V_bcalm];
    
    
    for (unitig_t i = 0; i < V_bcalm; i++) {
        global_indegree[i] = 0;
        global_outdegree[i] = 0;
        global_plusindegree[i] = 0;
        global_plusoutdegree[i] = 0;
        global_issinksource[i] = false;
        countedForLowerBound[i] = false;
    }
    
    int xc = 0;
    for(vector<edge_t> elist: adjList){
        for(edge_t e: elist){
            global_indegree[e.toNode] += 1;
            sortStruct[e.toNode].sortkey = sortStruct[e.toNode].sortkey + 1;
            if(e.right == true){
                global_plusindegree[e.toNode] += 1;
            }
            if(e.left == true){
                global_plusoutdegree[xc] += 1;
            }
        }
        global_outdegree[xc] = elist.size();
        xc++;
    }
    
    for(int i = 0; i<5; i++){
        for(int j = 0; j<5; j++){
            inOutCombo[make_pair(i,j)] = 0;
        }
    }
    
    stat.sink_count  = 0;
    
    for(int i = 0; i<V_bcalm; i++){
        pair<int, int> a;
        a = make_pair(global_plusindegree[i], global_plusoutdegree[i]);
        inOutCombo[a] = (inOutCombo.count(a)  ? inOutCombo[a] + 1 : 1  );
        if(global_plusoutdegree[i] == 0 && global_plusindegree[i] != 0){
            stat.sink_count++;
            global_issinksource[i] = 1;
            countedForLowerBound[i] = true;
        }
        if(global_plusindegree[i] == 0 && global_plusoutdegree[i] != 0){
            stat.source_count++;
            global_issinksource[i] = 1;
            countedForLowerBound[i] = true;
        }
        if(global_indegree[i] == 0){
            global_issinksource[i] = 1;
            stat.isolated_node_count++;
        }
        if(global_indegree[i] == 1){
            stat.onecount++;
        }
    }
    
    
    for (auto i = inOutCombo.begin(); i != inOutCombo.end(); i++) {
        globalStatFile << "PERCENT_DEGREE_"<<(i->first).first << "_" << (i->first).second <<  "=" << (i->second)*100.0/V_bcalm <<"%" << endl;
    }
    
    xc = 0; // current vertex while traversing the adjacency list
    for(vector<edge_t> elist: adjList){
        int neighborCount = 0;
        int spNeighborCount[2];
        spNeighborCount[0]=0;
        spNeighborCount[1]=0;
        stack<int> countedNodes;
        set<pair<int, bool> > countedSides;
        //if(true){
        if(FLG_NEWUB == true){
            //ENDPOINT SIDE UPPER BOUND - improved
            for(edge_t e_xy: elist){    //central node: all neighbors of x
                int y = e_xy.toNode;
                vector<edge_t> adjY = adjList[y];
                bool eligibleSp = true;
                
                //pair<int, bool> pairr;
                for(edge_t e_from_y : adjY){    // check if this neighbor is speacial
                    //pairr =make_pair(e_from_y.toNode, sRight(e_xy) );
                    if(e_from_y.toNode!=xc){
                        
                        if(sRight(e_xy) == sLeft(e_from_y)){
                            eligibleSp = false;
                            break;
                        }
                    }
                }
                
                if(eligibleSp){
                    spNeighborCount[sLeft(e_xy)]++;
                }
            }
            if(spNeighborCount[0]>1){
                stat.sharedparent_count += spNeighborCount[0] - 1 ;
            }
            if(spNeighborCount[1]>1){
                stat.sharedparent_count += spNeighborCount[1] - 1 ;
            }
        }
        if(FLG_NEWUB == false){
            //ENDPOINT SIDE UPPER BOUND
            for(edge_t e_xy: elist){
                int y = e_xy.toNode;
                vector<edge_t> adjY = adjList[y];
                bool eligible = true;
                pair<int, bool> pairr;
                for(edge_t e_from_y : adjY){
                    pairr =make_pair(e_from_y.toNode, sRight(e_xy) );
                    if(e_from_y.toNode!=xc){
                        
                        if(sRight(e_xy) == sLeft(e_from_y)){
                            eligible = false;
                            break;
                        }
                        
                    }
                    
                }
                
                if(eligible){
                    neighborCount++;
                }
            }
            
            if(global_issinksource[xc] == 1){
                if(neighborCount>1){
                    stat.sharedparent_count += neighborCount - 1 ;
                }
            }else{
                if(neighborCount>2){
                    stat.sharedparent_count += neighborCount - 2 ;
                }
            }
        }
        //sharedparent_count_wrong =sharedparent_count;
        
        //if(true){
        if(1==0){
            // OLDER UPPER BOUND CALC
            
            int neighborCount = 0;
            for(edge_t e_xy: elist){
                int y = e_xy.toNode;
                
                if(!countedForLowerBound[y]){
                    vector<edge_t> adjY = adjList[y];
                    bool eligible = true;
                    for(edge_t e_from_y : adjY){
                        if(e_from_y.toNode!=xc){
                            if(sRight(e_xy) == sLeft(e_from_y) ){
                                eligible = false;
                                break;
                            }
                        }
                    }
                    if(eligible){
                        countedForLowerBound[y] = true;
                        //global_priority[y] = 4;
                        neighborCount++;
                        countedNodes.push(y);
                    }
                }
            }
            
            if(global_issinksource[xc] == 1){
                if(neighborCount>1){
                    stat.sharedparent_count += neighborCount - 1 ;
                }else{
                    while(!countedNodes.empty()){
                        countedForLowerBound[countedNodes.top()] = false;
                        countedNodes.pop();
                    }
                }
            }else{
                if(neighborCount>2){
                    stat.sharedparent_count += neighborCount - 2 ;
                }else{
                    while(!countedNodes.empty()){
                        countedForLowerBound[countedNodes.top()] = false;
                        countedNodes.pop();
                    }
                }
            }
        }
        
        xc++;
    }
    
    delete [] global_indegree;
    delete [] global_outdegree;
    delete [] global_plusindegree;
    delete [] global_plusoutdegree;
    delete [] countedForLowerBound;
}


UST::~UST(){
    delete [] nodeSign;
    delete [] oldToNew;
    delete [] sortStruct;
}



int equal_files(){
    int result = system("diff /Users/Sherlock/Documents/bcl/bcl/a.txt /Users/Sherlock/Documents/bcl/bcl/output-my.txt >nul 2>nul");
    //They are different if result != 0
    return result;
}


#endif /* ust_h */
