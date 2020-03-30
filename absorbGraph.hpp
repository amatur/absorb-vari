//
//  absorbGraph.hpp

#ifndef absorbGraph_hpp
#define absorbGraph_hpp

#include "global.hpp"
#include <algorithm>



void hasherSmallerK(int node){
    for(int smallOverlap = K; smallOverlap >= SMALLK; smallOverlap--) {
        string seq = unitigs.at(node).sequence;
        //int SMALLK = 20;
        string pref = seq.substr(0, smallOverlap-1);
        string suf = seq.substr(seq.length()-(smallOverlap-1), smallOverlap-1);
        
        plusMapPre[pref].insert(node);
        //plusMapSuf[suf].insert(node);
        
        string rseq = reverseComplement(unitigs.at(node).sequence);
        string rpref = rseq.substr(0, smallOverlap-1);
        string rsuf = rseq.substr(rseq.length()-(smallOverlap-1), smallOverlap-1);
        
        minusMapPre[rpref].insert(node);
        //minusMapSuf[rsuf].insert(node);
        
    }
}
void addSmallerKEdge(int node){
    //int SMALLK= 20;
    
    for(int smallOverlap = K; smallOverlap >= SMALLK; smallOverlap--) {
        string seq = unitigs.at(node).sequence;
        string pref = seq.substr(0, smallOverlap-1);
        string suf = seq.substr(seq.length()-(smallOverlap-1), smallOverlap-1);
        
        string rseq = reverseComplement(unitigs.at(node).sequence);
        string rpref = rseq.substr(0, smallOverlap-1);
        string rsuf = rseq.substr(rseq.length()-(smallOverlap-1), smallOverlap-1);
        
        
        if (plusMapPre.find(suf) != plusMapPre.end()){
            set<int> x = plusMapPre[suf];
            for( int i : x){
                if(i!=node){
                    edge_t e, e2;
                    e.toNode = node;
                    e.left = false;
                    e.right = false;
                    e.isSmallK = 'y';
                    adjList[i].push_back(e);
                    
                    e2.toNode = i;
                    e2.left = true;
                    e2.right = true;
                     e2.isSmallK = 'y';
                    adjList[node].push_back(e2);
                }
            }
            
        }
        
        
        if (minusMapPre.find(suf) != minusMapPre.end()){
            set<int> x = minusMapPre[suf];
            for( int i : x){
                if(i!=node){
                    edge_t e, e2;
                    e.toNode = node;
                    e.left = true;
                    e.right = false;
                     e.isSmallK = 'y';
                    adjList[i].push_back(e);
                    
                    e2.toNode = i;
                    e2.left = true;
                    e2.right = false;
                     e2.isSmallK = 'y';
                    adjList[node].push_back(e2);
                }
            }
            
        }
        
        if (minusMapPre.find(rsuf) != minusMapPre.end()){
            set<int> x = minusMapPre[rsuf];
            for( int i : x){
                if(i!=node){
                    edge_t e, e2;
                    e.toNode = node;
                    e.left = true  ;
                    e.right = true;
                     e.isSmallK = 'y';
                    adjList[i].push_back(e);
                    
                    e2.toNode = i;
                    e2.left = false;
                    e2.right = false;
                     e2.isSmallK = 'y';
                    adjList[node].push_back(e2);
                }
            }
            
        }
        
        if (plusMapPre.find(rsuf) != plusMapPre.end()){
            set<int> x = plusMapPre[rsuf];
            for( int i : x){
                if(i!=node){
                    edge_t e, e2;
                    e.toNode = node;
                    e.left = false;
                    e.right = true;
                     e.isSmallK = 'y';
                    adjList[i].push_back(e);
                    
                    e2.toNode = i;
                    e2.left = false;
                    e2.right = true;
                     e2.isSmallK = 'y';
                    adjList[node].push_back(e2);
                }
            }
            
        }
    }
  
}


ofstream depthfile ("depthfile.txt");

void ccDFS(int start, bool* ccVisited, int &count, vector<int>* ccAdjList)
{
    int x=  countNewNode;
  ccVisited[start] = true;
  count++;

      //  vector<int> a=ccAdjList.at(start);
    
    
    vector<int>::iterator i;

       for(i = ccAdjList[start].begin(); i != ccAdjList[start].end(); ++i){
           if(!ccVisited[*i]){
               
           }
               //nodes.push_back(elem);
               //ccDFS(*i, ccVisited, count, ccAdjList);
       }
//    } catch (const std::exception& e) { // reference to the base of a polymorphic object
//         std::cout << e.what(); // information from length_error printed
//    }
    


//    for (int elem :a )
//    {
//        if(ccVisited[elem] == false){
//            //nodes.push_back(elem);
//            ccDFS(elem, ccVisited, count, ccAdjList);
//        }
//
//    }
    
    
}
//
void connectedCC_allAbsorb(){
    return;
    absorbGraphNumCC_allAbosrb = 0;
    vector<vector<int> > ccAdjList(countNewNode);
    bool* ccVisited =  new bool[countNewNode];
    for(int i = 0; i<countNewNode; i++){
        ccVisited[i] = (false);
    }
    
    
    for(int uid = 0; uid <adjList.size(); uid++){
        vector<edge_t> adju = adjList.at(uid);
        for (edge_t e : adju) {
            int absorberWalk = oldToNew[e.toNode].finalWalkId;
            int absorbedWalk = oldToNew[uid].finalWalkId;
            
            if(absorberWalk != absorbedWalk){
                ccAdjList[absorberWalk].push_back(absorbedWalk);
                ccAdjList[absorbedWalk].push_back(absorberWalk);
            }
        }
    }
    
    if(1==1){
        stack<int> stack;
        for(int ii=0;ii<countNewNode;ii++)
        {
            if(ccVisited[ii] == false &&  !obsoleteWalkId[ii])
            {
                //cout<<endl;
                stack.push(ii);
                //vector<int>nodes;
                //nodes.push_back(i);
                absorbGraphNumCC_allAbosrb++;
                
                while (!stack.empty())
                {
                    // Pop a vertex from stack and print it
                    int s = stack.top();
                    stack.pop();
                    
                    // Stack may contain same vertex twice. So
                    // we need to print the popped item only
                    // if it is not visited.
                    if (!ccVisited[s])
                    {
                        //cout << s << " ";
                        ccVisited[s] = true;
                    }
                    
                    // Get all adjacent vertices of the popped vertex s
                    // If a adjacent has not been visited, then push it
                    // to the stack.
                    for (auto i = ccAdjList[s].begin(); i != ccAdjList[s].end(); ++i)
                        if (!ccVisited[*i])
                            stack.push(*i);
                }
            }
            
        }
        cout<<"number of connected components (all): "<<absorbGraphNumCC_allAbosrb<<endl;
        delete [] ccVisited;
        
    }
}

//make the absorb graph with cycle
void sorterIndexAndAbsorbGraphMaker(vector<MyTypes::fourtuple>& sorter, map<int, int >& sorterIndexMap, map<int, vector<edge_t> > & absorbGraph, char*& absorbedCategory){
    int prevWalkId = -1;
    int lastWalkStartingIndex = -1;
    
    //list<int> *adj
    vector<int> *ccAdjList = new vector<int>[countNewNode];
    //vector<vector<int> > ccAdjList(countNewNode);
    //block for connected compo
   bool* ccVisited = new bool[countNewNode];
    
    bool* isAnAbsorbedWalk = new bool[countNewNode];
    
    for(int i = 0; i<countNewNode; i++){
        ccVisited[i]=false;
        isAnAbsorbedWalk[i]=false;
                  //set<int> myset;
                  //ccAdjList.push_back(myset);
              }
    
    
    for(int tup_i = 0; tup_i < sorter.size(); tup_i++){
        
        MyTypes::fourtuple tup = sorter[tup_i];
        int uid = get<0>(tup);
        new_node_info_t nd = oldToNew[uid];
        int finalWalkId = get<1>(tup);
        
        
        //for end abosrbing
        if(prevWalkId !=finalWalkId ){  //walk starting: 2 cases-> a) this is not walk ender b) ender
            lastWalkStartingIndex = tup_i;
            sorterIndexMap[finalWalkId] = lastWalkStartingIndex;
            
            vector<edge_t> adju = adjList.at(uid);
            for (edge_t e : adju) {
                int absorberWalk = oldToNew[e.toNode].finalWalkId;
                int absorbedWalk = nd.finalWalkId;

                
                //@change_for_endabsorb
                if(absorberWalk != absorbedWalk && oldToNew[e.toNode].pos_in_walk!=1 && isAnAbsorbedWalk[oldToNew[uid].finalWalkId] == false){
                //if(absorberWalk != absorbedWalk){
                    //the edge "enew" is from absorber uid to absorbed uid
                    assert(absorberWalk<countNewNode);
                    assert(absorbedWalk<countNewNode);

                    ccAdjList[absorberWalk].push_back(absorbedWalk);
                    ccAdjList[absorbedWalk].push_back(absorberWalk);
                    //
                    
                    edge_t enew;
                    if(e.left==e.right){
                        enew.left = (e.left?false:true);
                        enew.right = (e.right?false:true);
                    }else{
                        enew.left = e.left;
                        enew.right = e.right;
                    }
                    
                    enew.toNode = uid;
                    absorbGraph[e.toNode].push_back(enew);
                    isAnAbsorbedWalk[oldToNew[uid].finalWalkId]=true;
                    //cout<<absorberWalk<<"->absorbs->"<<absorbedWalk<<endl;
                    //TODO: need to fix for decoding
                    //cout<<"hey"<<endl;
                    //break;
                }
            }
        }else{
            if(ALLABSORBING && isAnAbsorbedWalk[oldToNew[uid].finalWalkId]==false){
                vector<edge_t> adju = adjList.at(uid);
                for (edge_t e : adju) {
                    int absorberWalk = oldToNew[e.toNode].finalWalkId;
                    int absorbedWalk = nd.finalWalkId;

                    
                    //@change_for_endabsorb
                    if(absorberWalk != absorbedWalk && oldToNew[e.toNode].pos_in_walk!=1){
                    //if(absorberWalk != absorbedWalk){
                        //the edge "enew" is from absorber uid to absorbed uid
                        assert(absorberWalk<countNewNode);
                        assert(absorbedWalk<countNewNode);

                        ccAdjList[absorberWalk].push_back(absorbedWalk);
                        ccAdjList[absorbedWalk].push_back(absorberWalk);
                        //
                        
                        edge_t enew;
                        if(e.left==e.right){
                            enew.left = (e.left?false:true);
                            enew.right = (e.right?false:true);
                        }else{
                            enew.left = e.left;
                            enew.right = e.right;
                        }
                        
                        enew.toNode = uid;
                        
                        
                        
//                      if(e.left == nodeSign[alluid] && e.right == nodeSign[e.toNode] )
//                          absorbedCategory[e.toNode] = '1';
//                      if(e.left == nodeSign[alluid] && e.right != nodeSign[e.toNode] )
//                          absorbedCategory[e.toNode] = '2';
//                      if(e.left != nodeSign[alluid] && e.right != nodeSign[e.toNode] )
//                                         absorbedCategory[e.toNode] = '3';
//                      if(e.left != nodeSign[alluid] && e.right == nodeSign[e.toNode] )
//                                             absorbedCategory[e.toNode] = '4';
                        //2 and 3
                        if((enew.left == nodeSign[e.toNode] && enew.right != nodeSign[uid]) || (enew.left != nodeSign[e.toNode] && enew.right  != nodeSign[uid] )){
                            absorbGraph[e.toNode].push_back(enew);
                                                       isAnAbsorbedWalk[oldToNew[uid].finalWalkId]=true;
                        }else{
                            absorbGraph[e.toNode].push_back(enew);
                            isAnAbsorbedWalk[oldToNew[uid].finalWalkId]=true;
                        }
                                                
                                          
                        
                        
                        
                    }
                }
            }
        }
        
        if(1==0){
            if(tup_i + 1 < sorter.size()){
                if(get<1>(sorter[tup_i + 1])!=finalWalkId && lastWalkStartingIndex!=tup_i){   //so it is the end of a walk, and it is not an isolated walk (with one vertex only)
                    vector<edge_t> adju = adjList.at(uid);
                    for (edge_t e : adju) {
                        int absorberWalk = oldToNew[e.toNode].finalWalkId;
                        int absorbedWalk = nd.finalWalkId;
                        if(absorberWalk != absorbedWalk){
                            //absorbGraph[e.toNode].push_back(enew);
                        }
                    }
                }
            }
        }
        prevWalkId = finalWalkId;
    }

    
    if(1==0){
        stack<int> stack;
        
          
          
           
           for(int ii=0;ii<countNewNode;ii++)
           {
               
                
             if(ccVisited[ii] == false &&  !obsoleteWalkId[ii])
             {
                 //cout<<endl;
                 stack.push(ii);
                 //vector<int>nodes;
                 //nodes.push_back(i);
                absorbGraphNumCC_endAbosrb++;
                int count=0;
                 
                 while (!stack.empty())
                                 {
                                     // Pop a vertex from stack and print it
                                     int s = stack.top();
                                     stack.pop();
                               
                                     // Stack may contain same vertex twice. So
                                     // we need to print the popped item only
                                     // if it is not visited.
                                     if (!ccVisited[s])
                                     {
                                         //cout << s << " ";
                                         ccVisited[s] = true;
                                     }
                               
                                     // Get all adjacent vertices of the popped vertex s
                                     // If a adjacent has not been visited, then push it
                                     // to the stack.
                                     for (auto i = ccAdjList[s].begin(); i != ccAdjList[s].end(); ++i)
                                         if (!ccVisited[*i])
                                             stack.push(*i);
                                 }
                //ccDFS(i,ccVisited,count, ccAdjList);
                //cout<<"This component has "<<count<<" nodes"<<"\n";
       //          for(int nn:nodes){
       //              cout<<nn<<",";
       //          }
                 //cout<<endl;
             }
               
           }
           cout<<"number of connected components (end): "<<absorbGraphNumCC_endAbosrb<<endl;
        delete [] ccVisited;

    }
       
    

}

//not  tested
vector<int> getAllUidsInWalk(vector<MyTypes::fourtuple>& sorter, map<int, int >& sorterIndexMap, int walkId){
    vector<int> uids;
    int it = sorterIndexMap[walkId];
    for(int i=it; i<sorter.size(); i++){
        int wid = get<1>(sorter[i]);
        if(wid!=walkId){
            break;
        }
    
        int uid = get<0>(sorter[i]);
        uids.push_back(uid);
    }
    
    return uids;
}



void absorbGraphIsCyclicUtil(int uid, char*& visited, map<int, vector<edge_t> >& absorbGraph, map<int, stack<edge_t> >& absorbGraphCycleRemoved, queue<int>& orderOfUnitigs, char*& absorbedCategory, vector<MyTypes::fourtuple>& sorter, map<int, int >& sorterIndexMap, char*& visitedUID, int depth )
{
    
    /**
     map<int, int > sorterIndexMap: map walkID to location of <walk start, end> in sorter
     //v = visited[oldToNew[uid].finalWalkId] is the absorber
     ****/
     
    //cout<<oldToNew[uid].finalWalkId<<"("<<uid<<"),";
    orderOfUnitigs.push(uid);
    visited[oldToNew[uid].finalWalkId] = 'g';
    visitedUID[uid]= 'g';
    
//    if(depth>THEDEPTH){
//        visited[oldToNew[uid].finalWalkId] = 'b';
//        return;
//    }
    
    
    vector<int> uidsNeighbors = getAllUidsInWalk(sorter,sorterIndexMap, oldToNew[uid].finalWalkId);
    
    for(int alluid:uidsNeighbors){
        vector<edge_t> adjv = absorbGraph[alluid];
      
        
        for(edge_t e : adjv){
            int i = e.toNode;
            if (visited[oldToNew[i].finalWalkId] == 'g'){
                //removed edge
            }else if (visited[oldToNew[i].finalWalkId] == 'w'){
                //add to the list                //find absorption category
                //call the absorber v[i]
                //.....orderOfUnitigs.push(i);
                
                
                absorbGraphCycleRemoved[alluid].push(e);
                visited[oldToNew[alluid].finalWalkId] = 'g';
                
                
                if(e.left == nodeSign[alluid] && e.right == nodeSign[e.toNode] )
                    absorbedCategory[e.toNode] = '1';
                if(e.left == nodeSign[alluid] && e.right != nodeSign[e.toNode] )
                    absorbedCategory[e.toNode] = '2';
                if(e.left != nodeSign[alluid] && e.right != nodeSign[e.toNode] )
                                   absorbedCategory[e.toNode] = '3';
                if(e.left != nodeSign[alluid] && e.right == nodeSign[e.toNode] )
                                       absorbedCategory[e.toNode] = '4';
                
                
                absorbGraphIsCyclicUtil(i, visited, absorbGraph, absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory, sorter, sorterIndexMap, visitedUID, depth+1);
            }
        }
    
    }
    visited[oldToNew[uid].finalWalkId]  = 'b';
    return;
}

bool sortPairPos1(const pair<int, int> &lhs, const pair<int, int> &rhs){
    return get<1>(lhs) < get<1>(rhs);
}

void  removeCycleFromAbsorbGraph(vector<MyTypes::fourtuple>& sorter,  map<int, int >& sorterIndexMap, map<int, vector<edge_t> >& absorbGraph, map<int, stack<edge_t> >& absorbGraphCycleRemoved,  queue<int>& orderOfUnitigss, char*& absorbedCategory )
{
    
    //map<int, int> indegree;
    vector<pair<int, int> > indegree(countNewNode, make_pair(0,0));
    vector<int> oldnode(countNewNode, 0);
    //vector<pair<int, int> > indegree(adjList.size(), make_pair(0,0));
    for(int i = 0; i< indegree.size(); i++){
        indegree[i] = make_pair(i, 0);
    }
    
    //map<int, vector<edge_t> >& absorbGraph
    for (std::map<int, vector<edge_t> >::iterator it=absorbGraph.begin(); it!=absorbGraph.end(); ++it){
        int x =  it->first ;
        vector<edge_t> adjx = it->second;
        for(edge_t e: adjx){
                   int walk = oldToNew[e.toNode].finalWalkId;
                   oldnode[walk] =e.toNode;
                   indegree[walk] = make_pair(indegree[walk].first, indegree[walk].second+1);
               }
    }

    
//    for (const auto& [x, adjx]: absorbGraph) {
//
//    }
    //sort(indegree.begin(),indegree.end(),sortPairPos1);
    
    
   // int Vuid = absorbGraph.size();
   // int Vwalkid = sorterIndexMap.size();  // only the number of walks participating
    char *visited = new char[countNewNode];
    char *visitedUID = new char[adjList.size()];
    
    //
    //for(int i = 0; i < Vwalkid; i++)
    for(int i = 0; i < countNewNode; i++)
    {
        if(obsoleteWalkId[i]){
            visited[i] = 'b';
        }else{
             visited[i] = 'w';
        }
       
    }
    
    for(int i = 0; i < adjList.size(); i++){
        visitedUID[i]='w';
    }
    
    deque<int> dq;
    for(int ii = 0; ii <adjList.size(); ii++){
        if(indegree[oldToNew[ii].finalWalkId].second==0){
            dq.push_front(ii);
        }else{
            dq.push_back(ii);
        }
    }
  
    // Call the recursive helper function to detect cycle in different DFS trees
//for(int i = 0; i <adjList.size(); i++){
    while(!dq.empty()){
        int i = dq.front();
        dq.pop_front();
        if(absorbGraph.count(i)>0){
            //for(int i = 0; i <adjList.size(); i++){
            int watch = visited[oldToNew[i].finalWalkId] ;
            if(visited[oldToNew[i].finalWalkId] == 'w'){
                //orderOfUnitigss.push(i);
                //cout<<"dfs: "<<": ";
                absorbGraphIsCyclicUtil(i, visited, absorbGraph, absorbGraphCycleRemoved, orderOfUnitigss, absorbedCategory,sorter, sorterIndexMap, visitedUID, 0);
                //cout<<endl;
            }
        }
    
    }
    delete[] visited;
    delete[] visitedUID;
    //map int, priority queue uid loc
}

string recursiveWalkStringMaker(int& startWalkIndex, vector<bool>& isItAPrintedWalk, vector<MyTypes::fourtuple>& sorter, map<int, int >& sorterIndexMap,  map<int, stack<edge_t> >& absorbGraphCycleRemoved, queue<int>& orderOfUnitigs, char*& absorbedCategory, int depth, int Kpass){

    
    
    bool isThisAbsorbedWalk = false;;
    string unitigString;
    string walkString = "";
    
    int prevWalkId=-1;
    //int SMALLK = 20;
    
    while(true){
        assert(startWalkIndex<sorter.size());
        MyTypes::fourtuple n = sorter[startWalkIndex];
        int finalWalkId = get<1>(n);
        
        //depthfile<<get<0>(n)<<" "<<depth<<endl;
        
       //@ABSORB
        if( MODE_ABSORPTION_NOTIP ){
            if (isItAPrintedWalk[finalWalkId]) return "";
        }
        
      
       int uid = get<0>(n);
       int pos_in_walk = get<2>(n);
       int isTip = get<3>(n);
       //cout<<uid<<" " <<finalWalkId<<" "<<pos_in_walk<<" "<<isTip<<endl;
        
        
        //ALLABSORBING
//        bool isNextNeighborAbsorbed14 = false;
//        int neighborUid = -1;
//        if(startWalkIndex+1<sorter.size()){
//            neighborUid = get<0>(sorter[startWalkIndex+1]);
//            if(ALLABSORBING){
//                isNextNeighborAbsorbed14 = (absorbedCategory[neighborUid] == '4' or absorbedCategory[neighborUid] == '1');
//            }
//        }else{
//
//        }
        
        
        
        if(nodeSign[uid] == false){
            unitigString =  reverseComplement(unitigs.at(uid).sequence);
        }else{
            unitigString =  (unitigs.at(uid).sequence);
        }
        
        stack<edge_t> stType12;
        stack<edge_t> stType34;
        if(absorbGraphCycleRemoved.count(uid) > 0){         /*populate two types of stacks*/
            stack<edge_t> st = absorbGraphCycleRemoved[uid];
            while(!st.empty()){
                edge_t st_top = st.top();
                if(absorbedCategory[st_top.toNode]=='3' or absorbedCategory[st_top.toNode]=='4'  ){
                    stType34.push(st_top);
                }else if(absorbedCategory[st_top.toNode]=='1' or absorbedCategory[st_top.toNode]=='2' ){
                    stType12.push(st_top);
                }else{
                    cout<<"errorrrrrrrrrrr"<<endl;
                    assert(false);
                }
                st.pop();
            }
        }
        
        
        //walk starting AND not absorbed
        if(pos_in_walk == 1 && absorbedCategory[uid]=='0' && MODE_ABSORPTION_NOTIP){
            walkString+= cutSuf(unitigString, Kpass);

            while(!stType34.empty()){
                int sindex = sorterIndexMap[oldToNew[stType34.top().toNode].finalWalkId];
                char isSmallK =stType34.top().isSmallK;
                stType34.pop();
                
                if(isSmallK=='y'){
                    //walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory, depth+1, SMALLK);
                }else{
                    walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory, depth+1, K);
                }
                
            }
            
            walkString+= suf(unitigString, Kpass);
           
            
            
            while(!stType12.empty()){
                int sindex = sorterIndexMap[oldToNew[stType12.top().toNode].finalWalkId];
                char isSmallK =stType12.top().isSmallK;
                stType12.pop();
                if(isSmallK=='y'){
                    //walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory, depth+1, SMALLK);
                }else{
                    walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory, depth+1, K);
                }
            }
        }
        
        
        //not walk starting AND not absorbed
        if(pos_in_walk != 1 && absorbedCategory[uid]=='0' && MODE_ABSORPTION_NOTIP){
            while(!stType34.empty()){
                int sindex = sorterIndexMap[oldToNew[stType34.top().toNode].finalWalkId];
                char isSmallK =stType34.top().isSmallK;
                stType34.pop();
                
                if(isSmallK=='y'){
                    //walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory, depth+1, SMALLK);
                }else{
                    walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory, depth+1, K);
                }
            }
            
            walkString+= cutPref(unitigString, Kpass);
            
            
            
            while(!stType12.empty()){
                int sindex = sorterIndexMap[oldToNew[stType12.top().toNode].finalWalkId];
                char isSmallK =stType12.top().isSmallK;
                stType12.pop();
                if(isSmallK=='y'){
                    //walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory, depth+1, SMALLK);
                }else{
                    walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory, depth+1, K);
                }
            }
        }
        
        //not walk starting AND absorbed
//        if(pos_in_walk != 1 && absorbedCategory[uid]!='0' && MODE_ABSORPTION_NOTIP){
//            assert(false);
//        }
        
        if(pos_in_walk == 1 && absorbedCategory[uid]!='0' && MODE_ABSORPTION_NOTIP){
        //if(absorbedCategory[uid]!='0' && MODE_ABSORPTION_NOTIP){
            isThisAbsorbedWalk=true;
            if(absorbedCategory[uid]=='2' or absorbedCategory[uid]=='3'){
                walkString+=cutSuf(unitigString, Kpass);
            }
            if(absorbedCategory[uid]=='1') walkString+="+";
            else if(absorbedCategory[uid]=='4') walkString+="-";
            while(!stType34.empty()){
                int sindex = sorterIndexMap[oldToNew[stType34.top().toNode].finalWalkId];
                char isSmallK =stType34.top().isSmallK;
                stType34.pop();
                
                if(isSmallK=='y'){
                    //walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory, depth+1, SMALLK);
                }else{
                    walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory, depth+1, K);
                }
            }
            
            if(absorbedCategory[uid]=='1' or absorbedCategory[uid]=='4' ){
                walkString+= cutPref(unitigString, Kpass);
            }
            
            while(!stType12.empty()){
                int sindex = sorterIndexMap[oldToNew[stType12.top().toNode].finalWalkId];
                char isSmallK =stType12.top().isSmallK;
                stType12.pop();
                
                if(isSmallK=='y'){
                    //walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory, depth+1, SMALLK);
                }else{
                    walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory, depth+1, K);
                }
            }
            if(absorbedCategory[uid]=='3') walkString+="+";
            else if(absorbedCategory[uid]=='2') walkString+="-";
        }

        if(pos_in_walk != 1 && absorbedCategory[uid]!='0' && MODE_ABSORPTION_NOTIP){
        //if(absorbedCategory[uid]!='0' && MODE_ABSORPTION_NOTIP){
            isThisAbsorbedWalk=true;
            string part1="";
            string part2="";
            if(absorbedCategory[uid]=='2' or absorbedCategory[uid]=='3'){
                int indexFromLast;
                walkString+=cutSuf(walkString, Kpass);
                //assert(false);
               // walkString+=cutSuf(walkString, Kpass);
                //cout<<suf(walkString, Kpass)<<endl;
                //walkString=cutSuf(walkString, Kpass)+cutSuf(unitigString, Kpass);
                //cout<<walkString<<endl;
                
                
//                int cntACGT=0;
//                string part3= "";
//                string sufK1 = "";
//                stack<char> checkStack;
//                for(indexFromLast=walkString.length()-1; indexFromLast >= 0;  ){
//                    if(walkString[indexFromLast]=='A' || walkString[indexFromLast]=='C' || walkString[indexFromLast]=='G' || walkString[indexFromLast]=='T'){
//                        cntACGT++;
//                        sufK1=walkString[indexFromLast--]+sufK1;
//                        if(cntACGT==Kpass) break;
//                    }else{
//                        while(1){
//                            part3=walkString[indexFromLast]+part3;
//                            if(walkString[indexFromLast]==']'){
//                                checkStack.push(']');
//                            }else if(walkString[indexFromLast]=='[' and !checkStack.empty()){
//                                checkStack.pop();
//                            }
//                            indexFromLast--;
//                            if(checkStack.empty()){
//                                break;
//                            }
//                        }
//
//                    }
//                }
//                 part1= walkString.substr(0, indexFromLast+1);
//                 part2= walkString.substr(indexFromLast+1, walkString.length()-indexFromLast);
//                cout<<walkString<<endl;
//                cout<<part1<<endl;
//                cout<<part2<<endl;
//                cout<<endl;
//                walkString = part1;
                
                
            }
            if(absorbedCategory[uid]=='1') walkString=cutSuf(walkString, Kpass)+"+";
            else if(absorbedCategory[uid]=='4') walkString=cutSuf(walkString, Kpass)+"-";
            //walkString = walkString + cutSuf(part2, Kpass);
            while(!stType34.empty()){
                int sindex = sorterIndexMap[oldToNew[stType34.top().toNode].finalWalkId];
                char isSmallK =stType34.top().isSmallK;
                stType34.pop();
                
                if(isSmallK=='y'){
                    //walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory, depth+1, SMALLK);
                }else{
                    walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory, depth+1, K);
                }
            }
            
            if(absorbedCategory[uid]=='1' or absorbedCategory[uid]=='4' ){
                walkString+= cutPref(unitigString, Kpass);
            }
            
            while(!stType12.empty()){
                int sindex = sorterIndexMap[oldToNew[stType12.top().toNode].finalWalkId];
                char isSmallK =stType12.top().isSmallK;
                stType12.pop();
                
                if(isSmallK=='y'){
                    //walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory, depth+1, SMALLK);
                }else{
                    walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory, depth+1, K);
                }
            }
            if(absorbedCategory[uid]=='3') walkString+="+";
            else if(absorbedCategory[uid]=='2') walkString+="-";
        }
    
        
        if( MODE_ABSORPTION_TIP && isTip == 0){//
            //walkString+= cutPref(unitigString, K);
            walkString = plus_strings(walkString, unitigString, K);
            
        }else if(isTip==1 && MODE_ABSORPTION_TIP){ //right R   R    ]   ]   ]   ]
            //cut prefix: correct
            if(0==0){
                unitigString = unitigString.substr(K - 1, unitigString.length() - (K - 1));
                if(walkString.length()<K){
                    cout<<"pos: "<<walkString.length()<<endl;
                }
                walkString += "(" + unitigString + ")";
            }
            if(1==0){
                tipFile<<">pref\n"<<unitigString<<endl;
            }
            
        }else if(isTip==2 && MODE_ABSORPTION_TIP){ //left L   L    [ [ [
            //cut suffix: correct
            if(0==0){
                unitigString = unitigString.substr(0, unitigString.length() - (K - 1));
                if(walkString.length()<K){
                    cout<<"pos: "<<walkString.length()<<endl;
                }
                walkString += "{" + unitigString + "}";
            }
            if(1==0){
                tipFile<<">suf\n"<<unitigString<<endl;
            }
        }
            
        if(startWalkIndex+1 == sorter.size()) {
            //print previous walk
            tipDebugFile<<">"<<finalWalkId << " " << uid<<" " <<finalWalkId<<" "<<pos_in_walk<<" "<<isTip<<endl;
            // V_tip_ustitch++;
           // C_tip_ustitch+=walkString.length();
//
            tipDebugFile<<walkString<<endl;
          //  tipFile<< walkString<<endl;
            
            if(MODE_ABSORPTION_NOTIP){
                isItAPrintedWalk[finalWalkId] = true;
            }
            
            break;
        }else if(get<1>(sorter[startWalkIndex+1]) != finalWalkId){
                        //print previous walk
                        tipDebugFile<<">"<<finalWalkId << " " << uid<<" " <<finalWalkId<<" "<<pos_in_walk<<" "<<isTip<<endl;
            //           V_tip_ustitch++;
            //            C_tip_ustitch+=walkString.length();
            //
                        tipDebugFile<<walkString<<endl;
             //           tipFile<< walkString<<endl;
                        
                        isItAPrintedWalk[finalWalkId] = true;
                        break;
        }
        startWalkIndex++;
    }
    if(isThisAbsorbedWalk) walkString="["+walkString+"]";
    return walkString;
}

void tipAbsorbedOutputter(vector<MyTypes::fourtuple>& sorter, map<int, int >& sorterIndexMap, map<int, stack<edge_t> >& absorbGraphCycleRemoved, queue<int>& orderOfUnitigs, char*& absorbedCategory ){
    /// START OUTPUTTING
    
    
    
    vector<bool> isItAPrintedWalk;
    for(int i = 0; i< countNewNode; i++) isItAPrintedWalk.push_back(false);

    tipFile.open("tipOutput.txt");
    tipDebugFile.open("tipDebug.txt");
    
    int depth = 0;
    while(!orderOfUnitigs.empty()){
        int it = sorterIndexMap[oldToNew[orderOfUnitigs.front()].finalWalkId];
        
        if(DBGFLAG==PRINTER){
            cout<<"order: "<<orderOfUnitigs.front()<<endl;
        }
        
        
        orderOfUnitigs.pop();
        depth = 0;
        string walkString = recursiveWalkStringMaker(it, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory, depth, K);

        //tipDebugFile<<">"<<finalWalkId << " " << uid<<" " <<finalWalkId<<" "<<pos_in_walk<<" "<<isTip<<endl;
        
        if(walkString!="" && MODE_ABSORPTION_NOTIP){
            V_tip_ustitch++;
            C_tip_ustitch+=walkString.length();
            
           // tipDebugFile<<walkString<<endl;
            tipFile<< ">\n"<< walkString<<endl;
            C_oneabsorb+=walkString.length();
            V_oneabsorb++;
            
            int brackets1 = std::count(walkString.begin(), walkString.end(), '[');
            int brackets2 = std::count(walkString.begin(), walkString.end(), ']');
             int stringPlus = std::count(walkString.begin(), walkString.end(), '+');
               int stringMinus = std::count(walkString.begin(), walkString.end(), '-');
            
            C_oneabsorb_ACGT +=walkString.length() -brackets1- brackets2 -stringPlus-stringMinus;
            C_oneabsorb_plusminus += stringPlus +stringMinus;
            C_oneabsorb_brackets +=brackets1+brackets2;
        }
    }
    
    //exit(1);
    
    int it = -1;
    while(true){
        if(orderOfUnitigs.empty()){
            it++;
            if(it>=sorter.size()) break;
        }else{
            assert(false);
            //n = sorter[sorterIndexMap[orderOfUnitigs.front()]];
            it = sorterIndexMap[oldToNew[orderOfUnitigs.front()].finalWalkId];
            orderOfUnitigs.pop();
        }
        //recursive call
        depth = 0;
        string walkString = recursiveWalkStringMaker(it, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory, depth, K);

        //tipDebugFile<<">"<<finalWalkId << " " << uid<<" " <<finalWalkId<<" "<<pos_in_walk<<" "<<isTip<<endl;
        
        if(walkString!=""){
            V_tip_ustitch++;
            C_tip_ustitch+=walkString.length();
            
           // tipDebugFile<<walkString<<endl;
            tipFile<<">\n"<< walkString<<endl;
            
            C_oneabsorb+=walkString.length();
            
            
           
            int brackets1 = std::count(walkString.begin(), walkString.end(), '[');
            int brackets2 = std::count(walkString.begin(), walkString.end(), ']');
             int stringPlus = std::count(walkString.begin(), walkString.end(), '+');
               int stringMinus = std::count(walkString.begin(), walkString.end(), '-');
            
            C_oneabsorb_ACGT +=walkString.length() -brackets1- brackets2 -stringPlus-stringMinus;
            C_oneabsorb_plusminus += stringPlus +stringMinus;
            C_oneabsorb_brackets +=brackets1+brackets2;
            V_oneabsorb++;
        }
        //
    }
    //delete [] isItAnAbsorbedWalk;
    tipDebugFile.close();
    tipFile.close();
}

//map<int, vector<edge_t> >
void print_absorb_graph(map<int,vector<edge_t> > const &m)
{
    cout<<"ABSORB GRAPH---------------"<<endl;
    for (auto const& pair: m) {
        std::cout << pair.first << ":  ";
        vector<edge_t> vec = pair.second;
        for(int i=0; i<vec.size(); i++){
            edge_t e = vec[i];
            cout<< e.toNode  << ", ";
        }
        cout<<"\n";
    }
    cout<<"ABSORB GRAPH END---------------"<<endl;
}

//map<int, vector<edge_t> >
void print_absorb_graph_acyclic(map<int,stack<edge_t> > const &m)
{
    cout<<"ABSORB GRAPH ACYCLIC---------------"<<endl;
    for (auto const& pair: m) {
        std::cout << pair.first << ":  ";
        stack<edge_t> vec = pair.second;
        while(!vec.empty()){
            edge_t e = vec.top();
            vec.pop();
            cout<< e.toNode  << ", ";
        }
        cout<<"\n";
    }
    cout<<"ABSORB GRAPH ACYCLIC END---------------"<<endl;
}



void makeGraphDotAbsorb(map<int, vector<edge_t> > adjList){
    FILE * fp;
    
    fp = fopen ("/Users/Sherlock/dot-visualize/graph.gv", "w+");
    
    fprintf(fp, "digraph d {\n");
    set<int> vertices;
    for (std::map<int, vector<edge_t> >::iterator it=adjList.begin(); it!=adjList.end(); ++it){
    int x =  it->first ;
    vector<edge_t> adjX = it->second;
        for(edge_t ex: adjX){
            vertices.insert(oldToNew[x].finalWalkId);
             vertices.insert(oldToNew[ex.toNode].finalWalkId);
            fprintf(fp, "%d -> %d[taillabel=\"%d\", headlabel=\"%d\"]\n", oldToNew[x].finalWalkId, oldToNew[ex.toNode].finalWalkId, ex.left, !ex.right);
            cout<<oldToNew[x].finalWalkId<<"->"<<oldToNew[ex.toNode].finalWalkId<<endl;
            
            //            pair<int, int> p;
            //            if(x < ex.toNode){
            //                p.first = x;
            //                p.second = ex.toNode;
            //            }else{
            //                p.second = x;
            //                p.first = ex.toNode;
            //            }
            //            edges.insert(p);
        }
    }
    for(int i = 0 ; i<countNewNode; i++){
        if(!obsoleteWalkId[i] && vertices.count(i)==0 ){
            vertices.insert(i);
        }
    }
    for(int x: vertices){
        if(false){
            fprintf(fp, "%d [label=\"%d\", color=\"red\"]\n", x, x);
        }else{
            fprintf(fp, "%d [label=\"%d\"]\n", x, x);
        }
        
    }
    
    //for all int in list
    //make a list of neighbors add them
    
    
    fprintf(fp, "}\n");
    
    fclose(fp);
}


void absorptionManager(vector<MyTypes::fourtuple> sorter) {

    if(DBGFLAG==PRINTER){
        //    int uid = get<0>(n);
        //    int finalWalkId = get<1>(n);
        //    int pos_in_walk = get<2>(n);
        for(auto i: sorter ){
            cout<<"walk "<<get<1>(i)<<":";
            cout<<"("<<get<2>(i)<<")";
             cout<<get<0>(i)<<"";
            cout<<endl;
        }
    }
    
    int totuids = adjList.size();
    char* absorbedCategory = new char[adjList.size()];
    for(int i=0; i< totuids; i++){
        absorbedCategory[i] = '0';
        
        //hasherSmallerK(i);
    }
    
    
    //new mode small K here
//    for(int i=0; i< totuids; i++){
//        addSmallerKEdge(i);
//       }
//    
    
    map<int, int > sorterIndexMap;
    map<int, vector<edge_t> > absorbGraph;
    if(MODE_ABSORPTION_NOTIP){
        sorterIndexAndAbsorbGraphMaker(sorter, sorterIndexMap, absorbGraph, absorbedCategory);
    }
    
    //exit(2);
    if(DBGFLAG==PRINTER){
        print_absorb_graph(absorbGraph);
    }
    
    //makeGraphDotAbsorb(absorbGraph);
    
    map<int, stack<edge_t> > absorbGraphCycleRemoved;
    

    queue<int> orderOfUnitigs;
     if(MODE_ABSORPTION_NOTIP ){
    removeCycleFromAbsorbGraph(sorter,  sorterIndexMap, absorbGraph, absorbGraphCycleRemoved,  orderOfUnitigs,  absorbedCategory);//last two are output
         
    absorbGraph.clear();
     }
    cout<<"cycle removed from ab graph!"<<endl;
 
    
    if(DBGFLAG==PRINTER){
        for (int i=0; i<totuids; i++) {
            cout<<i<<"->"<<absorbedCategory[i]<<endl;
        }
         print_absorb_graph_acyclic(absorbGraphCycleRemoved);
    }
    tipAbsorbedOutputter(sorter, sorterIndexMap, absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory);
    
    cout<<"tip outputter done!"<<endl;
    delete [] absorbedCategory;
    
    depthfile.close();
}

#endif /* absorbGraph_hpp */
