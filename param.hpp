//
//  param.hpp
//  UST
//
//  Created by Amatur Rahman on 15/6/20.
//  Copyright © 2020 medvedevgroup. All rights reserved.
//

#ifndef param_h
#define param_h

#include<string>

using namespace std;
class Param{
public:
    bool VALIDATE = true;
    
//    string APP_PATH_DSK="bin/dsk";
//    string APP_PATH_DSK2ASCII="bin/dsk2ascii";
//    string APP_PATH_BCALM="bin/bcalm";
    
    string DSK_PATH="/Volumes/exFAT/work/dsk-v2.3.0-bin-Darwin/bin/";
    //string DSK_PATH="/home/aur1111/w/dsk/build/bin";
    ///home/aur1111/w/dsk/build/bin
    string UNITIG_FILE;
    string OUTPUT_FILENAME;
    
    string VERSION = "2.0";
    
};

class ESSTipParam:public Param{
public:
    string ofileTipOutput = "fa.esstip";
    
};

class ESSParam : public Param{
public:
    bool GETLOWERBOUND_CC = false;
    bool PROFILE_ONLY_ABS = false;
    bool EARLYABSORB = true;
    bool VALIDATE = true;
    
    /*FILENAMES*/
    string ofileTipOutput = "fa.ess";
};
#endif /* param_h */
