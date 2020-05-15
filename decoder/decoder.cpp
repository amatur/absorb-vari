//
//  decoder.h
//  bcl
//
//  Created by Amatur Rahman on 28/11/19.
//  Copyright © 2019 psu. All rights reserved.
//

//#ifndef decoder_h
//#define decoder_h
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
//

#include "decoder.hpp"

// string GetStdoutFromCommand(string cmd) {
//
//     string data;
//     FILE * stream;
//     const int max_buffer = 256;
//     char buffer[max_buffer];
//     cmd.append(" 2>&1");
//
//     stream = popen(cmd.c_str(), "r");
//     if (stream) {
//     while (!feof(stream))
//     if (fgets(buffer, max_buffer, stream) != NULL) data.append(buffer);
//     pclose(stream);
//     }
//     return data;
// }


int main(int argc, char** argv) {
    // int pid = (int) getpid();
    // cout<<pid<<endl;

    const char* nvalue = "" ;
    int K=0;
    string filename = "ust_ess_abs.txt";
    bool tip_mode = 0;



        int c ;

            while( ( c = getopt (argc, argv, "i:k:t:") ) != -1 )
            {
                switch(c)
                {
                    case 'i':
                        if(optarg) nvalue = optarg;
                        break;
                    case 't':
                        if(optarg) {
                            tip_mode = static_cast<bool>(std::atoi(optarg));
                        }
                        break;
                    case 'k':
                        if(optarg) {
                            K = std::atoi(optarg) ;
                            if(K<=0){
                                fprintf(stderr, "Error: Specify a positive k value.\n");
                                exit(EXIT_FAILURE);
                            }
                        }else{
                            fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name>\n",
                                    argv[0]);
                            exit(EXIT_FAILURE);
                        }
                        break;
                    default: //
                        fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-path>\n",
                                argv[0]);
                        exit(EXIT_FAILURE);

                }
            }


    if(K==0 || strcmp(nvalue, "")==0){
        fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name>\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if(tip_mode){
        filename = "ust_ess_tip.txt";
    }

    string pathname = string(nvalue);
    pathname+="/"+filename;
    //int K=31, string ENCODED_FILE= "tipOutput.txt"

    if(tip_mode){
        decodeTip(K, pathname);
        cout<<"UST-TIP decoding done!"<<endl;
    }else{
        decodeOneAbsorb(K, pathname);
        cout<<"UST-ABS decoding done!"<<endl;
    }


    // cout<<GetStdoutFromCommand("top -p "+ std::to_string(pid)+" -b -n1 | tail -1 |  awk '{ print $5}'")<<endl;
    // cout<<"printed memory requirement"<<endl;
    //

    return 0;
}



//#endif /* decoder_h */
