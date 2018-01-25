/*
 * CSCI 191T
 * T,TH 2:00PM
 * ASSIGNMENT 3
 * 3/13/2017
 *
 *******************************************************************************
 *                          GLOBAL DOCUMENTATION
 *
 * This program finds the median string given DNA sequences and the length of
 * the motif. The program will find the 5 best median strings and store them
 * into an array using a priority queue method. with each median string the
 * array will also store the total distance, motif consensus string and
 * position as well as the consensus score. The program starts off
 * opening a text file given a file path and then begins to read the file and
 * store the label and sequences into an unordered map. The program will then to
 * generate all l-mers by using a search space of 4^L, where the L is a length of
 * 6 from AAAAAA~TTTTTT. As every l-mers is generetad we begin to find the total
 * distance and store the total distance and l-mers string into the struct array
 * of 5 using the priority queue method. Once all l-mers total distances have
 * been computed as well as there total distances and the 5 best l-mers has
 * been stored, the program will then begin to compute and store the motif
 * consensus position and string for all five median strings. While the motif
 * consensus strings is being stored the program will generete a profile and
 * store/update the profile using a 2D array. Once all motif consensus strings
 * have been generated and the profile is up to date then the consensus score
 * will be computed and stored into  the median struct with respectively with
 * the median string. Once all information for the five median strings are
 * created, then the function will return the struct array and display all 
 * the information along with the runtime output.
 
 *******************************************************************************
 
                        How to compile and run
 * This program was written in C++ using xCode. In order to run this program, 
 * the file path must must first be updated on line 241. Once the file path 
 *is updated then you can simply run the program using xcode. If you would like 
 * to run this program using the command line on terminal
 * then you will have to comment out line 242 and uncomment line 2.., then 
 * you can use
 ************* g++ CSCI191T_Assign3.cpp < textfile ***********************
 * Also your text file, and program should be in the same directory and then 
 *your terminal should be in the directory of your files in order to run.
 *******************************************************************************
 *                      Program Code Begins Here
*/


#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <queue>

using std::unordered_map;
using namespace std;

ifstream ifs; //input file stream
int numOfSeq;   //number of sequences stored

//struct which will hold the motif string and position
struct motif{
    string m_Str;   //motif string
    int m_pos;    //motif position
};

//struct which will hold all information
//holds, medain string, total distance,
//motif consensus score/string
//vector of motif positions and strings
struct median {
    string med_Str;  //median string
    string consStr; //consensus string
    int consScore;  //consunsus score
    int t_d;        //total_distance
    vector<motif> s_t;    //vector of motif positions/string
}med[5];

//findTotalDist() contains three parameters,
//the 'v' string, unordered map which contains labels and sequences from file,
//the size of our V(l-mer)
//function will return the total distance which is a sum of all the total
//hamming distances
int findTotalDist(string v, unordered_map<string,string> lns, int s){
    string seq;
    int totalH = 0;
    
    for(auto it : lns)
    {
        int hammingDist =999999;
        seq = it.second;
        for(int i=0; i<=seq.size()-s+1 && hammingDist !=0; i++){
            int hCounter = 0;
            if(v[0] != seq[i]){
                ++hCounter;
            }
            if(v[1] != seq[i+1]){
                ++hCounter;
            }
            if(v[2] != seq[i+2]){
                ++hCounter;
            }
            if(v[3] != seq[i+3]){
                ++hCounter;
            }
            if(v[4] != seq[i+4]){
                ++hCounter;
            }
            if(v[5] != seq[i+5]){
                ++hCounter;
            }
            if(hammingDist>hCounter){
                hammingDist = hCounter;
            }
        }
        totalH += hammingDist;
    }
    return totalH;
}

//findMotif() takes in three parameters,
// an array med_five holds the best five median strings
// unordered_map lands holds the labels and sequences
// sz is the size of the l-mers
//function will generate motif strings/positions from all sequences
//will compute a profile matrix in order to accumulate the
// consensus score and string
//will return a motif struct array
void findMotif(median med_five [],unordered_map<string, string> lands, int sz){
    string seq;
    string v;
    string sStr;
    string lst;
    int sPos=0;
    int profile [4][sz];    //matrix profile
    vector<motif> st;
    
    //loop through the 5 best l-mers
    for(int i=0; i<5; i++){
        //create 2d array with all counts as 0
        for(int q=0; q < 4; q++){
            for(int p=0;p<sz; p++){
                profile[q][p] =0;
            }
        }//end of 2d array initialization
        
        v=med_five[i].med_Str;
        int cnt=0;
        //loop through all sequences to find best positions and string
        for(auto it : lands){
            int hammingDist =999999;
            seq = it.second;
            //loop through all possible starting positions in from 0~n-l+1
            for(int i=0; i<=seq.size()-sz+1 && hammingDist !=0; i++){
                int hCounter = 0;\
                if(v[0] != seq[i]){
                    ++hCounter;
                }
                if(v[1] != seq[i+1]){
                ++hCounter;
                }
                if(v[2] != seq[i+2]){
                ++hCounter;
                }
                if(v[3] != seq[i+3]){
                ++hCounter;
                }
                if(v[4] != seq[i+4]){
                    ++hCounter;
                }
                if(v[5] != seq[i+5]){
                    ++hCounter;
                }
                //if found better hamming distance, store motif position and string
                if(hammingDist>hCounter){
                    hammingDist = hCounter;
                    sPos = i;
                    sStr= seq.substr(i,sz);
                }
            }//end of one sequence loop
            //store the best fit motif position and string
            med[i].s_t.push_back(motif());
            med[i].s_t[cnt].m_Str=sStr;
            med[i].s_t[cnt++].m_pos = sPos;
            
            //update profile matrix for every hamming distance string
            for(int k =0; k<sz; k++){
                if(sStr[k] == 'A'){
                    profile[0][k] +=1;
                }
                else if(sStr[k] == 'C'){
                    profile[1][k] +=1;
                }
                else if(sStr[k] == 'G'){
                    profile[2][k] +=1;
                }
                else if(sStr[k] == 'T'){
                    profile[3][k] +=1;
                }
            }
        }//end of all sequences loop
        
        int a=0;
        int c=0;
        int g=0;
        int t=0;
        
        //acccumulate the motif consesus score and string
        //store respectively into the pq struct of best 5 l-mers
        for(int k = 0; k<sz; k++){
            a=profile[0][k];
            c=profile[1][k];
            g=profile[2][k];
            t=profile[3][k];
            if(a>=c&&a>=g&&a>=t){
                med[i].consStr += 'A';
                med[i].consScore += a;
            }
            else if(c>=a&&c>=g&&c>=t){
                med[i].consStr += 'C';
                med[i].consScore += c;
            }
            else if(g>=a&&g>=c&&g>=t){
                med[i].consStr += 'G';
                med[i].consScore += g;
            }
            else if(t>=a&&t>=c&&t>=g){
                med[i].consStr += 'T';
                med[i].consScore += t;
            }
            
        }
    }//end of all 5 best l-mers loop
    return;
}

int main(int argc, const char * argv[]) {
    unordered_map<string, string> labelAndSeq;  //unordered_map for labels and sequences
    string vStr= "AAAAAA";  //starting l-mers string
    int total_dist= 0;      //total distance
    ifs.open("/Users/carlos/Desktop/csci191/HMP-part.fa");
    //ifs.open(argv[1]);
    //open file and store all labels and sequences into unordered_map until end of file
    while(!ifs.eof()){
        string seq;
        string label;
        getline(ifs,label);
        getline(ifs,seq);
        labelAndSeq[label] = seq;
        numOfSeq++;
    }
    ifs.close();
    //generate all l-mers from AAAAAA-TTTTTT
    for(int i = 0; i<4; i++){
        if(i == 0)
            vStr[0] ='A';
        else if(i == 1)
            vStr[0]='C';
        else if(i == 2)
            vStr[0]='G';
        else if (i ==3)
            vStr[0]='T';
        
        for(int j = 0; j<4; j++){
            if(j == 0)
                vStr[1] ='A';
            else if(j == 1)
                vStr[1]='C';
            else if(j ==2)
                vStr[1]='G';
            else if (j == 3)
                vStr[1]='T';
            for(int k = 0; k<4; k++){
                if(k == 0)
                    vStr[2] ='A';
                else if(k == 1)
                    vStr[2]='C';
                else if(k ==2)
                    vStr[2]='G';
                else if (k == 3)
                    vStr[2]='T';
                for(int l = 0; l<4; l++){
                    if(l == 0)
                        vStr[3] ='A';
                    else if(l == 1)
                        vStr[3]='C';
                    else if(l ==2)
                        vStr[3]='G';
                    else if(l == 3)
                        vStr[3]='T';
                    for(int m = 0; m<4; m++){
                        if(m == 0)
                            vStr[4] ='A';
                        else if(m == 1)
                            vStr[4]='C';
                        else if(m ==2)
                            vStr[4]='G';
                        else if(m == 3)
                            vStr[4]='T';
                        
                        for(int n = 0; n<4; n++){
                            if(n == 0)
                                vStr[5] ='A';
                            else if(n == 1)
                                vStr[5]='C';
                            else if(n ==2)
                                vStr[5]='G';
                            else if(n ==3)
                                vStr[5]='T';
                            //call and store total distance for given V(l-mer)
                            total_dist = findTotalDist(vStr, labelAndSeq,6);
                            string holdV= vStr;
                            string temp_v="";
                            //priority queue, keep best 5
                            for(int p = 0; p < 5 && med[5].t_d != total_dist;p++){
                                int temp_td =0;
                                if(med[p].med_Str == "" ){
                                    med[p].med_Str = holdV;
                                    med[p].t_d = total_dist;
                                    break;
                               }
                                else if(med[p].t_d>total_dist){
                                    temp_td = med[p].t_d;
                                    temp_v = med[p].med_Str;
                                    med[p].t_d = total_dist;
                                    med[p].med_Str = holdV;
                                    total_dist = temp_td;
                                    holdV = temp_v;
                                }
                                else if(holdV==temp_v && med[p].t_d==total_dist){
                                    temp_td = med[p].t_d;
                                    temp_v = med[p].med_Str;
                                    med[p].t_d = total_dist;
                                    med[p].med_Str = holdV;
                                    total_dist = temp_td;
                                    holdV = temp_v;
                                }
                            }
                        }
                    }
                }
            }
        }
    }//end of l-mers
    
    //call and store the motif score, string, positions/strings for best 5 l-mers
    findMotif(med,labelAndSeq,6);
    
    //display all info such as
    //best 5 l0-mers, consensus score/string, total distance, motif position/string
    for(int i=0; i< 5; i++){
        cout<<"median string: "<<med[i].med_Str<<"(tot_dist= "<<med[i].t_d<<")"<<endl;
        cout<<"motif consensus string: "<<med[i].consStr<<"(consensus_score = "<<med[i].consScore<<")"<<endl;
        cout<<"motif positions/string s=(s1..st): "<<endl;
        for(int p=0; p<numOfSeq; p++){
            cout<<med[i].s_t[p].m_pos<<"("<<med[i].s_t[p].m_Str<<")  ";
        }
        cout<<endl<<endl;
    }
    
    //display run-time output
    cout<<"run time in seconds: "<<clock()/ CLOCKS_PER_SEC<<endl;
    return 0;
}
