//
//  main.cpp
//  Prog6
//
//  Carlos Reyes
//  CSCI191T/ SPRING 2017
//  5/10/2017
/*************************************************************
 *              GLOBAL DOCUMENTATION
 * This program implements a branch and bound algorithm of
 * finding a median string and motif string by using an
 * L-mer tree. Using the branch and bound algorithm the time
 * complexity for obtaining the best median string is shortened
 * because we no longer search a space of 4^L. The algorithm
 * bypasses any nodes that are not within a range frame that
 * is suitable to be a good motif. The program will obtain the
 * five best median and motif strings and ouput the consensus
 * score/string along with the best positions and strings from
 * the sequences. The L-mers size will be a user a input size
 * and the file will also be given.
 *
 ************************************************************
 *              How to compile and run
 * This program was written in C++. In order to run the program through the
 * terminal, first you must compile the program using this command line
 * g++ main.cpp -o main.out
 * Once the file has been compiled, it can then be executed using the
 * command line
 * ./main.out HMP-part.fa > Prog6_sol.txt
 *            ^^^^^^^^^^^   ^^^^^^^^^^^^^
 *                  ^             ^
 *                  |             |
 *        (NameOfFile) (NameOfOutputfile)
 * The database sequence file must be the first parameter.
 * The second parameter that is in front of '>' is the output file which
 * the solution will be written to.
 *
 ****************************************************************************
*/

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <sstream>
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
// stuct that will hold all median information
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
int total_dist(string v, unordered_map<string,string> lns, int l){
    string seq;
    int totalH = 0;
    
    for(auto it : lns)
    {
        int hammingDist =999999;
        seq = it.second;
        for(int i=0; i<=seq.size()-l+1 && hammingDist !=0; i++){
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

//takes in the v string, level,l_mers, and 4
// the 4 is for the letter possibilites
//returns new v string and updates level
//works with leaf nodes but skips useless subtrees
string bypass(string s, int &i, int l, int k){
    for(int j = i; j >= 0; j--){
        if(s[j]-48 < k){
            s[j] +=1;
            i=j+1;
            return s;
        }
        s[j] = '1';
    }
    i =0;
    return s;
}
//takes in the v string, level, l-mer, and 4
// the 4 is for the letter possibilites
//iterates through through the tree nodes
//traverses through left child and moves on
//returns new v string and updates the level
string next_vertex(string s, int &i, int l, int k){
    if(i<l){
        s[i] = '1';
        i += 1;
        return s;
    }
    else{
        for(int j = l; j >=0; j--){
            if(s[j] -48 < k){
                s[j] += 1;
                i = j+1;
                return s;
            }
            s[j] ='1';
        }
    }
    i = 0;
    return s;
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

//implements the branch and bound algorithm
//Will ask for user input for l-mer
//Will open the database sequence file
//stores it into unordered map
//then begins to implements the branch and bound
//will call upon all other functions
int main(int argc, const char * argv[]) {
    int l_mer;
    int level;
    int pos=0;
    int optimistic_dist =0;
    string v_string;
    string prefix_str;
    unordered_map<string, string> labelAndSeq;  //unordered_map for labels and sequences
    cout<<"Input L_mers length size: "<<endl;
    cin >> l_mer; //keyboard input
    --l_mer;    //for base 0 purposes
    level = l_mer;
    ifs.open(argv[1]);
    //open file and store all labels and sequences into unordered_map until end of file
    while(!ifs.eof()){
        string seq;
        string label;
        getline(ifs,label);
        getline(ifs,seq);
        labelAndSeq[label] = seq;
        numOfSeq++;
    }
    ifs.close();    //close file
    //creates prefix string size and initializes v_string
    for(int i=0; i<=l_mer; i++){
        v_string += "1";
        prefix_str +=" ";
    }
    
    while(level>0){
        if(level<l_mer){
            for(int i = 0; i<=l_mer; i++){
                if(v_string[i] == '1'){
                    prefix_str[i]='A';
                }
                else if(v_string[i] == '2'){
                    prefix_str[i]='C';
                }
                else if(v_string[i] == '3'){
                    prefix_str[i]='G';
                }
                else if(v_string[i] == '4'){
                    prefix_str[i]='T';
                }
            }
            
            if(med[pos].med_Str ==""){
                while(med[pos].med_Str !="" && pos != 4){
                    ++pos;
                }
            }
            
            optimistic_dist =  total_dist(prefix_str, labelAndSeq, level);
            if(optimistic_dist >= med[pos].t_d){
                v_string = bypass(v_string, level, l_mer, 4);
            }
            else{
                v_string = next_vertex(v_string, level, l_mer, 4);
            }
        }
        else{
            for(int i = 0; i<=l_mer; i++){
                if(v_string[i] == '1'){
                    prefix_str[i]='A';
                }
                else if(v_string[i] == '2'){
                    prefix_str[i]='C';
                }
                else if(v_string[i] == '3'){
                    prefix_str[i]='G';
                }
                else if(v_string[i] == '4'){
                    prefix_str[i]='T';
                }
            }
            
            int tot_dist = total_dist(prefix_str, labelAndSeq, level);
            
            cout<<prefix_str<<endl;
            
            if(med[pos].med_Str ==""){
                while(med[pos].med_Str !="" && pos != 4){
                    ++pos;
                }
            }
            
            //update priority queue if total_distance > worst queue td
            string holdV= prefix_str;
            string temp_v="";
            //priority queue, keep best 5
            for(int p = 0; p < 5;p++){
                int temp_td =0;
                if(med[p].med_Str == "" ){
                    med[p].med_Str = holdV;
                    med[p].t_d = tot_dist;
                    break;
                }
                else if(med[p].t_d>tot_dist){
                    temp_td = med[p].t_d;
                    temp_v = med[p].med_Str;
                    med[p].t_d = tot_dist;
                    med[p].med_Str = holdV;
                    tot_dist = temp_td;
                    holdV = temp_v;
                }
                else if(holdV==temp_v && med[p].t_d==tot_dist){
                    temp_td = med[p].t_d;
                    temp_v = med[p].med_Str;
                    med[p].t_d = tot_dist;
                    med[p].med_Str = holdV;
                    tot_dist = temp_td;
                    holdV = temp_v;
                }
            }
            
            v_string = next_vertex(v_string, level, l_mer, 4);
        }
        
    }
    cout<<endl;
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
    return 0;
}
