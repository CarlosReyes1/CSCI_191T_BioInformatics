 //
// Carlos Reyes
// SID: 108939200
// CSCI 191T
// PROGRAM 5
// 4/24/2017
//
/************************************************************************
 *                      Global Documentation
 * This short program checks the mappability of gene sequences which falls
 * within its scope. In this program three different read lengths are
 * implemented. The three read lengths are 50,70,100. Upon reading
 * reading the three read lengths with the corresponding gene ID from
 * the gene sequences, the percentage of mappability is provided.
 * The mappibility is computed by taking the length of gene sequence
 * using the end position minus the start position which is n. Then
 * we use n-l+1 to produce the total number of tiles(tot), where l is the
 * the read length. Once tot is computed we can divide it with the read count
 * which is accumulated by adding the number of times a gene ID appears
 * in the reads file that are within the scope, and multiply it by 100
 * to gives the percentage.
 *************************************************************************
 *                  How to run the program
 * This program was written in C++. In order to run the program through the
 * terminal, first you must compile the program using this command line
 * g++ main.cpp -o main.out
 * Once the file has been compiled, it can then be executed using the
 * command line
 * ./main.out reads50.fa HG38-refseq-annot-chr1-250 > Prog5_sol.txt
 *            ^^^^^^^^^^ ^^^^^^^^^^^^^^^^^^^^^^^^^    ^^^^^^^^^^^^^
 *                  ^             ^                   ^
 *                  |             |                   |
 *        (NameOfReadsFile) (NameOfAnnotfile)     (NameOfOutputfile)
 * The reads file must be the first parameter and the annotation file must be
 * the second paramater in order for this program to function properly.
 * The third parameter that is in front of '>' is the output file which
 * will the solution will be written to.
 *
 ****************************************************************************
 *                      Program Starts Here
*/

#include <iostream>
#include <string>
#include <unordered_map>
#include <fstream>
#include <queue>
#include <iomanip>

using std::unordered_map;
using namespace std;

ifstream ifr; //input file reads
ifstream ifga; //input file gene annotation

//opens the reads file and stores the unique gene ID and number of
//appearances that fall within its scope
//Once the entire file has been read then the file closes
//Then the annotation file is opened three compenents are
//read from the file. The first two components are the
//starting and ending positions. The third component is the
//geneID. The geneID is then searched throughout the
//unorded map created. If not found the geneID will print
//along with a zero for the mappibility. If the geneID
//is found then the mappability is computed using this
//equation read_count(from unorded map) / TOT(total number
//of tiles) *100 and the geneID and mappibility percentage
//is displayed. Will do this till the end of file is reached
//Once complete the annotation file will be closed and program
//will terminate.
int main(int argc, const char * argv[]) {
    //unordered map for the reads file
    //key is string for unique gene ID
    //value is an int to count number of appearances of gene ID
    unordered_map<string , int> mymap;
    string seq;
    ifr.open("/Users/carlos/Desktop/csci191/Prog5/reads50.fa");
    ifga.open(argv[2]);
    string gene;
    string chr = "";
    string endPos = "";
    string startPos = "";
    string line = "";
    string rstart="";
    string rend="";
    string nID="";
    double mappability;

    while(!(ifr.eof())){
        int period=0;
        string readID = "";
        string geneID = "";
        chr ="";
        startPos="";
        endPos="";
        nID="";
        rend="";
        rstart="";
        line ="";
        ifr>>chr;
        ifr>>startPos;
        ifr>>endPos;
        ifr>>readID;
        getline(ifr, line);
        for(int i =5; i<readID.length() && period <3;i++){
            geneID = readID.substr(0,5);
            if(readID[i]=='.'){
                ++period;
                continue;
            }
            if(period ==0){
                rstart += readID[i];
            }
            else if(period ==1)
                rend += readID[i];
            else
                nID += readID[i];
        }
        
        geneID += rstart+"."+rend+"."+nID;
        int sp = stoi(startPos);
        int rs =stoi(rstart);
        int ep =stoi(endPos);
        int re = stoi(rend);
        
        if(( sp >= rs) && (ep <= re)){
            mymap[geneID]++;
        }
    }
    ifr.close();
    
    cout<<setw(20)<<left<<"geneID"<<setw(20)<<"readLen50"<<endl;
    
    //cout<<setw(20)<<left<<"geneID"<<setw(20)<<"readLen70"<<endl;
    
    //cout<<setw(20)<<left<<"geneID"<<setw(20)<<"readLen100"<<endl;
    while(!(ifga.eof())){
        chr ="";
        startPos="";
        endPos="";
        nID="";
        line="";
        int read_count=0;
        ifga>>chr;
        ifga>>startPos;
        ifga>>endPos;
        ifga>>nID;
        getline(ifga, line);
        int n = stoi(endPos)-stoi(startPos);
        string s2 = nID;
        for (unordered_map<string,int>::const_iterator it = mymap.begin(); it != mymap.end(); ++it ){
            string s1 = it->first;
            if (s1.find(s2) != string::npos) {
                read_count = it->second;
                break;
            }
        }
        
        cout<<setw(20)<<left<<nID;
        if (read_count==0) {
            cout<<right<<setw(9)<<0<<endl;
        }
        else if(n-50+1<=0){
            cout<<right<<setw(9)<<1<<endl;
        }
        else{
            mappability = (double(read_count)/double(n-50+1))*100;
            cout<<right<<setw(9)<<setprecision(4)<<mappability<<endl;
        }
        
        /*
         if (read_count==0) {
         cout<<setw(20)<<0;
         }
         else{
         mappability = double(read_count)/double(n-70+1)*100;
         cout<<setw(20)<<setprecision(4)<<mappability;
         }
         */
        
        /*
         if (read_count==0) {
         cout<<setw(20)<<0<<endl;
         }
         else{
         mappability = double(read_count)/double(n-100+1)*100;
         cout<<setw(20)<<setprecision(4)<<mappability<<endl;
         }
         */
    }
    ifga.close();
    return 0;
}
