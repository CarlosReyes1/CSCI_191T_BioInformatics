 /* Carlos Reyes
 * CSCI 191T
 * M,W 12:30
 * FEB. 14, 2017
 * Program 2
 *****************************************************************************
                `       Global Documentation
 This program provides all possible solutions to a partial digest problem via the 
 practical algorithm. Given a list, this program will output the possible 
 solutions. This program too, will retrieve a set of pair wise distances from a 
 file and produce DNA restriction sites. This program is written in C++ using 
 XCODE as the IDE, which will compile and run the program. However, to run a 
 different file you must change the file path in line 194 . The file being 
 given must be formatted with spaces between number, for example "1 1 1 2 2 2."
 
 *****************************************************************************
*/
                       //Program begins here
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <iostream>
#include <stdlib.h>

using namespace std;

vector<int> list1;  //L
vector<int> list2;  //X
int m;  //max
ifstream ifs;

int find_max(vector<int>); void Place(vector<int>, vector <int>);

//Recursive funtion in order to compute, store, and display sets of the
//restriction sites from given a set of pair-wise distances of restriction
//sites. Once list1 is empty then we hit the recursion terminal case
//In Place() we will compute all the distances using two different cases.
//The first case is using max first.
//The second case is using Max - y(new max in list1)
void Place(vector<int> l, vector<int> x){
    vector<int> tempList1;  //temp vector for list1
    vector<int> dist_vec;   //dist vector
    int y;    //y is max in list1
    int m_y;
    int dist;
    bool all_found_flag;
    int pos;
    
    //terminal condition for recursion
    if(l.empty()){
        sort(x.begin(), x.end());
        cout<<"X = {";
        for(int i =0; i< x.size(); i++){
            if(i==x.size()-1){
                cout<<x[i]<<"}"<<endl;
            }
            else
                cout<<x[i]<<", ";
        }
        return;
    }
    //y is max of list1
    y = find_max(l);
    
    //build distance vector for list2
    for(int i = 0; i<x.size(); i++){
        dist= (int)abs(y-x[i]);     //distance computed
        dist_vec.push_back(dist);   //push dist to distance vector
    }
    //copy list1 to temp list
    tempList1.resize(l.size());
    for(int i =0; i< l.size(); i++){
        tempList1[i]= l[i];
    }
    all_found_flag=1;   //default true
    
    for(int j=0; j< dist_vec.size();j++){
        vector <int>::iterator i = tempList1.begin ();
        i = find (tempList1.begin (),tempList1.end (), dist_vec[j]);
        pos = (int)distance(tempList1.begin(), i);
        
        if(dist_vec[j] != tempList1[pos]){
            all_found_flag=0;
            break;
        }
        else{
            tempList1.erase(tempList1.begin()+pos);
        }
    }//end of j=0~dist_vec.size
    if(all_found_flag){
        //add y to list2
        x.push_back(y);
        
        //copy tempList to list1
        l.resize(tempList1.size());
        for(int k = 0; k < tempList1.size(); k++){
            l[k] = tempList1[k];
        }
        //recursive call with list1, list2
        Place(l, x);
        
        //pos<-index of y in list2
        vector <int>::iterator it = x.begin();
        it = find (x.begin (),x.end (), y);
        pos = (int)distance(x.begin(), it);
        
        //delete y from list2
        x.erase(x.begin()+pos);
        
        //restore list1
        for(int i =0; i<dist_vec.size(); i++){
            l.push_back(dist_vec[i]);
        }
    }//end of if(all_found_flag == true)
    
    //flush tempList and dist_vect
    tempList1.resize(0);
    dist_vec.resize(0);
    m_y = m-y;
    
    //build dist_vec from list2(x)
    for(int j=0; j<x.size(); j++){
        dist = abs(m_y - x[j]);
        dist_vec.push_back(dist);
    }
    //copy list1 to tempList
    tempList1.resize(l.size());
    for(int k=0; k<l.size(); k++){
        tempList1[k] = l[k];
    }
    all_found_flag=1;   //default true
    
    for(int i=0; i<dist_vec.size(); i++){
        //pos<- index of dist_vec[i] from tempList
        vector <int>::iterator it = tempList1.begin ();
        it = find (tempList1.begin (),tempList1.end (), dist_vec[i]);
        pos = (int)distance(tempList1.begin(), it);
        
        //if not found in list set flag to false and break
        //else delete tempList1[pos]
        if(dist_vec[i]!=tempList1[pos]){
            all_found_flag=0;
            break;
        }
        else{
            tempList1.erase(tempList1.begin()+pos);
        }
    }//end of for i=0~dist_vec.size()
    
    //if all_found_flag == true
    if(all_found_flag){
        //push m_y to list2 (x)
        x.push_back(m_y);
        //copy tempList1 to list1
        l.resize(tempList1.size());
        for(int k = 0; k < tempList1.size(); k++){
            l[k] = tempList1[k];
        }//end of copy for loop
        
        //recursive call with list1(l), list2(x)
        Place(l, x);
        //pos<- index of m_y in list2(x)
        vector <int>::iterator it = x.begin();
        it = find (x.begin (),x.end (), m_y);
        pos = (int)distance(x.begin(), it);
        //delete list2[pos]
        x.erase(x.begin()+pos);
        //restore list1 from dist_vec
        for(int i =0; i<dist_vec.size(); i++){
            l.push_back(dist_vec[i]);
        }
    }//end of if all_found_flag==true
    
    return;
}//end of Place(l,x)

int find_max(vector<int> l){
    int max=0;
    for(int i =0; i<l.size(); i++){
        if(l[i]> max)
            max = l[i];
    }
    
    return max;
}

int main(int argc, const char * argv[]) {
    string c;
    
    ifs.open("/Users/carlos/Desktop/csci191/L_prog2.txt");
    while(!(ifs.eof())){
        ifs>>c;
        list1.push_back(stoi(c));
    }
    ifs.close();
    
    cout<<"List1 = {";
    for(int i=0; i<list1.size(); i++){
        if(i == list1.size()-1)
            cout<<list1[i]<<"}"<<endl;
        else
            cout<<list1[i]<<", ";
    }
    m = find_max(list1);
    list1.erase(max_element(list1.begin(),list1.end()));
    list2.push_back(0);
    list2.push_back(m);
    
    //call place funtion
    Place(list1,list2);
    return 0;
}