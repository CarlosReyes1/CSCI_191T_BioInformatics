//
//  Carlos Reyes
//  Student ID: 108939200
//  Date: 4/3/2016
//  main.cpp
//  Csci191T_Prog4
//
/*
 ****************************************************************************
 *                                  Global Documentation
 * This program is an implementation of the Smith-Waterman local alignment
 * algorithm. The program uses linear gap scroing and includes a trace-back
 * operation in order to display the best local alignment.The program will take
 * in two files, the first being the query seqeunce file and the second being
 * the database sequence file. Upon opening the query sequence file, the
 * program will read the first line from the query file which is the label.
 * The label is not needed for this program. After reading the label the
 * program uses a generic method of reading the sequence in case there are
 * multiple lines for the sequence. The program will loop through the file
 * while reading every line and storing/adding the sequence string to a string
 * variable. Once the file has read the entire query sequence the file will
 * close and then open the database file. Once the database file is opened the
 * program will store and display the label, then the program will read the
 * sequence, display its length, and call upon the smith-waterman function.
 * The Smith-Waterman function creates a struct matrix which will hold the score
 * and direction for the best local alignment. The linear gap scoring will be
 * determined by the diagonal match case of +2 and the mismatch case of -1,
 * insertion or deletion case of -1. Once the scoring matrix is complete
 * then the program will trace-back operation in order to create the the best
 * local alignment strings fromt the query sequence and database sequence.
 * The sequences will be displayed with the score, starting and end positions,
 * and indel or match/mismatch.
 *****************************************************************************
 *                          How to run the program
 * This program was written in C++. In order to run the program through the
 * terminal, first you must compile the program using this command line
 * g++ main.cpp -o main.out
 * Once the file has been compiled, it can then be executed using the 
 * command line
 * ./main.out Prog4-query.fa Prog4-database.fa > Prog4_sol.txt
 *            ^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^
 *                  ^               ^                   ^
 *                  |               |                   |
 *          (NameOfQueryfile) (NameOfDBfile)     (NameOfOutputfile)
 * The Query file must the first parameter and the Database file must be 
 * the second paramater in order for this program to function properly.
 * the third parameter that is in front of '>' is the output file which
 * will the solution will be written to.
 *
 ****************************************************************************
 *                      Program Starts Here
*/

#include <iostream>
#include <string>
#include <fstream>

using namespace std;

ifstream q_seq_file; //sequence input file stream
ifstream db_seq_file;   //db sequence file stream

// Struct the will hold the direction and score for the scoring matrix
struct v_tTable{
    int dir, score;
};

// Smith-Waterman function which will take in two paramaters and does not
// return anything.
// first paramater is the query seqence
// second paramter is database sequence
// will geneterate a struct matrix where the row is the query sequence size+1
// the column is the db sequence size+1
// first row and column will be zeros
// Then the score matrix will be computed using the gap, match, and mistmatch
// scoring. The max value will be stored along with its direction
// Once entire matrix is computed the trace-operation will be used and
// the best local alignment strings will be stored and displayed. Once
// completed the function will return back to main
void s_w(string q_sq, string db_sq){
    v_tTable v_t [q_sq.length()+1][db_sq.length()+1];
    int best_sw_score = 0, best_sw_i=0, best_sw_j=0;
    int end_i, end_j;
    string q_str = "";
    string db_str = "";
    string center_str ="";
    
    // Zero out first row
    for(int i = 0; i < q_sq.length()+1; i++){
        v_t[i][0].score =0;
        v_t[i][0].dir =0;
    }
    // Zero out first column
    for(int i = 0; i< db_sq.length()+1; i++){
        v_t[0][i].score =0;
        v_t[0][i].dir =0;
    }
    
    // Compute Scoring Matrix using Smith-Waterman algorithm
    for(int i= 1; i<q_sq.length()+1; i++){
        for(int j = 1; j<db_sq.length()+1; j++){
            int diagonal=0, upper=0, left=0;
            int direction=0;
            int max_val=0;
            if(q_sq[i-1] == db_sq[j-1]){
                diagonal=v_t[i-1][j-1].score+2;
            }
            else if(q_sq[i-1] != db_sq[j-1]){
                diagonal=v_t[i-1][j-1].score-1;
            }
            
            left = v_t[i][j-1].score-1;
            upper = v_t[i-1][j].score-1;
            
            if(diagonal<=0 && left <=0 && upper<=0){
                max_val=0;
                direction =0;
            }
            else if(diagonal >= left && diagonal >= upper){
                max_val = diagonal;
                direction =1;
            }
            else if(upper >= left && upper > diagonal){
                max_val = upper;
                direction = 2;
            }
            else if(left > diagonal && left > upper){
                max_val = left;
                direction =3;
            }
            v_t[i][j].score = max_val;
            v_t[i][j].dir = direction;
            if(v_t[i][j].score > best_sw_score){
                best_sw_score = v_t[i][j].score;
                best_sw_i =i;
                best_sw_j = j;
            }
        }
    }
    // store ending indexes for display purposes
    end_i = best_sw_i;
    end_j = best_sw_j;
    
    // Trace-back operation
    while(v_t[best_sw_i][best_sw_j].score !=0){
        if(v_t[best_sw_i][best_sw_j].dir==1){
            q_str = q_sq[best_sw_i-1] + q_str ;
            db_str = db_sq[best_sw_j-1] + db_str;
            if(db_sq[best_sw_j-1]==q_sq[best_sw_i-1])
               center_str = '|' + center_str ;
            else
               center_str =' ' + center_str;
            best_sw_i--;
            best_sw_j--;
        }
        else if(v_t[best_sw_i][best_sw_j].dir==2){
            q_str = q_sq[best_sw_i-1] + q_str ;
            db_str = '-' + db_str;
            center_str =' ' + center_str;
            best_sw_i--;
        }
        else if(v_t[best_sw_i][best_sw_j].dir==3){
            q_str = '-' + q_str ;
            db_str = db_sq[best_sw_j-1] + db_str;
            center_str = ' ' + center_str;
            best_sw_j--;
        }
    }
    
    //Displaying vital information
    cout<<"Optimum Smith-Waterman score = "<<best_sw_score<<endl;
    cout<<"Query : alignment start index = "<<best_sw_i;
    cout<<", end index = "<<end_i<<endl;
    cout<<"DB seq: alignment start index = "<<best_sw_j;
    cout<<", end index = "<<end_j<<endl;
    cout<<q_str<<endl;
    cout<<center_str<<endl;
    cout<<db_str<<endl<<endl;

    return;
}

// Opens two files, first the query file then the database file
// Displays the db label then calls upon the Smith-Waterman function
// Will do this for every sequence in the DB file until we hit
// the end of the file.
int main(int argc, const char * argv[]) {
    string query_label, query_seq;
    string db_label, db_seq;
    string str1;
    //opens query file
    q_seq_file.open(argv[1]);
    
    //gets the label from query file
    getline(q_seq_file, query_label);
    //reads all lines after and stores the sequence into one string
    while(!q_seq_file.eof()){
        getline(q_seq_file, str1);
        query_seq += str1;
    }
    //closes file
    q_seq_file.close();
    //opends the database file
    db_seq_file.open(argv[2]);
    //loops till end of file is reached
    while(!db_seq_file.eof()){
        //reads and stores label from db file
        getline(db_seq_file, db_label);
        //reads and stores sequence from db file
        getline(db_seq_file, db_seq);
        //display db label
        cout<<db_label<<"       len = "<<db_seq.length()<<endl;
        //calls smith-waterman function and passes the query sequence
        //and the db sequence
        s_w("TTAAC", "GTTAC");
    }
    //closes db file
    db_seq_file.close();
    //program completed
    return 0;
}