//===========================================================================
//                                                                           
// File: removecomments.C                                                    
//                                                                           
// Created: Thu Feb 10 13:45:40 2005                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: removecomments.C,v 1.1 2005-12-13 12:08:25 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include <iostream>

using namespace std;

inline bool start_token_found(char c1, char c2, bool& c_start, bool& cpp_start)
{
    c_start = (c1 == '/' && c2 == '*');
    cpp_start = (c1 == '/' && c2 == '/');
    return c_start || cpp_start;
}

inline bool all_whitespace(const string& s) 
{
    for (int i = 0; i < s.size(); ++i) {
	char c = s[i];
	if (c != ' ' && c != '\t' && c != '\n') return false;
    }
    return true;
}

int main()
{
    bool in_c_comment = false;
    bool in_cpp_comment = false;

    while (!cin.eof()) {
	string cur_line;
	bool line_contained_end_token = false;
	char c = ' ';
	while (c != '\n' && !cin.eof()) {
	    char cprev = c;
	    cin.read(&c, 1);

	    // looking for end token of comment block
	    if (in_c_comment) {
		if (cprev == '*' && c == '/') {
		    in_c_comment = false;
		    line_contained_end_token = true;
		    cprev = ' '; // set cprev to avoid recognition of "false" token 
		    cin.read(&c, 1);
		}
	    } else if (in_cpp_comment) {
		if (c == '\n') {
		    in_cpp_comment = false;
		    line_contained_end_token = true;
		}
	    }

	    bool not_in_comment = ! (in_c_comment || in_cpp_comment);
	    
	    if (not_in_comment) {
		if (start_token_found(cprev, c, in_c_comment, in_cpp_comment)) {
		    cur_line.resize(cur_line.size() - 1); // remove last character
		} else {
		    // business as usual
		    cur_line += c;
		}
	    }
	}
	if (!line_contained_end_token || !all_whitespace(cur_line)){
	    cout << cur_line;
	}
    }
    cout << '\n';
};


