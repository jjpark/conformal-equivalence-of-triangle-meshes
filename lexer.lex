/* Comments go here */
%{
#include <iostream>
#include <list>
#include "data.h"
/* Include tokens and variables from the parser (NUM, OPEN, CLOSE, yylval) */
#include "parser.tab.hpp"
using namespace std;
%}
%option noyywrap


DIGIT [0-9]
INT \-?[0-9]+
FLOAT \-?[0-9]*\.?[0-9]+|\-?[0-9]+\.?[0-9]*|\-?[0-9]\.?[0-9]*e\-[0-9]*

%%
v  								return V;
f                               return F;
{FLOAT}							{ yylval.fval = atof(yytext); return NUMBER; }
[ \t\r\n\f\v]+   				/* remove whitespace */
.           { cerr << "Unexpected character: " << yytext << endl;
			  yyterminate();
            }
%%

/* Extra code goes here */
