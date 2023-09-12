#include<stdio.h>
#include<stdlib.h>
#include<string.h>

extern void ffoo( char* );

int main( int argc, char* argv[] )
{
	char c_string[64];
	memset(c_string,0,sizeof(c_string));
	strcpy(c_string,argv[1]);
	ffoo(c_string);

	strcpy(c_string,argv[2]);
	ffoo(c_string);

	return 0;
}
