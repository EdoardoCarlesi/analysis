/* General Main for GADGET-2 Analysis Routines */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

int main(int argc, char *argv[]){
	
	fprintf(stderr, "\n");
	printf(" =========================================== \n");
	printf(" = GADGET - 2 ANALYSIS ROUTINES COLLECTION = \n" );
	printf(" =========================================== \n");
	
	int value=1;

	do {
	fprintf(stderr, "\n\nChoose the type of your analysis: \n\n");
	fprintf(stderr, "(1) Fit haloes' densities to a NFW profile \n"); 
	fprintf(stderr, "(2) Generate power spectra \n"); 
	fprintf(stderr, "(3) Compute the growth factor \n"); 
	fprintf(stderr, "(4) Compute halo mass function (theoretical and numerical) \n"); 
	fprintf(stderr, "(5) Get halo statistical properties \n"); 
	fprintf(stderr, "(6) Compute number densities\n"); 
	fprintf(stderr, "(7) Halo and subhalo evolution\n"); 
	fprintf(stderr, "(8) Compare cross correlated haloes \n"); 
	fprintf(stderr, "(9) Get subhalo properties\n"); 
	fprintf(stderr, "(10) Run a test\n"); 
	fprintf(stderr, "\nAny other key to exit.\n");

	char input[10];
	char str[250];
	scanf("%s", input);
	int check = atoi(input);

{
	if(check==1 || check==2 ||  check==3 || check==4 || check==5 || check==6 
			|| check==7 || check==8 || check==9  || check==10) {
	value=1;
	sprintf(str, "./scripts/execute.sh %d", check);
	system(str);
	} else {
	value =0;
	} }
	} while (value == 1);
	
	fprintf(stderr, "\n Execution terminated. \n\n");
	return 0;
}
