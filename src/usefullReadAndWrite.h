#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


using namespace std;

void readFile(char *nameFile, double * profil);
void writeFileD(char *nameFile, double profil);
void writeFileI(char *nameFile, int *profil, int nb, int jump);
void copyD(double *oldD, double *newD, int nb);
void traceback(char *nameFile, char *outFile, int nb, int Kmax);

