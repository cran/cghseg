#include "usefullReadAndWrite.h"

using namespace std;

void readFile(char *nameFile, double * profil){
    ifstream inFile (nameFile);
    string line;
    int linenum = 0;
	
	int i=0;
    while (getline (inFile, line))
    {
        linenum++;
        istringstream linestream(line);
        string item;
        while (getline (linestream, item, ','))
        {
			profil[i] = atof(item.c_str());
            i++;
        }
    }
	inFile.close();
}

void writeFileD(char *nameFile, double profil){
	ofstream myfile;
	myfile.open(nameFile, ios::out | ios::app);
	myfile.precision(10);
	myfile << profil << std::endl;
	myfile.close();
}

void writeFileI(char *nameFile, int *profil, int nb, int jump){
	ofstream myfile;
	char *buffer = new char[40];
	int i =0;
	
	sprintf (buffer, "%s-%d ", nameFile, jump);
	//std::cout << jump << " , " << buffer << std::endl;
	myfile.open(buffer, ios::out | ios::app);
	i=0;
	while(i < jump){
		myfile << 0 << std::endl;
		i++;
	}
	while(i < nb){
		myfile << profil[i] << std::endl;
		i++;
	}
	myfile.close();
	delete(buffer);
}

void copyD(double *oldD, double *newD, int nb){
	int i=0;
	while(i < nb){
		newD[i] = oldD[i];
		i++;
	}
}

void traceback(char *nameFile, char *outFile, int nb, int Kmax){
	
    string line;
	int **res= new int* [Kmax];
	int i,j;
	i=0;
	while(i < Kmax){
		res[i] = new int [Kmax];
		i++;	
	}
	i=0;
	while(i < Kmax){
		j=0;
		while(j < Kmax){
		if(i == j){
			res[i][j] = nb;
		} else{
			res[i][j] = 0;
		}
		j++;
		}
		i++;	
	}
	int *profil= new int [nb];
	
	i=Kmax;
	ifstream inFile;
	char *buffer = new char[40];

	while(i > 1 ){
	buffer[0]='\0';
	sprintf (buffer, "%s-%d ", nameFile, i-1);
	//std::cout << i << " , " << buffer << " , "<< nameFile << std::endl;
	inFile.open(buffer);
	j=0;
		
    	string item;
		j=0;
    	while (j < nb){
				getline (inFile, item);
				profil[j] = atoi(item.c_str());
           	 j++;
        	}
		j=Kmax-1;
		while(j >= i-1){
			//std::cout << "j : " << j << ", i : " << i << ", res : "<< res[j][i-1] - 1 << ", profil : " << profil[(res[j][i-1] - 1)] << std::endl;
			res[j][i-2]	= profil[(res[j][i-1] - 1)];
			j--;
		}
		i--;
	inFile.close();
	} 
	delete(buffer);
	
	/* ecriture */
	ofstream myfile;
	myfile.open(outFile, ios::out | ios::app);
	i=0;
	while(i < Kmax){
		j=0;
		while(j < Kmax-1){
			myfile << res[i][j] << ";";
		j++;
		}
		myfile << res[i][j] << std::endl;
		i++;
	}
	
	delete(profil);
	i=0;
	while(i < Kmax){
		delete(res[i]);
		i++;	
	}
	delete(res);
	
}
