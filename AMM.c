#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int data_WIDTH = 500;    //data # of lines
double data[500][0];
int database_WIDTH = 6332;    //database # of lines
double database[6332][0];

int main (int argc, char *argv[]) {

    double *data[data_WIDTH];
    int idx1 = 0;
    int j1 = 0;
    char *buffer1 = NULL;
    size_t len1 = 0;
    ssize_t read1;
    char *ptr1 = NULL;

//Read the data file
    FILE *fp1;
    fp1 = fopen ("110217-rock brain10252017-3-Top500mz.csv", "r");      //open MALDI exported data file , read only
    if (!fp1) {
        fprintf (stderr, "failed to open file for reading\n");
        return 1;
    }

    while ((read1 = getline (&buffer1, &len1, fp1)) != -1) {

        data[idx1] = malloc (sizeof (data));

        for (j1 = 0, ptr1 = buffer1; j1 < 1; j1++, ptr1++){
            data[idx1][j1] = (double)strtod(ptr1, &ptr1);
		    }
        idx1++;
    }

    fclose (fp1);

  //Read the database file
    double *database[database_WIDTH];
    int idx2 = 0;
    int j2 = 0;
    char *buffer2 = NULL;
    size_t len2 = 0;
    ssize_t read2;
    char *ptr2 = NULL;

    FILE *fp2;
    fp2 = fopen ("lipid database_v1.csv", "r");      //open database file , read only
    if (!fp2) {
        fprintf (stderr, "failed to open file for reading\n");
        return 1;
    }

    while ((read2 = getline (&buffer2, &len2, fp2)) != -1) {

        database[idx2] = malloc (sizeof (database));

        for (j2 = 0, ptr2 = buffer2; j2 < 3; j2++, ptr2++){
        
            database[idx2][j2] = (double)strtod(ptr2, &ptr2);
          
		}  
		    	
        idx2++;
    }

    fclose (fp2);

	printf("M+H in data, M+H in database, mass tolerance (ppm), Acyl C# (m), Acyl C=C # (n)\n");
	for (int i = 0; i < data_WIDTH; i++){
		for (int j = 0; j < database_WIDTH; j++){
			double mzij = data[i][0] - database[j][0];
			if (fabs(mzij)/data[i][0] <= 0.00001){
				printf("%f, %f, %f, %d, %d\n", data[i][0], database[j][0], fabs(mzij)/database[j][0]*1000000.0, (int)database[j][1], (int)database[j][2]);
			}
		}
	}
    return 0;
}
