#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int ARRAY_WIDTH = 212;
double array[212][3];

int main (int argc, char *argv[]) {

    double *array[ARRAY_WIDTH];
    int idx = 0;
    int j = 0;
    char *buffer = NULL;
    size_t len = 0;
    ssize_t read;
    char *ptr = NULL;

    FILE *fp;
    fp = fopen ("peptide-T-APM.csv", "r");      //open file , read only
    if (!fp) {
        fprintf (stderr, "failed to open file for reading\n");
        return 1;
    }

    while ((read = getline (&buffer, &len, fp)) != -1) {

        array[idx] = malloc (sizeof (array));

        for (j = 0, ptr = buffer; j < 3; j++, ptr++){
            array [idx][j] = (double)strtod(ptr, &ptr);
		}
        idx++;
    }

    fclose (fp);
    
    printf("m/z      retentiontime        M+H\n");

	for (int i = 0; i < ARRAY_WIDTH; i++){
		for (int j = i + 1; j < ARRAY_WIDTH; j++){
			double mzij = array[j][0] - array[i][0];
			double rtij = array[j][1] - array[i][1];
			double mhij = array[j][2] - array[i][2];
			if (fabs(mhij)/array[j][2] <= 0.000001 && fabs(rtij) >= 1.0 && fabs(mzij)/array[j][0] <= 0.000001){
				printf("%f, %f, %f, | %f, %f, %f\n", array[i][0], array[i][1], array[i][2], array[j][0], array[j][1], array[j][2]);
			}
		}
	}
    return 0;
}
	
		