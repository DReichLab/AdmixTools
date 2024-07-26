/*
 * MapFile.h
 *
 *  Created on: Jun 21, 2021
 *      Author: csf
 */

#ifndef MAPFILE_H_
#define MAPFILE_H_

#include <stdbool.h>

typedef struct {	//struct containing pointer to genome data file
	char* snps;
	int rowSize;
	int colSize;
        int hdrSize ;
} genFile;

genFile* constructMap(char* fileName, long rowSize, long colSize, long hdrSize);

//helper function for checking correctness of array, not currently used in main
bool checkArray(long rowSize, long colSize, char* array);
int getgtypemap(genFile *data, int row, int col) ; 

void destroyMap(genFile* data);

#endif /* MAPFILE_H_ */
