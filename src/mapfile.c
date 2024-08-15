/*
 * MapFile.c
 *
 *  Created on: Jun 21, 2021
 *      Author: Cokie; modified by Nick
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>
#include <stdint.h>
#include <errno.h>

#include "mapfile.h"
#include "admutils.h"  

// Memory mapping by Cokie Parke 

//function returns value of single 2-bit snp, given
//row and col location and genFile object containing data matrix
//return int value
int getgtypemap(genFile *data, int snpRow, int snpCol) 
{
	char* snps = data->snps;
        static long ncall = 0 ; 
        long rowsize, offset ; 
        char *buff ;  

        ++ncall ;
        rowsize = data -> rowSize ;  // rowsize is size of a row NOT number of rows 
        offset = data -> hdrSize ;  

        buff = snps + offset + snpRow*rowsize ;  // start of row  
        return rbuff((unsigned char *) buff, snpCol) ;   

} 

//memory maps file with given filename, returns
//genFile struct containing pointer representing file
genFile* constructMap(char* filename, long numrows, long rowSize, long hdrSize) {
	genFile* geno;
        long siz ; 

	//open file to be mapped
	int mmfd = open(filename, O_RDWR);

	if (mmfd <  0) fatalx("(constructMap) failure to open %s\n", filename) ;

	//get file length
	struct stat status;
	if (fstat(mmfd, &status) != 0) fatalx("fstat fails\n") ;
	long fileLength = status.st_size;

	siz = numrows * rowSize + hdrSize ; 
	if (siz != fileLength) fatalx("(constructMqp) file size mismatch\n") ;

	//memory map file
	void *mmfaddr = mmap(NULL, fileLength, PROT_READ, MAP_PRIVATE, mmfd, 0);
	if (mmfaddr == (void*)-1) {
		printf("error number: %d/n", errno);
	}

	//close file
	close(mmfd);

	//build genFile struct
	geno = malloc(sizeof(genFile));
	geno->colSize = numrows ;
	geno->rowSize = rowSize;
        geno -> hdrSize = hdrSize ;
	geno->snps = (char*) mmfaddr;

        printf("mmfaddr: %p  fileLength: %ld\n", geno -> snps, fileLength) ;;
        fflush(stdout) ;

	return geno;
}

//unmaps data file contained in given genFile struct
//frees dynamically allocated struct container
void destroyMap(genFile* data) {
	void* mmfaddr = (void*) data->snps;
	int fileLength = data->colSize * data->rowSize + data -> hdrSize;
	munmap(mmfaddr, fileLength);
	free(data);
}
