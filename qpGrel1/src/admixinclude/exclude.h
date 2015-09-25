#ifdef __cplusplus
extern "C" {
#endif
#ifndef _EXCLUDE_
#define _EXCLUDE_

#include <admutils.h>
#include <mcio.h>
#include <stdio.h>

/*  file name parameter : xregionname */
/*  HW filter parameter : nhwfilter (-1 means no-filter) */
/*  maximum number of regions : 1000
    closed intervals in physical position include endpoints */ 
/*  read file and set ignore flag for SNPs */
void excluderegions(char *xregionname, SNP **snps, int nsnps, char *deletesnpoutname);
void hwfilter(SNP **snps, int nsnps, int nindiv, double nhwfilter, char *deletesnpoutname);

#endif
#ifdef __cplusplus
}
#endif
