#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <time.h>
#include <omp.h>

SEXP recode(SEXP g1, SEXP g2, SEXP unique_list, SEXP ncpu)
{
	// OMP STUFF
	omp_set_dynamic(0);
	omp_set_num_threads(asInteger(ncpu));
	// nrow IS THE NUMBER OF SAMPLE
	// g1 AND g2 SHOULD HAVE SAME DIM
	int nrow=INTEGER(GET_DIM(g1))[0];
	int loci=INTEGER(GET_DIM(g1))[1];
	SEXP output=PROTECT(allocMatrix(RAWSXP, nrow, loci));
	Rbyte *h_g1=RAW(g1); Rbyte *h_g2=RAW(g2); 
	Rbyte *h_output=RAW(output);
	#pragma omp parallel for
	for (int i=0; i<loci; i++)
	{
		int offset=i*nrow;
		Rbyte larger_allele;
		int num_allele=LENGTH(VECTOR_ELT(unique_list, i));
		// IF FIXED, len_unique IS 1, THEN ALL GENOTYPE IS CODED AS 0
		if (num_allele==1)
		{
			# pragma omp simd
			for (int j=0; j<nrow; j++)
			{
				h_output[offset+j]=0;
			}
		}
		// ELSE, SMALLER HOMOZYGOTE=0, HETEROZYGOTE=1, LARGER HOMOZYGOTE=2
		else
		{
			larger_allele=RAW(VECTOR_ELT(unique_list, i))[1];
			for (int j=0; j<nrow; j++)
			{
				h_output[offset+j]=(h_g1[offset+j]==larger_allele)+(h_g2[offset+j]==larger_allele);
			}
		}
	}
	UNPROTECT(1);
	return output;
}
