#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include "blas.h"

struct Desc{
	int n;
	int dimension;
	float* data;
	float* x;
	float* y;
};

struct Match {
	int id1,id2;
};

extern "C" {
void sgesvd_(char *jobu, char *jobvt, int *m, int *n, float *A, int *lda, 
	     float *S, float *U, int *ldu, float *VT, int *ldvt,
	     float *work, int *lwork, int *info);
}

void pca(Desc desc, int newDimension, float* basis) {
	char jobu = 'A', jobvt = 'N';
	int M = desc.dimension, N = desc.n;
	int l = M>N?N:M;
	int info,lwork=5*l;
	float* A = new float[M*N];
	memcpy(A,desc.data,M*N*sizeof(float));
	float* work = new float[lwork];
	float* S = new float[l];
	float* U = new float[M*M];
	sgesvd_(&jobu, &jobvt, &M, &N, A, &M, 
	     S, U, &M, NULL, &N, work, &lwork, &info);
	for (int i=0;i<M;i++)
		memcpy(basis + i*newDimension, U + i*M,newDimension * sizeof(float));
	delete[] work;
	delete[] S;
	delete[] U;
}

void changeBasis(Desc* desc, int newDimension, float* basis) {
	char transA = 'N';
	char transB = 'N';
	int M = newDimension, N = desc->n, K = desc->dimension;
	float alpha = 1,beta = 0;
	desc->dimension = newDimension;
	float* newData = new float[desc->n * newDimension];

	sgemm_(&transA,&transB,&M,&N,&K,
		&alpha,basis,&M,desc->data,&K,&beta,newData,&M);
	delete[] desc->data;
	desc->data = newData;
}

Desc parseKeyFile(char* filename) {
	Desc res;
	FILE* f = fopen(filename,"r");
	char buffer[128];
	fgets(buffer,128,f);
	char* c = buffer;
	res.n = strtol(c,&c,10);
	res.dimension = strtol(c,&c,10);
	res.data = new float[res.n * res.dimension];
	res.x = new float[res.n];
	res.y = new float[res.n];
	float* v = res.data;
	for (int i=0;i<res.n;i++) {
		float u;
		fscanf(f,"%f",res.y + i);
		fscanf(f,"%f",res.x + i);
		fscanf(f,"%f",&u);
		fscanf(f,"%f",&u);
		for (int j=0;j<res.dimension;j++)
			fscanf(f,"%f",v++);
	}
	fclose(f);
	return res;
}

std::vector<Match> findMatch(Desc source,Desc target,float ratio) {
	std::vector<Match> res;
	for (int i=0;i<source.n;i++) {
		int firstID=0,secondID=0;
		float firstScore,secondScore;
		for (int j=0;j<target.n;j++) {
			float score = 0;
			float *v1 = source.data+i*source.dimension;
			float *v2 = target.data+j*source.dimension;
			for (int k=0;k<source.dimension;k++,v1++,v2++)
				score += (*v1-*v2) * (*v1-*v2);
			if (j==0) 
				firstScore = secondScore = score;
			else if (score < firstScore) {
				firstScore = score;
				firstID = j;
			} else if (score < secondScore) {
				secondScore = score;
				secondID = j;
			}
		}
		if (firstScore < ratio * ratio * secondScore) {
			Match m = {i,firstID};
			res.push_back(m);
		}
	}
	return res;
} 

int main(int argc,char* argv[]) {
	if (argc < 5) {
		printf("%s ratio target id_start id_end\n",argv[0]);
		return 1;
	}

	float ratio = atof(argv[1]);
	int id_start = atoi(argv[3]);
	int id_end = atoi(argv[4]);

	Desc targetDesc = parseKeyFile(argv[2]);
	printf("Loaded %d (dim %d) keypoints from %s\n",targetDesc.n,targetDesc.dimension,argv[2]);
	
	int newDimension = 64;
	float* pca_basis = new float[newDimension * targetDesc.dimension];
	pca(targetDesc,newDimension,pca_basis);
	changeBasis(&targetDesc,newDimension,pca_basis);

	char buffer[128];
	struct timespec tic,toc;
	clock_gettime(CLOCK_MONOTONIC,&tic);
	for (int i=id_start;i<=id_end;i++) {
		sprintf(buffer,"%d.key",i);
		Desc sourceDesc = parseKeyFile(buffer);
		changeBasis(&sourceDesc,newDimension,pca_basis);
		std::vector<Match> match = findMatch(sourceDesc,targetDesc,ratio);
		sprintf(buffer,"%d.site.match",i);
		FILE* output = fopen(buffer,"w");
		for (size_t i=0;i<match.size();i++) {
			fprintf(output,"0 %f %f 1 %f %f\n",
				sourceDesc.x[match[i].id1],sourceDesc.y[match[i].id1],
				targetDesc.x[match[i].id2],targetDesc.y[match[i].id2]);
		}
		printf("Wrote %lu matches to %s\n",match.size(),buffer);
		fclose(output);

		delete[] sourceDesc.x;
		delete[] sourceDesc.y;
		delete[] sourceDesc.data;
	}
	clock_gettime(CLOCK_MONOTONIC,&toc);
	printf("Profiling: %fs for %d images\n",toc.tv_sec - tic.tv_sec + 0.000000001 * toc.tv_nsec - 0.000000001 * tic.tv_nsec,id_end-id_start+1);
	delete[] targetDesc.x;
	delete[] targetDesc.y;
	delete[] targetDesc.data;
	delete[] pca_basis;
}
