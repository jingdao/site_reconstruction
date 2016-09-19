#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <vector>
#include <string.h>
#include <math.h>

typedef struct {
	float x,y;
} Coordinate;

int sortByXY(const void* v1, const void* v2) {
	Coordinate* c1 = (Coordinate*)v1;
	Coordinate* c2 = (Coordinate*)v2;
	if (c1->x == c2->x)
		return c1->y - c2->y;
	else
		return c1->x - c2->x;
}

int cross(Coordinate O, Coordinate A, Coordinate B) {
	return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}

std::vector<Coordinate> convexHull(std::vector<Coordinate> P) {
	int n = P.size(), k = 0;
	std::vector<Coordinate> H(2*n);
	if (n==0)
		return H;

	// Sort points lexicographically
	qsort(P.data(),n,sizeof(Coordinate),sortByXY);

	// Build lower hull
	for (int i = 0; i < n; ++i) {
		while (k >= 2 && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
		H[k++] = P[i];
	}

	// Build upper hull
	for (int i = n-2, t = k+1; i >= 0; i--) {
		while (k >= t && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
		H[k++] = P[i];
	}

	H.resize(k-1);
	return H;
}

float polygonArea(std::vector<Coordinate> P) {
	float sum=0;
	for (size_t i=0;i<P.size()-1;i++) {
		sum += P[i].x * P[i+1].y;
		sum -= P[i].y * P[i+1].x;
	}
	sum += P[P.size()-1].x * P[0].y;
	sum -= P[P.size()-1].y * P[0].x;
	return fabs(sum) / 2;
}

std::vector<Coordinate> polygonCombine(std::vector<Coordinate> P1, std::vector<Coordinate> P2) {
	std::vector<Coordinate> res;
	for (size_t i=0;i<P1.size();i++)
		res.push_back(P1[i]);
	for (size_t i=0;i<P2.size();i++)
		res.push_back(P2[i]);
	return res;
}

int main(int argc,char* argv[]) {
	if (argc < 3) {
		printf("%s target_point.txt label_point.txt\n",argv[0]);
		return 1;
	}
	
	FILE* target = fopen(argv[1],"r");
	if (!target) {
		printf("Cannot open %s\n",argv[1]);
		return 1;
	}
	FILE* label = fopen(argv[2],"r");
	if (!label) {
		printf("Cannot open %s\n",argv[2]);
		return 1;
	}
	char buffer[512];
	std::vector< std::vector<Coordinate> > target_box;
	std::vector< std::vector<Coordinate> > label_box;
	int numObjects;
	while (fgets(buffer,512,target)) {
		float x,y;
		int numMatch;
		char* c = buffer;
		int n = strtol(c,&c,10);
		numObjects = n/4;
		std::vector<Coordinate> box;
		for (int i=0;i<n;i++) {
			x = strtod(c,&c);
			y = strtod(c,&c);
			numMatch = strtol(c,&c,10);
			Coordinate p = {x,y};
			box.push_back(p);
			if (i%4==3) {
				target_box.push_back(box);
				box.clear();
			}
		}
	}
	while (fgets(buffer,512,label)) {
		float x,y;
		int numMatch;
		char* c = buffer;
		int n = strtol(c,&c,10);
		std::vector<Coordinate> box;
		if (n==0) {
			for (int i=0;i<numObjects;i++)
				label_box.push_back(box);
			continue;
		}
		for (int i=0;i<n;i++) {
			x = strtod(c,&c);
			y = strtod(c,&c);
			numMatch = strtol(c,&c,10);
			Coordinate p = {x,y};
			box.push_back(p);
			if (i%4==3) {
				label_box.push_back(box);
				box.clear();
			}
		}
	}

	float acc_avg=0,acc_stddev=0;
	float acc_min,acc_max;
	float prec_avg=0,prec_stddev=0;
	float prec_min,prec_max;
	int count=0;
	for (size_t i=0;i<target_box.size();i++) {
		if (label_box[i].size()==0)
			continue;
		float A1 = polygonArea(convexHull(target_box[i]));
		float A2 = polygonArea(convexHull(label_box[i]));
		float A3 = polygonArea(convexHull(polygonCombine(target_box[i],label_box[i])));
		float acc = A3 > A1 + A2 ? 0 : (A1 + A2 - A3) / A2;
		float prec = A3 > A1 + A2 ? 0 : (A1 + A2 - A3) / A1;
		if (acc > 1) acc = 1;
		if (prec > 1) prec = 1;
		acc_avg += acc;
		acc_stddev += acc*acc;
		if (i==0 || acc < acc_min) acc_min = acc;
		if (i==0 || acc > acc_max) acc_max = acc;
		prec_avg += prec;
		prec_stddev += prec*prec;
		if (i==0 || prec < prec_min) prec_min = prec;
		if (i==0 || prec > prec_max) prec_max = prec;
		printf("Box %lu: %.2f %.2f %8.2f %.4f %.4f\n",i,A1,A2,A3,acc,prec);
		count++;
	}
	acc_avg /= count;
	acc_stddev = sqrt(acc_stddev/count - acc_avg*acc_avg);
	prec_avg /= count;
	prec_stddev = sqrt(prec_stddev/count - prec_avg*prec_avg);
	printf("Accuracy Min %.4f Max %.4f Average %.4f +- %.4f Overlap\n",acc_min,acc_max,acc_avg,acc_stddev);
	printf("Precision Min %.4f Max %.4f Average %.4f +- %.4f Overlap\n",prec_min,prec_max,prec_avg,prec_stddev);

	fclose(target);
	fclose(label);

}
