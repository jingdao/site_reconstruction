#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <vector>
#define RANSAC_ITERS 1000
#define RANSAC_THRESHOLD 10

struct Point {
	float x,y,z;
};

struct Param {
	float theta,k,Tx,Ty;
};

int main(int argc, char* argv[]) {
	if (argc<3) {
		printf("%s target_point.txt camera_location_2d.txt [1.site.match ..]\n",argv[0]);
		return 1;
	}

	srand(time(NULL));
	char buffer[128];
	std::vector<float> x_target,y_target;
	std::vector<int> match_target;
	std::vector<int> numObjects;
	FILE* target_point = fopen(argv[1],"r");
	if (!target_point) {
		printf("Cannot open %s\n",argv[1]);
		return 1;
	}
	while (fgets(buffer,128,target_point)) {
		float x,y;
		int numMatch;
		char* c = buffer;
		int n = strtol(c,&c,10);
		numObjects.push_back(n);
		for (int i=0;i<n;i++) {
			x = strtod(c,&c);
			y = strtod(c,&c);
			numMatch = strtol(c,&c,10);
			x_target.push_back(x);
			y_target.push_back(y);
			match_target.push_back(numMatch);
		}
	}
	fclose(target_point);

	FILE* camera_location = fopen(argv[5],"w");
	int totalObjects=0;
	for (size_t k=0;k<numObjects.size();k++) {
		//get keypoint matches from file
		std::vector<float> x_src,y_src,x_dst,y_dst;
		sprintf(buffer,"%lu.site.match",k+1);
		FILE* siteMatch = fopen(buffer,"r");
		while (fgets(buffer,128,siteMatch)) {
			int id1,id2;
			float x1,x2,y1,y2;
			if (sscanf(buffer,"%d %f %f %d %f %f",&id1,&x1,&y1,&id2,&x2,&y2)==6) {
				x_src.push_back(x1);
				y_src.push_back(y1);
				x_dst.push_back(x2);
				y_dst.push_back(y2);
			}
		}
		fclose(siteMatch);
		//use RANSAC to obtain best parameters
		int maxInliers=0;
		Param bestParam;
		for (int i=0;i<RANSAC_ITERS;i++) {
			Param p;
			int id1 = rand() % x_src.size();
			int id2 = rand() % x_src.size();
			float x1 = x_src[id1];
			float y1 = y_src[id1];
			float x1p = x_dst[id1];
			float y1p = y_dst[id1];
			float x2 = x_target[id2];
			float y2 = y_target[id2];
			float x2p = x_dst[id1];
			float y2p = y_dst[id1];
			float r = (y2 - y1) / (x2 - x1);
			float q = (y2p - y1p) / (x2p - x1p);
			p.theta = atan((r - q) / (r*q + 1));
			p.k = (x2p - x1p) / (cos(p.theta)*(x2-x1) + sin(p.theta)*(y2-y1));
			p.Tx = x1p - k*(cos(p.theta)*x1 + sin(p.theta)*y1);
			p.Ty = y1p - k*(-sin(p.theta)*x1 + cos(p.theta)*y1);
			int numInliers=0;
			for (int j=0;j<x_src.size();j++) {
				float dx = p.k*(cos(p.theta)*x_src[j] + sin(p.theta)*y_src[j]) + p.Tx - x_dst[j];
				float dy = p.k*(-sin(p.theta)*x_src[j] + cos(p.theta)*y_src[j]) + p.Ty - y_dst[j];
				if (dx*dx + dy*dy < RANSAC_THRESHOLD)
					numInliers++;
			}
			if (numInliers > maxInliers) {
				maxInliers = numInliers;
				bestParam = p;
			}
		}
		//write transformed coordinates to file
		fprintf(camera_location,"%d",numObjects[k]);
		for(int i=0;i<numObjects[k];i++) {
			float x = x_target[totalObjects + i];
			float y = y_target[totalObjects + i];
			float xp = bestParam.k*(cos(bestParam.theta)*x + sin(bestParam.theta)*y) + bestParam.Tx;
			float yp = bestParam.k*(-sin(bestParam.theta)*x + cos(bestParam.theta)*y) + bestParam.Ty;
			fprintf(camera_location," %f %f %d",xp,yp,match_target[totalObjects+i]);
		}
		fprintf(camera_location,"\n");
		totalObjects += numObjects[k];
	}
	fclose(camera_location);
	return 1;

}
