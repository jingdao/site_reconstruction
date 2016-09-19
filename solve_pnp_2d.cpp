#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <vector>
#define RANSAC_ITERS 10000
#define RANSAC_THRESHOLD 50

typedef struct {
	unsigned char r,g,b;
} Color;

struct Point {
	float x,y,z;
};

struct Param {
	float theta,k,Tx,Ty;
};

typedef struct {
	int width,height;
	unsigned char* data;
} Image;

void drawLine(Image* surf, float x1, float y1, float x2, float y2, Color color ) {
	if (x1 < 1) x1 = 1;
	if (x1 > surf->width - 1) x1 = surf->width - 1;
	if (x2 < 1) x2 = 1;
	if (x2 > surf->width - 1) x2 = surf->width - 1;
	if (y1 < 1) y1 = 1;
	if (y1 > surf->height - 1) y1 = surf->height - 1;
	if (y2 < 1) y2 = 1;
	if (y2 > surf->height - 1) y2 = surf->height - 1;
	// Bresenham's line algorithm
	const bool steep = (fabs(y2 - y1) > fabs(x2 - x1));
	if(steep) {
		std::swap(x1, y1);
		std::swap(x2, y2);
	}
	if(x1 > x2) {
		std::swap(x1, x2);
		std::swap(y1, y2);
	}

	const float dx = x2 - x1;
	const float dy = fabs(y2 - y1);
	float error = dx / 2.0f;
	const int ystep = (y1 < y2) ? 1 : -1;
	int y = (int)y1;
	const int maxX = (int)x2;
	unsigned char* dst;

	for(int x=(int)x1; x<maxX; x++) {
		if(steep)
			dst = surf->data + (x * surf->width + y)*3;
		else
			dst = surf->data + (y * surf->width + x)*3;
		dst[0] = color.r;
		dst[1] = color.g;
		dst[2] = color.b;
		error -= dy;
		if(error < 0) {
			y += ystep;
			error += dx;
		}
	}
}

void drawBox(Image* img, float* x,float* y,Color color) {
	drawLine(img,x[0],y[0],x[1],y[1],color);
	drawLine(img,x[0],y[0],x[2],y[2],color);
	drawLine(img,x[1],y[1],x[3],y[3],color);
	drawLine(img,x[2],y[2],x[3],y[3],color);
}

Image imclone(Image image) {
	Image im = {image.width,image.height,NULL};
	im.data = new unsigned char[im.width*im.height*3];
	memcpy(im.data,image.data,im.width*im.height*3);
	return im;
}

bool loadImage(char* name,Image *image,bool color) {
	char buffer[128];
	FILE* pgm = fopen(name,"r");
	if (!pgm) {
		printf("%s not found\n",name);
		return false;
	}
	fgets(buffer,128,pgm); //P5 or P6
	do {
		fgets(buffer,128,pgm);
	} while (buffer[0]=='#'); //remove comments
	char *c = buffer;
	image->width = strtol(c,&c,10);
	image->height = strtol(c,&c,10);
	fgets(buffer,128,pgm); //255
	if (image->data)
		delete[] image->data;
	image->data = new unsigned char[image->width*image->height*3];
	if (color) {
		fread(image->data,1,image->width*image->height*3,pgm);
//		for (int i=image->width*image->height;i>=0;i--) {
//			unsigned char tmp = image->data[i*3+2];
//			image->data[i*3+2] = image->data[i*3];
//			image->data[i*3] = tmp;
//		}
	} else {
		fread(image->data,1,image->width*image->height,pgm);
		for (int i=image->width*image->height;i>=0;i--) {
			image->data[i*3+2] = image->data[i*3+1] = image->data[i*3] = image->data[i];
		}
	}
	return true;
}

void writeImage(char* name,Image *image,bool color) {
	char buffer[128];
	FILE* f = fopen(name,"w");
	if (color) {
		fprintf(f,"P6\n%d %d\n255\n",image->width,image->height);
		fwrite(image->data,1,image->width*image->height*3,f);
	} else {
		fprintf(f,"P5\n%d %d\n255\n",image->width,image->height);
		fwrite(image->data,1,image->width*image->height,f);
	}
	fclose(f);
}

int main(int argc, char* argv[]) {
	if (argc<4) {
		printf("%s target_point.txt camera_location_2d.txt site.ppm [1.site.match ..]\n",argv[0]);
		return 1;
	}

	srand(time(NULL));
	char buffer[1024];
	std::vector<float> x_target,y_target;
	std::vector<int> match_target;
	std::vector<int> numObjects;
	FILE* target_point = fopen(argv[1],"r");
	if (!target_point) {
		printf("Cannot open %s\n",argv[1]);
		return 1;
	}
	while (fgets(buffer,1024,target_point)) {
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

	Image site = {0,0,NULL};
	if (!loadImage(argv[3],&site,true)) {
		return 1;
	}
	Color yellow = {255,255,0};
	Color red = {255,0,0};
	Color blue = {0,0,255};

	FILE* camera_location = fopen(argv[2],"w");
	int totalObjects=0;
	float score=0;
	for (size_t k=0;k<numObjects.size();k++) {
		//get keypoint matches from file
		std::vector<float> x_src,y_src,x_dst,y_dst;
		sprintf(buffer,"%lu.site.match",k+1);
		FILE* siteMatch = fopen(buffer,"r");
		if (!siteMatch) {
			fprintf(camera_location,"0\n");
			totalObjects += numObjects[k];
			continue;
		}
		while (fgets(buffer,1024,siteMatch)) {
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
			float x2 = x_src[id2];
			float y2 = y_src[id2];
			float x2p = x_dst[id2];
			float y2p = y_dst[id2];
			if (x1==x2 || x2p==x1p)
				continue;
			float r = (y2 - y1) / (x2 - x1);
			float q = (y2p - y1p) / (x2p - x1p);
			p.theta = atan((r - q) / (r*q + 1));
			p.k = (x2p - x1p) / (cos(p.theta)*(x2-x1) + sin(p.theta)*(y2-y1));
			p.Tx = x1p - p.k*(cos(p.theta)*x1 + sin(p.theta)*y1);
			p.Ty = y1p - p.k*(-sin(p.theta)*x1 + cos(p.theta)*y1);
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
		Image raster = imclone(site);
		float box_x[4];
		float box_y[4];
		//write transformed coordinates to file
		fprintf(camera_location,"%d",numObjects[k]);
		for(int i=0;i<numObjects[k];i++) {
			float x = x_target[totalObjects + i];
			float y = y_target[totalObjects + i];
			float xp = bestParam.k*(cos(bestParam.theta)*x + sin(bestParam.theta)*y) + bestParam.Tx;
			float yp = bestParam.k*(-sin(bestParam.theta)*x + cos(bestParam.theta)*y) + bestParam.Ty;
			fprintf(camera_location," %f %f %d",xp,yp,match_target[totalObjects+i]);
			box_x[i%4] = xp;
			box_y[i%4] = yp;
			if (i%4 == 3)
				drawBox(&raster,box_x,box_y,red);
		}
		fprintf(camera_location,"\n");
		totalObjects += numObjects[k];
		sprintf(buffer,"%lu-site.ppm",k+1);
		printf("%s: %d/%lu inliers (%f)\n",buffer,maxInliers,x_src.size(),1.0*maxInliers/x_src.size());
		score += 1.0*maxInliers/x_src.size();
		writeImage(buffer,&raster,true);
		delete[] raster.data;
	}
	printf("PnP score: %f\n",score / numObjects.size());
	fclose(camera_location);
	delete[] site.data;
	return 1;

}
