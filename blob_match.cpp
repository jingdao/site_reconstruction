#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <vector>
#define K_PARAM 10
#define MIN_BLOB 1000
#define BLOB_RADIUS 5

typedef struct {
	int width,height;
	unsigned char* data;
} Image;

typedef struct {
	unsigned char r,g,b;
} Color;

typedef struct {
	int x,y;
} Coordinate;

void drawKeyPoint(Image image, int x, int y, int size, Color color) {
	for (int i=x-size+1;i<=x+size-1;i++) {
		for (int j=y-size+1;j<=y+size-1;j++) {
			if (i<0 || i>=image.width || j<0 || j>=image.height)
				continue;
			unsigned char* dst = image.data + (j*image.width + i) * 3; 
			*dst++ = color.r;
			*dst++ = color.g;
			*dst++ = color.b;
		}
	}
}

void writePPM(char* filename,unsigned char* imageBuffer,int width,int height) {
	FILE* f = fopen(filename,"w");
	if (!f)
		return;
	fprintf(f,"P6\n%d %d\n255\n",width,height);
	fwrite(imageBuffer,1,width*height*3,f);
	printf("Saved to %s\n",filename);
	fclose(f);
}

Image imclone(Image image) {
	Image im = {image.width,image.height,NULL};
	im.data = new unsigned char[im.width*im.height*3];
	memcpy(im.data,image.data,im.width*im.height*3);
	return im;
}

void convertLAB(Image image) {
	unsigned char *src = image.data;
	for (int i=0;i<image.height;i++) {
		for (int j=0;j<image.width;j++) {
			float r = (*src++) / 255.0;
			float g = (*src++) / 255.0;
			float b = (*src++) / 255.0;
			r = (r > 0.04045 ? pow((r + 0.055) / 1.055, 2.4) : r / 12.92) * 100.0;
			g = (g > 0.04045 ? pow((g + 0.055) / 1.055, 2.4) : g / 12.92) * 100.0;
			b = (b > 0.04045 ? pow((b + 0.055) / 1.055, 2.4) : b / 12.92) * 100.0;
			float X = r * 0.4124 + g * 0.3576 + b * 0.1805;
			float Y = r * 0.2126 + g * 0.7152 + b * 0.0722;
			float Z = r * 0.0193 + g * 0.1192 + b * 0.9505;

			float var_X = X / 95.047;
			float var_Y = Y / 100;
			float var_Z = Z / 108.883;

			if ( var_X > 0.008856 ) var_X = pow(var_X,( 1.0/3 ));
			else var_X = (903.3*var_X + 16) / 116;
			if ( var_Y > 0.008856 ) var_Y = pow(var_Y,( 1.0/3 ));
			else var_Y = (903.3*var_Y + 16) / 116;
			if ( var_Z > 0.008856 ) var_Z = pow(var_Z,( 1.0/3 ));
			else var_Z = (903.3*var_Z + 16) / 116;

			float L = 2.56 * (( 116 * var_Y ) - 16);
			float A = 1.388 * 500 * ( var_X - var_Y ) + 119.624;
			float B = 1.26 * 200 * ( var_Y - var_Z ) + 135.932;
//			printf("%f %f %f\n",L,A,B);
			src[-3] = L < 0 ? 0 : L > 255 ? 255 : (unsigned char) L;
			src[-2] = A < 0 ? 0 : A > 255 ? 255 : (unsigned char) A;
			src[-1] = B < 0 ? 0 : B > 255 ? 255 : (unsigned char) B;
		}
	}
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
	} else {
		fread(image->data,1,image->width*image->height,pgm);
		for (int i=image->width*image->height;i>=0;i--) {
			image->data[i*3+2] = image->data[i*3+1] = image->data[i*3] = image->data[i];
		}
	}
	return true;
}

int getDiff(Color c1, Color c2) {
	int d=0;
//	d += (c1.r - c2.r) * (c1.r - c2.r);
	d += (c1.g - c2.g) * (c1.g - c2.g);
	d += (c1.b - c2.b) * (c1.b - c2.b);
	return d;
}

void getKMeans(Image image, Color *palette, int k) {
	std::vector<Color> colors;
	for (int i=0;i<image.height;i++) {
		unsigned char* src = image.data + i * image.width * 3;
		for (int j=0;j<image.height;j++) {
			Color c;
			c.r = *src++;
			c.g = *src++;
			c.b = *src++;
			colors.push_back(c);
		}
	}
	if (colors.size()==0)
		return;
	for (int i=0;i<k;i++) {
		palette[i] = colors[rand() % colors.size()];
	}
	int* newR = new int[k];
	int* newG = new int[k];
	int* newB = new int[k];
	int* count = new int[k];
	int* match = new int[colors.size()]();
	bool updated = true;
	while (updated) {
		updated = false;
		memset(newR,0,k*sizeof(int));
		memset(newG,0,k*sizeof(int));
		memset(newB,0,k*sizeof(int));
		memset(count,0,k*sizeof(int));
		for (size_t i=0;i<colors.size();i++) {
			int minD=255*255*3,minID;
			for (int j=0;j<k;j++) {
				int d = getDiff(colors[i],palette[j]);
				if (d < minD) {
					minD = d;
					minID = j;
				}
			}
			newR[minID] += colors[i].r;
			newG[minID] += colors[i].g;
			newB[minID] += colors[i].b;
			count[minID]++;
			if (minID != match[i]) {
				match[i] = minID;
				updated = true;
			}
		}
		for (int j=0;j<k;j++) {
			if (count[j] == 0) {
				palette[j] = colors[rand() % colors.size()];
			} else {
				palette[j].r = newR[j] / count[j];
				palette[j].g = newG[j] / count[j];
				palette[j].b = newB[j] / count[j];
			}
		}
	}
	delete[] newR;
	delete[] newG;
	delete[] newB;
	delete[] count;
	delete[] match;
}

std::vector<Coordinate> getBlobs(Image image, Color *palette, int k, bool updateColor) {
	//label each pixel from 0 to k-1
	int** label = new int*[image.width];
	for (int i=0;i<image.width;i++)
		label[i] = new int[image.height]();
	unsigned char *src = image.data;
	for (int i=0;i<image.height;i++) {
		for (int j=0;j<image.width;j++) {
			Color c;
			c.r = *src++;
			c.g = *src++;
			c.b = *src++;
			int minD=255*255*3;
			for (int l=0;l<k;l++) {
				int d = getDiff(c,palette[l]);
				if (d < minD) {
					minD = d;
					label[j][i] = l;
				}
			}
			if (updateColor) {
				src[-3] = palette[label[j][i]].r;
				src[-2] = palette[label[j][i]].g;
				src[-1] = palette[label[j][i]].b;
			}
		}
	}
	//flood fill to detect blobs
	std::vector<Coordinate> blobs;
	bool** visited = new bool*[image.width];
	for (int i=0;i<image.width;i++)
		visited[i] = new bool[image.height]();
	for (int i=0;i<image.height;i++) {
		for (int j=0;j<image.width;j++) {
			if (visited[j][i])
				continue;
			Coordinate center = {};
			int numPoints = 0;
			int targetLabel = label[j][i];
			std::vector<Coordinate> stack;
			Coordinate p = {j,i};
			stack.push_back(p);
			while (stack.size() > 0) {
				Coordinate c = stack.back();
				stack.pop_back();
				if (c.x < 0 || c.x >= image.width || c.y < 0 || c.y >= image.height)
					continue;
				if (visited[c.x][c.y])
					continue;
				if (label[c.x][c.y] == targetLabel) {
					visited[c.x][c.y] = true;
					center.x += c.x;
					center.y += c.y;
					numPoints++;
					for (int m = -BLOB_RADIUS; m<= BLOB_RADIUS; m++) {
						for (int n = -BLOB_RADIUS; n <= BLOB_RADIUS; n++) {
							Coordinate r = {c.x + m, c.y + n};
							stack.push_back(r);
						}
					}
				}
			}
			if (numPoints > MIN_BLOB) {
				center.x /= numPoints;
				center.y /= numPoints;
				blobs.push_back(center);
			}
		}
	}
	for (int i=0;i<image.width;i++)
		delete[] visited[i];
	delete[] visited;
	for (int i=0;i<image.width;i++)
		delete[] label[i];
	delete[] label;
	return blobs;
}

int main(int argc, char* argv[]) {

	if (argc < 3) {
		printf("./%s source.pcd target.ppm\n",argv[0]);
		return 1;
	}
	srand(0);

	Image baseimage = {0,0,NULL};
	Image modimage;
	if (!loadImage(argv[2],&baseimage,true))
		return 1;
	Image labimage = imclone(baseimage);
	convertLAB(labimage);

	char buffer[256];
	Color yellow = {255,255,0};
	Color red = {255,0,0};
	Color blue = {0,0,255};
	std::vector< std::vector<Coordinate> > region;
	Color* palette = new Color[K_PARAM];

//	getKMeans(labimage,palette,K_PARAM);
	getKMeans(baseimage,palette,K_PARAM);
//	std::vector<Coordinate> blobs = getBlobs(labimage,palette,K_PARAM,false);
	std::vector<Coordinate> blobs = getBlobs(baseimage,palette,K_PARAM,true);
	for (size_t i=0;i<blobs.size();i++)
		drawKeyPoint(baseimage,blobs[i].x,blobs[i].y,3,red);
	sprintf(buffer,"out.ppm");
	writePPM(buffer,baseimage.data,baseimage.width,baseimage.height);

	delete[] palette;
	delete[] baseimage.data;
	delete[] labimage.data;

}
