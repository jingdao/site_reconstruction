#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <SDL/SDL.h>
#define FOREGROUND 1
#define K_PARAM 3
#define EXPAND 3.0
#define OUTLIER_RATIO 0.05
#define BOX_COLOR blue

SDL_Surface *screen;
bool mouseDrag = false;
int previousX,previousY;
int targetIndex = 0;

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

void drawRect(SDL_Surface *surf,SDL_Rect rect,Color color) {
	int top = rect.y < 0 ? 0 : rect.y ;
	int bottom = rect.y + rect.h >= surf->h ? surf->h-1 : rect.y + rect.h;
	int left = rect.x < 0 ? 0 : rect.x;
	int right = rect.x + rect.w >= surf->w ? surf->w-1 : rect.x + rect.w;
	//TOP
	int i=left;
	unsigned char* c = (unsigned char*)surf->pixels + top * surf->pitch + left*3;
	while (i++ < right) {
		*c++ = color.b;
		*c++ = color.g;
		*c++ = color.r;
	}	
	//BOTTOM
	i=left;
	c = (unsigned char*)surf->pixels + bottom * surf->pitch + left*3;
	while (i++ < right) {
		*c++ = color.b;
		*c++ = color.g;
		*c++ = color.r;
	}	
	//LEFT
	i=top;
	c = (unsigned char*)surf->pixels + top * surf->pitch + left*3;
	while (i++ < bottom) {
		c[0] = color.b;
		c[1] = color.g;
		c[2] = color.r;
		c += surf->pitch;
	}	
	//RIGHT
	i=top;
	c = (unsigned char*)surf->pixels + top * surf->pitch + right*3;
	while (i++ < bottom) {
		c[0] = color.b;
		c[1] = color.g;
		c[2] = color.r;
		c += surf->pitch;
	}
}

void drawLine(SDL_Surface* surf, float x1, float y1, float x2, float y2, Color color ) {
	if (x1 < 1) x1 = 1;
	if (x1 > surf->w - 1) x1 = surf->w - 1;
	if (x2 < 1) x2 = 1;
	if (x2 > surf->w - 1) x2 = surf->w - 1;
	if (y1 < 1) y1 = 1;
	if (y1 > surf->h - 1) y1 = surf->h - 1;
	if (y2 < 1) y2 = 1;
	if (y2 > surf->h - 1) y2 = surf->h - 1;
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
			dst = (unsigned char*)surf->pixels + x * surf->pitch + y*3;
		else
			dst = (unsigned char*)surf->pixels + y * surf->pitch + x*3;
		dst[0] = color.b;
		dst[1] = color.g;
		dst[2] = color.r;
		error -= dy;
		if(error < 0) {
			y += ystep;
			error += dx;
		}
	}
}

void drawBox(SDL_Surface *surf, float* box, Color color) {
	drawLine(surf,box[0],box[1],box[2],box[3],color);
	drawLine(surf,box[0],box[1],box[4],box[5],color);
	drawLine(surf,box[2],box[3],box[6],box[7],color);
	drawLine(surf,box[4],box[5],box[6],box[7],color);
}

SDL_Rect getVarRect(int x1,int y1,int x2,int y2) {
	int left = x1 < x2 ? x1 : x2;
	int right = x1 > x2 ? x1 : x2;
	int top = y1 < y2 ? y1 : y2;
	int bottom = y1 > y2 ? y1 : y2;
	SDL_Rect r = {(short)left,(short)top,(unsigned short)(right-left),(unsigned short)(bottom-top)};
	return r;
}

void imgcpy(Image image,SDL_Surface* surf) {
	unsigned char* src = image.data, *dst = (unsigned char*)surf->pixels;
	for (int i=0;i<image.height;i++) {
		memcpy(dst,src,image.width*3);
		src += image.width*3;
		dst += surf->pitch;
	}
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

int getDiff(Color c1, Color c2) {
	int d=0;
//	d += (c1.r - c2.r) * (c1.r - c2.r);
	d += (c1.g - c2.g) * (c1.g - c2.g);
	d += (c1.b - c2.b) * (c1.b - c2.b);
	return d;
}

int sortByX(const void* v1, const void* v2) {
	Coordinate* c1 = (Coordinate*)v1;
	Coordinate* c2 = (Coordinate*)v2;
	return c1->x - c2->x;
}

int sortByY(const void* v1, const void* v2) {
	Coordinate* c1 = (Coordinate*)v1;
	Coordinate* c2 = (Coordinate*)v2;
	return c1->y - c2->y;
}

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

unsigned char getThreshold(Image image, SDL_Rect* rect) {
	std::vector<unsigned char> intensity;
	unsigned char min=255,max=0;
	unsigned char* src = image.data + rect->y * image.width + rect->x;
	for (int i=0;i<rect->h;i++) {
		for (int j=0;j<rect->w;j++) {
			unsigned char val = *src++;
			if (val < min) min = val;
			else if (val > max) max = val;
			intensity.push_back(val);
		}
		src += image.width - 1;
	}
	if (intensity.size()==0)
		return 0;
	std::sort(intensity.begin(),intensity.end());
	unsigned char median = intensity[intensity.size()/2];
	return median;
}

unsigned char getKmeansThreshold(Image image, SDL_Rect* rect) {
	std::vector<int> intensity;
	unsigned char* src = image.data + rect->y * image.width + rect->x;
	for (int i=0;i<rect->h;i++) {
		for (int j=0;j<rect->w;j++) {
			unsigned char val = *src++;
			intensity.push_back(val);
		}
		src += image.width - 1;
	}
	if (intensity.size()==0)
		return 0;
	int low = rand() % 128;
	int high = rand() % 128 + 128;
	int mid = (low + high) / 2;
	int newLow,newHigh,numLow,numHigh;
	for (int n=0;n<5;n++) {
		newLow = 0;
		newHigh = 0;
		numLow = 0;
		numHigh = 0;
		for (size_t i=0;i<intensity.size();i++) {
			if (intensity[i] < mid) {
				newLow += intensity[i];
				numLow++;
			} else {
				newHigh += intensity[i];
				numHigh++;
			}
		}
		low = newLow / numLow;
		high = newHigh / numHigh;
		mid = (low + high) / 2;
	}
	return mid;
}

void getKMeans(Image image, SDL_Rect* rect, Color *palette, int k) {
	std::vector<Color> colors;
	for (int i=0;i<rect->h;i++) {
		unsigned char* src = image.data + ((rect->y+i) * image.width + rect->x)*3;
		for (int j=0;j<rect->w;j++) {
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
//		palette[i].r = rand() % 256;
//		palette[i].g = rand() % 256;
//		palette[i].b = rand() % 256;
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
//				palette[j].r = rand() % 256;
//				palette[j].g = rand() % 256;
//				palette[j].b = rand() % 256;
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

int findTargetColor(Image image, SDL_Rect* rect, Color *palette, int k) {
	int* count = new int[k]();
	int maxCount=0,target;
	int dh = rect->h/2;
	int dw = rect->w/2;
	for (int i=0;i<dh;i++) {
		unsigned char* src = image.data + ((rect->y+(rect->h-dh)/2+i) * image.width + (rect->w-dw)/2+rect->x)*3;
		for (int j=0;j<dw;j++) {
			Color c;
			c.r = *src++;
			c.g = *src++;
			c.b = *src++;
			int minD=255*255*3,minID;
			for (int j=0;j<k;j++) {
				int d = getDiff(c,palette[j]);
				if (d < minD) {
					minD = d;
					minID = j;
				}
			}
			count[minID]++;
			if (count[minID] > maxCount) {
				maxCount = count[minID];
				target = minID;
			}
		}
	}
	delete[] count;
	return target;
}

int findBackgroundColor(Image image, SDL_Rect* rect, Color *palette, int k) {
	int* count = new int[k]();
	int maxCount=0,target;
	for (int i=0;i<rect->h;i++) {
		unsigned char* src = image.data + ((rect->y+i) * image.width + rect->x)*3;
		for (int j=0;j<rect->w;j++) {
			Color c;
			c.r = *src++;
			c.g = *src++;
			c.b = *src++;
			int minD=255*255*3,minID;
			for (int j=0;j<k;j++) {
				int d = getDiff(c,palette[j]);
				if (d < minD) {
					minD = d;
					minID = j;
				}
			}
			count[minID]++;
			if (count[minID] > maxCount) {
				maxCount = count[minID];
				target = minID;
			}
		}
	}
	delete[] count;
	return target;
}

void regionGrow(Image image, std::vector<Coordinate> *region, Coordinate seed, unsigned char threshold) {
	region->clear();
	bool** visited = new bool*[image.width];
	for (int i=0;i<image.width;i++)
		visited[i] = new bool[image.height]();
	std::vector<Coordinate> stack;
	while (image.data[(seed.y*image.width+seed.x)*3] > threshold) {
		seed.x += rand() % 5 - 2;
		seed.y += rand() % 5 - 2;
	}
	stack.push_back(seed);
	while (stack.size() > 0) {
		Coordinate c = stack.back();
		stack.pop_back();
		if (c.x < 0 || c.x >= image.width || c.y < 0 || c.y >= image.height)
			continue;
		if (visited[c.x][c.y])
			continue;
		visited[c.x][c.y] = true;
//		printf("%d %d %hhu\n",c.x,c.y,image.data[c.y*image.width+c.x]);
		if (image.data[(c.y*image.width+c.x)*3] > threshold)
			continue;
		region->push_back(c);
		Coordinate left = {c.x - 1, c.y};
		Coordinate right = {c.x + 1, c.y};
		Coordinate top = {c.x, c.y - 1};
		Coordinate bottom = {c.x, c.y + 1};
		stack.push_back(left);
		stack.push_back(right);
		stack.push_back(top);
		stack.push_back(bottom);
	}
	for (int i=0;i<image.width;i++)
		delete[] visited[i];
	delete[] visited;
}

void regionGrow2(Image image, std::vector<Coordinate> *region, Coordinate seed, int threshold) {
	region->clear();
	bool** visited = new bool*[image.width];
	for (int i=0;i<image.width;i++)
		visited[i] = new bool[image.height]();
	std::vector<Coordinate> stack;
	stack.push_back(seed);
	int seed_intensity = image.data[(seed.y*image.width+seed.x)*3];
	while (stack.size() > 0) {
		Coordinate c = stack.back();
		stack.pop_back();
		if (c.x < 0 || c.x >= image.width || c.y < 0 || c.y >= image.height)
			continue;
		if (visited[c.x][c.y])
			continue;
		visited[c.x][c.y] = true;
		int intensity = image.data[(c.y*image.width+c.x)*3];
		int distance = (intensity - seed_intensity) * (intensity - seed_intensity);
//		printf("%d %d %d %d\n",c.x,c.y,seed_intensity,intensity);
//		if (image.data[c.y*image.width+c.x] > threshold)
		if (distance > threshold)
			continue;
		region->push_back(c);
		Coordinate left = {c.x - 1, c.y};
		Coordinate right = {c.x + 1, c.y};
		Coordinate top = {c.x, c.y - 1};
		Coordinate bottom = {c.x, c.y + 1};
		stack.push_back(left);
		stack.push_back(right);
		stack.push_back(top);
		stack.push_back(bottom);
	}
	for (int i=0;i<image.width;i++)
		delete[] visited[i];
	delete[] visited;
}

void colorRegion(Image image, std::vector<Coordinate> *region, SDL_Rect* rect, Color* palette, int k, int targetID) {
	region->clear();
	for (int i=0;i<rect->h;i++) {
		unsigned char* src = image.data + ((rect->y+i) * image.width + rect->x)*3;
		for (int j=0;j<rect->w;j++) {
			Color c;
			c.r = *src++;
			c.g = *src++;
			c.b = *src++;
			int minD=255*255*3,minID;
			for (int l=0;l<k;l++) {
				int d = getDiff(c,palette[l]);
				if (d < minD) {
					minD = d;
					minID = l;
				}
			}
			if (minID == targetID) {
				Coordinate r = {j+rect->x,i+rect->y};
				region->push_back(r);
			}
		}
	}
}

void colorRegionInBox(Image image, std::vector<Coordinate> *region, float* box, Color* palette, int k, int targetID) {
	region->clear();
	float minX = box[0];
	float maxX = box[0];
	float minY = box[1];
	float maxY = box[1];
	for (int i=1;i<4;i++) {
		float x = box[i*2];
		float y = box[i*2+1];
		if (x < minX) minX = x;
		if (x > maxX) maxX = x;
		if (y < minY) minY = y;
		if (y > maxY) maxY = y;
	}
	if (minX < 0) minX = 0;
	if (maxX > image.width - 1) maxX = image.width - 1;
	if (minY < 0) minY = 0;
	if (maxY > image.height - 1) maxY = image.height - 1;
	float ux = box[2] - box[0];
	float uy = box[3] - box[1];
	float vx = box[4] - box[0];
	float vy = box[5] - box[1];
	float mu = ux*ux + uy*uy;
	float mv = vx*vx + vy*vy;
	for (int i=(int)minY; i<= (int)maxY; i++ ) {
		unsigned char* src = image.data + (i * image.width + (int)minX) *3;
		for (int j=(int)minX; j<= (int)maxX; j++ ) {
			Color c;
			c.r = *src++;
			c.g = *src++;
			c.b = *src++;
			float wu = ((j - box[0])*ux + (i - box[1])*uy) / mu;
			float wv = ((j - box[0])*vx + (i - box[1])*vy) / mv;
			if (wu>0 && wu<1 && wv>0 && wv<1) {
				int minD=255*255*3,minID;
				for (int l=0;l<k;l++) {
					int d = getDiff(c,palette[l]);
					if (d < minD) {
						minD = d;
						minID = l;
					}
				}
				if (minID == targetID) {
					Coordinate r = {j,i};
					region->push_back(r);
				}
			}
		}
	}

}

void floodFillRegion(Image image, std::vector<Coordinate> *region, SDL_Rect* rect, Color* palette, int k, int targetID) {
	region->clear();
	bool** visited = new bool*[rect->w];
	for (int i=0;i<rect->w;i++)
		visited[i] = new bool[rect->h]();
	for (int i=0;i<rect->h;i++) {
		for (int j=0;j<rect->w;j++) {
			if (visited[j][i])
				continue;
			std::vector<Coordinate> marked;
			std::vector<Coordinate> stack;
			Coordinate p = {j+rect->x,i+rect->y};
			stack.push_back(p);
			while (stack.size() > 0) {
				Coordinate c = stack.back();
				int dx = c.x - rect->x;
				int dy = c.y - rect->y;
				stack.pop_back();
				if (dx < 0 || dx >= rect->w || dy < 0 || dy >= rect->h)
					continue;
				if (visited[dx][dy])
					continue;
				visited[dx][dy] = true;
				Color q;
				q.r = image.data[(c.y*image.width+c.x)*3];
				q.g = image.data[(c.y*image.width+c.x)*3+1];
				q.b = image.data[(c.y*image.width+c.x)*3+2];
				int minD=255*255*3,minID;
				for (int l=0;l<k;l++) {
					int d = getDiff(q,palette[l]);
					if (d < minD) {
						minD = d;
						minID = l;
					}
				}
#if FOREGROUND
				if (minID == targetID) {
#else
				if (minID != targetID) {
#endif
					marked.push_back(c);
					for (int m = -2; m<= 2; m++) {
						for (int n = -2; n <= 2; n++) {
							Coordinate r = {c.x + m, c.y + n};
							stack.push_back(r);
						}
					}
				}
			}
			if (marked.size() > region->size())
				*region = marked;
		}
	}
	for (int i=0;i<rect->w;i++)
		delete[] visited[i];
	delete[] visited;

}

void floodFillRegionInBox(Image image, std::vector<Coordinate> *region, float* box, Color* palette, int k, int targetID) {
	region->clear();
	float minX = box[0];
	float maxX = box[0];
	float minY = box[1];
	float maxY = box[1];
	for (int i=1;i<4;i++) {
		float x = box[i*2];
		float y = box[i*2+1];
		if (x < minX) minX = x;
		if (x > maxX) maxX = x;
		if (y < minY) minY = y;
		if (y > maxY) maxY = y;
	}
	if (minX < 0) minX = 0;
	if (maxX > image.width - 1) maxX = image.width - 1;
	if (minY < 0) minY = 0;
	if (maxY > image.height - 1) maxY = image.height - 1;
	float ux = box[2] - box[0];
	float uy = box[3] - box[1];
	float vx = box[4] - box[0];
	float vy = box[5] - box[1];
	float mu = ux*ux + uy*uy;
	float mv = vx*vx + vy*vy;
	int width = (int)maxX - (int)minX;
	int height = (int)maxY - (int)minY;
	bool** visited = new bool*[width];
	for (int i=0;i<width;i++)
		visited[i] = new bool[height]();
	for (int i=0;i<height;i++) {
		for (int j=0;j<width;j++) {
			if (visited[j][i])
				continue;
			std::vector<Coordinate> marked;
			std::vector<Coordinate> stack;
			Coordinate p = {j+(int)minX,i+(int)minY};
			stack.push_back(p);
			while (stack.size() > 0) {
				Coordinate c = stack.back();
				stack.pop_back();
				int dx = c.x - (int)minX;
				int dy = c.y - (int)minY;
				if (dx < 0 || dx >= width || dy < 0 || dy >= height)
					continue;
				if (visited[dx][dy])
					continue;
				visited[dx][dy] = true;
				float wu = ((c.x - box[0])*ux + (c.y - box[1])*uy) / mu;
				float wv = ((c.x - box[0])*vx + (c.y - box[1])*vy) / mv;
				if (wu<0 || wu>1 || wv<0 || wv>1)
					continue;
				Color q;
				q.r = image.data[(c.y*image.width+c.x)*3];
				q.g = image.data[(c.y*image.width+c.x)*3+1];
				q.b = image.data[(c.y*image.width+c.x)*3+2];
				int minD=255*255*3,minID;
				for (int l=0;l<k;l++) {
					int d = getDiff(q,palette[l]);
					if (d < minD) {
						minD = d;
						minID = l;
					}
				}
#if FOREGROUND
				if (minID == targetID) {
#else
				if (minID != targetID) {
#endif
					marked.push_back(c);
					for (int m = -2; m<= 2; m++) {
						for (int n = -2; n <= 2; n++) {
							Coordinate r = {c.x + m, c.y + n};
							stack.push_back(r);
						}
					}
				}
			}
			if (marked.size() > region->size())
				*region = marked;
		}
	}
	for (int i=0;i<width;i++)
		delete[] visited[i];
	delete[] visited;

}

void highlightRegion(SDL_Surface *surf, std::vector<Coordinate> *region, Color color) {
	for (size_t i=0;i<region->size();i++) {
		unsigned char* c = (unsigned char*)surf->pixels + region->at(i).y * surf->pitch + region->at(i).x*3;
		*c++ = color.b;
		*c++ = color.g;
		*c = color.r;
	}
}

void binarize(Image image, Image output, unsigned int threshold) {
	unsigned char* src = image.data;
	unsigned char* dst = output.data;
	for (int i=0;i<image.width*image.height;i++) {
		if (*src > threshold) {
			*dst++ = 0xFF;
			*dst++ = 0xFF;
			*dst++ = 0xFF;
		} else {
			*dst++ = 0;
			*dst++ = 0;
			*dst++ = 0;
		}
		src += 3;
	}
}

void sobel(Image image, int threshold) {
	int* src = new int[image.width*image.height];
	int width = image.width;
	int height = image.height;
	for (int i=0;i<image.width*image.height;i++)
		src[i] = image.data[i*3];
	for (int i=1;i<height-1;i++) {
		for (int j=1;j<width-1;j++) {
			int gx=0,gy=0;
			gx -= 2 * src[i*width+j-1];
			gx += 2 * src[i*width+j+1];
			gx -= src[(i-1)*width+j-1];
			gx += src[(i-1)*width+j+1];
			gx -= src[(i+1)*width+j-1];
			gx += src[(i+1)*width+j+1];
			gy -= 2 * src[(i-1)*width+j];
			gy += 2 * src[(i+1)*width+j];
			gy -= src[(i-1)*width+j-1];
			gy += src[(i+1)*width+j-1];
			gy -= src[(i-1)*width+j+1];
			gy += src[(i+1)*width+j+1];
			unsigned char* dst = image.data + (i*width+j)*3;
			if (abs(gx) + abs(gy) > threshold) {
				*dst++ = 0xFF;
				*dst++ = 0xFF;
				*dst++ = 0xFF;
			} else {
				*dst++ = 0;
				*dst++ = 0;
				*dst++ = 0;
			}
		}
	}
	delete[] src;
}

void mooreTracing(Image image, std::vector<Coordinate> *region, SDL_Rect *rect) {
	region->clear();
	bool* border = new bool[image.width*image.height];
	bool inside = false;
	int pos = 0;
	int neighborhood[8][2] = {
		{-1,6},
		{-3-image.width,6},
		{-2-image.width,0},
		{-1-image.width,0},
		{1,2},
		{3+image.width,2},
		{2+image.width,4},
		{1+image.width,4}
	};
	for (int y=rect->y+1;y<rect->y+rect->h-1;y++) {
		for (int x=rect->x+1;x<rect->x+rect->w-1;x++) {
			pos = x + y*image.width;
			if (border[pos] && !inside)
				inside = true;
			else if (image.data[pos*3] == 0xFF && inside)
				continue;
			else if (image.data[pos*3] == 0 && inside)
				inside = false;
			else if (image.data[pos*3] == 0xFF && !inside) {
				border[pos] = true;
				Coordinate c = {pos%image.width, pos/image.width};
				region->push_back(c);
				int checkId = 0;
				int checkPos,nextId,startPos = pos;
				int stop=0,numVisited=0;
				while (true) {
					checkPos = pos + neighborhood[checkId][0];
					nextId = neighborhood[checkId][1];
					if (image.data[checkPos*3] == 0xFF) {
						if (checkPos == startPos) {
							stop++;
							if (nextId==0 || stop >= 3) {
								inside = true;
								break;
							}
						}
						checkId = nextId;
						pos = checkPos;
						numVisited = 0;
						border[checkPos] = true;
						Coordinate c = {checkPos%image.width, checkPos/image.width};
						region->push_back(c);
					} else {
						checkId = (checkId + 1) % 8;
						if (numVisited > 8) {
							numVisited = 0;
							break;
						} else 
							numVisited++;
					}
				}
			}
		}
	}
	delete[] border;
}

void getPCA(std::vector<Coordinate> *cloud, float *box) {
	if (cloud->size() == 0)
		return;
	qsort(cloud->data(),cloud->size(),sizeof(Coordinate),sortByX);
	qsort(cloud->data()+(int)(cloud->size()*OUTLIER_RATIO),
			cloud->size()-(int)(cloud->size()*OUTLIER_RATIO),sizeof(Coordinate),sortByY);
	int inlier_size = cloud->size() - (int)(cloud->size()*OUTLIER_RATIO*4);
	double cov[4] = {}; //column major
	Coordinate center = {};
	for (int i = cloud->size()*OUTLIER_RATIO*2; i < cloud->size()-(int)(cloud->size()*OUTLIER_RATIO*2); i++) {
		center.x += cloud->at(i).x;
		center.y += cloud->at(i).y;
	}
	center.x /= inlier_size; 
	center.y /= inlier_size;
	for (int j = cloud->size()*OUTLIER_RATIO*2; j<cloud->size()-(int)(cloud->size()*OUTLIER_RATIO*2); j++) {
		float deltaP[2] = {
			(float) cloud->at(j).x- center.x,
			(float) cloud->at(j).y- center.y,
		};
		cov[0] += deltaP[0] * deltaP[0];
		cov[1] += deltaP[1] * deltaP[0];
		cov[2] += deltaP[0] * deltaP[1];
		cov[3] += deltaP[1] * deltaP[1];
	}
	cov[0] /= inlier_size * inlier_size;
	cov[1] /= inlier_size * inlier_size;
	cov[2] /= inlier_size * inlier_size;
	cov[3] /= inlier_size * inlier_size;
	float trace = cov[0] + cov[3];
	float det = cov[0] * cov[3] - cov[1] * cov[2];
	float L1 = trace / 2 + sqrt(trace*trace / 4 - det);
	float L2 = trace / 2 - sqrt(trace*trace / 4 - det);
	float minScale[2], maxScale[2];
	float v[4] = {};
	if (cov[2] != 0) {
		v[0] = L1 - cov[3];
		v[1] = L2 - cov[3];
		v[2] = v[3] = cov[2];
	}
	else if (cov[1] != 0) {
		v[0] = v[1] = cov[1];
		v[2] = L1 - cov[0];
		v[3] = L2 - cov[0];
	}
	else {
		v[0] = v[3] = 1;
	}
	float m1 = sqrt(v[0] * v[0] + v[2] * v[2]);
	float m2 = sqrt(v[1] * v[1] + v[3] * v[3]);
	v[0] /= m1;
	v[2] /= m1;
	v[1] /= m2;
	v[3] /= m2;
	for (int j = cloud->size()*OUTLIER_RATIO*2; j<cloud->size()-(int)(cloud->size()*OUTLIER_RATIO*2); j++) {
		for (int i = 0; i<2; i++) {
			float dotProduct =
				cloud->at(j).x * v[i * 2] +
				cloud->at(j).y * v[i * 2 + 1];
			if (j == (int)(cloud->size()*OUTLIER_RATIO*2) || dotProduct < minScale[i])
				minScale[i] = dotProduct;
			if (j == (int)(cloud->size()*OUTLIER_RATIO*2) || dotProduct > maxScale[i])
				maxScale[i] = dotProduct;
		}
	}
	float bbCenter[2] = {0,0};
	for (int i = 0; i<2; i++) {
		bbCenter[0] += (minScale[i] + maxScale[i]) / 2 * v[i * 2];
		bbCenter[1] += (minScale[i] + maxScale[i]) / 2 * v[i * 2 + 1];
	}
	for (int i = 0; i<4; i++) {
		float coords[2];
		for (int j = 0; j<2; j++) {
			coords[j] = bbCenter[j];
			for (int axis = 0; axis<2; axis++) {
				float sign = (i & 1 << axis) ? 1 : -1;
				coords[j] += sign * (maxScale[axis]-minScale[axis]) / 2 * v[axis * 2 + j];
			}
			box[i*2 + j] = coords[j];
		}
	}
}

void expandBox(float* box, float expand_u, float expand_v) {
	float x0 = box[0];
	float y0 = box[1];
	float ux = box[2] - box[0];
	float uy = box[3] - box[1];
	float vx = box[4] - box[0];
	float vy = box[5] - box[1];
	box[0] = x0 + (1-expand_u)/2 * ux + (1-expand_v)/2 * vx;
	box[1] = y0 + (1-expand_u)/2 * uy + (1-expand_v)/2 * vy;
	box[2] = x0 + (1+expand_u)/2 * ux + (1-expand_v)/2 * vx;
	box[3] = y0 + (1+expand_u)/2 * uy + (1-expand_v)/2 * vy;
	box[4] = x0 + (1-expand_u)/2 * ux + (1+expand_v)/2 * vx;
	box[5] = y0 + (1-expand_u)/2 * uy + (1+expand_v)/2 * vy;
	box[6] = x0 + (1+expand_u)/2 * ux + (1+expand_v)/2 * vx;
	box[7] = y0 + (1+expand_u)/2 * uy + (1+expand_v)/2 * vy;
}

void minBoundRect(std::vector<Coordinate> *cloud, float *box, bool keepSize) {
	if (cloud->size() == 0)
		return;
	std::vector<float> edgeAngles;
	for (size_t i=0;i<cloud->size()-1;i++) {
		float theta = atan2(cloud->at(i+1).y-cloud->at(i).y,cloud->at(i+1).x-cloud->at(i).x);
		edgeAngles.push_back(theta);
	}
	float ux = box[2] - box[0];
	float uy = box[3] - box[1];
	float vx = box[4] - box[0];
	float vy = box[5] - box[1];
	float mu = sqrt(ux*ux + uy*uy);
	float mv = sqrt(vx*vx + vy*vy);
	float scale_u = mu / EXPAND;
	float scale_v = mv / EXPAND;
	float minArea;
	for (size_t i=0;i<edgeAngles.size();i++) {
		float theta = edgeAngles[i];
		float R[4] = {(float)cos(theta), (float)sin(theta), (float)-sin(theta), (float)cos(theta)};
		float xmin,xmax,ymin,ymax;
		for (size_t j=0;j<cloud->size();j++) {
			float x = R[0] * cloud->at(j).x + R[1] * cloud->at(j).y;
			float y = R[2] * cloud->at(j).x + R[3] * cloud->at(j).y;
			if (j==0 || x < xmin) xmin = x;
			if (j==0 || x > xmax) xmax = x;
			if (j==0 || y < ymin) ymin = y;
			if (j==0 || y > ymax) ymax = y;
		}
		float area = (xmax-xmin) * (ymax-ymin);
		if (i==0 || area < minArea) {
			minArea = area;
			box[0] = R[0] * xmin + R[2] * ymin;
			box[1] = R[1] * xmin + R[3] * ymin;
			box[2] = R[0] * xmax + R[2] * ymin;
			box[3] = R[1] * xmax + R[3] * ymin;
			box[4] = R[0] * xmin + R[2] * ymax;
			box[5] = R[1] * xmin + R[3] * ymax;
			box[6] = R[0] * xmax + R[2] * ymax;
			box[7] = R[1] * xmax + R[3] * ymax;
		}
	}
	if (keepSize) {
		ux = box[2] - box[0];
		uy = box[3] - box[1];
		vx = box[4] - box[0];
		vy = box[5] - box[1];
		mu = sqrt(ux*ux + uy*uy);
		mv = sqrt(vx*vx + vy*vy);
		if (mu<mv && scale_u<scale_v || mu>mv && scale_u>scale_v)
			expandBox(box,scale_u / mu, scale_v / mv);
		else
			expandBox(box,scale_v / mu, scale_u / mv);
	}
}

void writeBoxToFile(std::vector<float> *box, FILE* f) {
	fprintf(f,"%lu",box->size() / 2);
	for (size_t i=0;i<box->size();i++) {
		fprintf(f," %f",box->at(i));
		if (i%2 == 1)
			fprintf(f," 1");
	}
	fprintf(f,"\n");
}

SDL_Rect boxToRect(SDL_Surface *surf, float *box,float expand) {
	float minX = box[0];
	float maxX = box[0];
	float minY = box[1];
	float maxY = box[1];
	for (int i=1;i<4;i++) {
		float x = box[i*2];
		float y = box[i*2+1];
		if (x < minX) minX = x;
		if (x > maxX) maxX = x;
		if (y < minY) minY = y;
		if (y > maxY) maxY = y;
	}
	SDL_Rect r = {
		(short) ((1+expand)/2 * minX + (1-expand)/2 * maxX),
		(short) ((1+expand)/2 * minY + (1-expand)/2 * maxY),
		(unsigned short) (expand * (maxX - minX)),
		(unsigned short) (expand * (maxY - minY))
	};
	if (r.x < 0) r.x = 0;
	if (r.x >= surf->w ) r.x = surf->w - 1;
	if (r.y < 0) r.y = 0;
	if (r.y >= surf->h ) r.y = surf->h - 1;
	if (r.x + r.w >= surf->w) r.w = surf->w - 1 - r.x;
	if (r.y + r.h >= surf->h) r.h = surf->h - 1 - r.y;
	return r;
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

bool loadImageByIndex(int index, Image *image,bool color) {
	char buffer[128];
	if (color)
		sprintf(buffer,"%d.ppm",index);
	else
		sprintf(buffer,"%d.pgm",index);
	return loadImage(buffer,image,color);
}

void writeImageByIndex(int index,SDL_Surface *surf,bool color) {
	char buffer[128];
	if (color)
		sprintf(buffer,"%d-seg.ppm",index);
	else
		sprintf(buffer,"%d-seg.pgm",index);
	FILE* f = fopen(buffer,"w");
	if (color) {
		fprintf(f,"P6\n%d %d\n255\n",surf->w,surf->h);
		unsigned char* src = (unsigned char*) surf->pixels;
		for (int i=0;i<surf->h;i++) {
			fwrite(src,1,surf->w*3,f);
			src += surf->pitch;
		}
	} else {
		fprintf(f,"P5\n%d %d\n255\n",surf->w,surf->h);
		unsigned char* src = (unsigned char*) surf->pixels;
		for (int i=0;i<surf->h;i++) {
			fwrite(src,1,surf->w,f);
			src += surf->pitch;
		}
	}
	fclose(f);
}

int main(int argc, char* argv[]) {

	if (argc < 3) {
		printf("%s target.ppm target_point.txt [1.ppm ..]\n",argv[0]);
		return 1;
	}
	srand(0);

	SDL_Init(SDL_INIT_VIDEO);
	SDL_WM_SetCaption("imsegment",NULL);
	Image baseimage = {0,0,NULL};
	Image modimage;
	bool use_color = strncmp(argv[1]+strlen(argv[1])-4,".ppm",4)==0;
	if (!loadImage(argv[1],&baseimage,use_color))
		return 1;
	Image labimage = imclone(baseimage);
	convertLAB(labimage);
	screen = SDL_SetVideoMode(baseimage.width,baseimage.height,24,SDL_SWSURFACE);

	float cx = 0.5 * (baseimage.width-1);
	float cy = 0.5 * (baseimage.height-1);

	Color yellow = {255,255,0};
	Color red = {255,0,0};
	Color blue = {0,0,255};
	std::vector<SDL_Rect> rectList;
	SDL_Rect currentRect;
	std::vector< std::vector<Coordinate> > region;
	std::vector<Color> palette;
	std::vector<float> box;
	std::vector<int> targetList;
	FILE* target_point = fopen(argv[2],"w");

	imgcpy(baseimage,screen);
	SDL_Flip(screen);

	SDL_Event event;
	while (true) {
		while (SDL_PollEvent(&event)) {
			switch(event.type){
				case SDL_KEYDOWN:
					switch( event.key.keysym.sym ){
						case 'b':
						if (loadImageByIndex(--targetIndex,&baseimage,use_color)) {
							delete[] labimage.data;
							labimage = imclone(baseimage);
							convertLAB(labimage);
							imgcpy(baseimage,screen);
							for (size_t i=0;i<rectList.size();i++) {
								float* currentBox = box.data() + i * 8;
								Color* currentPalette = palette.data() + i * K_PARAM;
								std::vector<Coordinate> *currentRegion = region.data() + i;
								int target = targetList[i];
								expandBox(currentBox,EXPAND,EXPAND);
								floodFillRegionInBox(labimage,currentRegion,currentBox,currentPalette,K_PARAM,target);
								highlightRegion(screen,currentRegion,BOX_COLOR);
								*currentRegion = convexHull(*currentRegion);
								minBoundRect(currentRegion,currentBox,true);
								drawBox(screen,currentBox,BOX_COLOR);
							}
							writeBoxToFile(&box,target_point);
							writeImageByIndex(targetIndex,screen,use_color);
							SDL_Flip(screen);
						}
						break;
						case 'n':
						if (loadImageByIndex(++targetIndex,&baseimage,use_color)) {
							delete[] labimage.data;
							labimage = imclone(baseimage);
							convertLAB(labimage);
							imgcpy(baseimage,screen);
							for (size_t i=0;i<rectList.size();i++) {
								float* currentBox = box.data() + i * 8;
								Color* currentPalette = palette.data() + i * K_PARAM;
								std::vector<Coordinate> *currentRegion = region.data() + i;
								int target = targetList[i];
								expandBox(currentBox,EXPAND,EXPAND);
								floodFillRegionInBox(labimage,currentRegion,currentBox,currentPalette,K_PARAM,target);
								highlightRegion(screen,currentRegion,BOX_COLOR);
								*currentRegion = convexHull(*currentRegion);
								minBoundRect(currentRegion,currentBox,true);
								drawBox(screen,currentBox,BOX_COLOR);
							}
							writeBoxToFile(&box,target_point);
							writeImageByIndex(targetIndex,screen,use_color);
							SDL_Flip(screen);
						}
						break;
						case 'm': //bulk process
                        struct timespec tic,toc;
                        clock_gettime(CLOCK_MONOTONIC,&tic);
						while (loadImageByIndex(++targetIndex,&baseimage,use_color)) {
							delete[] labimage.data;
							labimage = imclone(baseimage);
							convertLAB(labimage);
							imgcpy(baseimage,screen);
							for (size_t i=0;i<rectList.size();i++) {
								float* currentBox = box.data() + i * 8;
								Color* currentPalette = palette.data() + i * K_PARAM;
								std::vector<Coordinate> *currentRegion = region.data() + i;
								int target = targetList[i];
								expandBox(currentBox,EXPAND,EXPAND);
								floodFillRegionInBox(labimage,currentRegion,currentBox,currentPalette,K_PARAM,target);
								highlightRegion(screen,currentRegion,BOX_COLOR);
								*currentRegion = convexHull(*currentRegion);
								minBoundRect(currentRegion,currentBox,true);
								drawBox(screen,currentBox,BOX_COLOR);
							}
							writeBoxToFile(&box,target_point);
							writeImageByIndex(targetIndex,screen,use_color);
							SDL_Flip(screen);
						}
                        clock_gettime(CLOCK_MONOTONIC,&toc);
                        printf("imsegment: %fs\n",toc.tv_sec - tic.tv_sec + 0.000000001 * toc.tv_nsec - 0.000000001 * tic.tv_nsec);
						break;
						case 'v':
						break;
						case SDLK_UP:
						break;
						case SDLK_DOWN:
						break;
						default:
						break;
					}
					break;
				case SDL_MOUSEBUTTONDOWN:
					if (event.button.button == SDL_BUTTON_LEFT) {
						mouseDrag = true;
						previousX = event.button.x;
						previousY = event.button.y;
					}
					break;
				case SDL_MOUSEMOTION:
					if (mouseDrag) {
						imgcpy(baseimage,screen);
						currentRect = getVarRect(previousX,previousY,event.motion.x,event.motion.y);
						for (int i=0;i<rectList.size();i++)
							drawRect(screen,rectList[i],yellow);
						drawRect(screen,currentRect,yellow);
						SDL_Flip(screen);
					}
					break;
				case SDL_MOUSEBUTTONUP:
					if (mouseDrag) {
						mouseDrag = false;
						rectList.push_back(currentRect);
						palette.resize(rectList.size() * K_PARAM);
						box.resize(rectList.size() * 8);
						region.resize(rectList.size());
						Color* currentPalette = palette.data() + palette.size() - K_PARAM;
						std::vector<Coordinate> *currentRegion = region.data() + region.size() - 1;
						float* currentBox = box.data() + box.size() - 8;
						getKMeans(labimage,&currentRect,currentPalette,K_PARAM);
#if FOREGROUND
						targetList.push_back(findTargetColor(labimage,&currentRect,currentPalette,K_PARAM));
#else
						targetList.push_back(findBackgroundColor(labimage,&currentRect,currentPalette,K_PARAM));
#endif
						floodFillRegion(labimage,currentRegion,&currentRect,currentPalette,K_PARAM,targetList.back());
						highlightRegion(screen,currentRegion,BOX_COLOR);
						*currentRegion = convexHull(*currentRegion);
						minBoundRect(currentRegion,currentBox,false);
						drawBox(screen,currentBox,red);
					} else {
						imgcpy(baseimage,screen);
						rectList.clear();
						region.clear();
						box.clear();
						palette.clear();
						targetList.clear();
						currentRect.w = 0;
						currentRect.h = 0;
					}
					SDL_Flip(screen);
					break;
				case SDL_QUIT:
					exit(0);
					break;
			}
		}
		usleep(1000);
	}

	fclose(target_point);
	delete[] baseimage.data;

}
