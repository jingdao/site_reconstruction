#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <SDL/SDL.h>
#define K_PARAM 3

SDL_Surface *screen;
bool mouseDrag = false;
int previousX,previousY;
int targetIndex = 1;
float cx,cy;

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

typedef struct {
	int size,capacity;
	void* data;
} List;

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
	for (int i=0;i<k;i++) {
		palette[i].r = rand() % 256;
		palette[i].g = rand() % 256;
		palette[i].b = rand() % 256;
	}
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
				int d=0;
				d += (colors[i].r - palette[j].r) * (colors[i].r - palette[j].r);
				d += (colors[i].g - palette[j].g) * (colors[i].g - palette[j].g);
				d += (colors[i].b - palette[j].b) * (colors[i].b - palette[j].b);
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
				palette[j].r = rand() % 256;
				palette[j].g = rand() % 256;
				palette[j].b = rand() % 256;
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
				int d=0;
				d += (c.r - palette[l].r) * (c.r - palette[l].r);
				d += (c.g - palette[l].g) * (c.g - palette[l].g);
				d += (c.b - palette[l].b) * (c.b - palette[l].b);
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
		for (int i=image->width*image->height;i>=0;i--) {
			unsigned char tmp = image->data[i*3+2];
			image->data[i*3+2] = image->data[i*3];
			image->data[i*3] = tmp;
		}
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

int main(int argc, char* argv[]) {

	if (argc < 2) {
		printf("./match_image target.pgm [1.pgm ..]\n");
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
	screen = SDL_SetVideoMode(baseimage.width,baseimage.height,24,SDL_SWSURFACE);

	cx = 0.5 * (baseimage.width-1);
	cy = 0.5 * (baseimage.height-1);

	Color yellow = {255,255,0};
	Color red = {255,0,0};
	Color blue = {0,0,255};
	Color palette[K_PARAM];
	List rectList = {0,8,NULL};
	rectList.data = calloc(8,sizeof(SDL_Rect));
	SDL_Rect* currentRect = (SDL_Rect*) rectList.data;
	unsigned char threshold;
	std::vector<Coordinate> region;
	std::vector<Coordinate> region2;
	std::vector<Coordinate> region3;
	Coordinate seed;

	imgcpy(baseimage,screen);
	SDL_Flip(screen);

	SDL_Event event;
	while (true) {
		while (SDL_PollEvent(&event)) {
			switch(event.type){
				case SDL_KEYDOWN:
					switch( event.key.keysym.sym ){
						case 'b':
						if (loadImageByIndex(targetIndex--,&baseimage,use_color)) {
							delete[] modimage.data;
							modimage = imclone(baseimage);
							binarize(baseimage,modimage,threshold);
							imgcpy(modimage,screen);
							SDL_Flip(screen);
						}
						break;
						case 'n':
						if (loadImageByIndex(targetIndex++,&baseimage,use_color)) {
							delete[] modimage.data;
							modimage = imclone(baseimage);
							binarize(baseimage,modimage,threshold);
							imgcpy(modimage,screen);
							SDL_Flip(screen);
						}
						break;
						case 'm':
						break;
						case 'v':
						break;
						case SDLK_UP:
						threshold++;
						delete[] modimage.data;
						modimage = imclone(baseimage);
						binarize(baseimage,modimage,threshold);
						imgcpy(modimage,screen);
						SDL_Flip(screen);
						break;
						case SDLK_DOWN:
						threshold--;
						delete[] modimage.data;
						modimage = imclone(baseimage);
						binarize(baseimage,modimage,threshold);
						imgcpy(modimage,screen);
						SDL_Flip(screen);
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
						rectList.size++;
					}
					break;
				case SDL_MOUSEMOTION:
					if (mouseDrag) {
						imgcpy(baseimage,screen);
						*currentRect = getVarRect(previousX,previousY,event.motion.x,event.motion.y);
						for (int i=0;i<rectList.size;i++)
							drawRect(screen,((SDL_Rect*)rectList.data)[i],yellow);
						SDL_Flip(screen);
					}
					break;
				case SDL_MOUSEBUTTONUP:
					if (mouseDrag) {
						mouseDrag = false;
						threshold = getKmeansThreshold(baseimage,currentRect);
						getKMeans(baseimage,currentRect,palette,K_PARAM);
//						seed.x = currentRect->x + currentRect->w/2;
//						seed.y = currentRect->y + currentRect->h/2;
//						sobel(baseimage,80);
//						regionGrow(baseimage,&region,seed,threshold);
//						mooreTracing(baseimage,&region,currentRect);
//						highlightRegion(screen,&region,red);
//						modimage = imclone(baseimage);
//						binarize(baseimage,modimage,threshold);
//						imgcpy(modimage,screen);
						colorRegion(baseimage,&region,currentRect,palette,K_PARAM,0);
						colorRegion(baseimage,&region2,currentRect,palette,K_PARAM,1);
						colorRegion(baseimage,&region3,currentRect,palette,K_PARAM,2);
						highlightRegion(screen,&region,red);
						highlightRegion(screen,&region2,yellow);
						highlightRegion(screen,&region3,blue);
						currentRect++;
					} else {
						imgcpy(baseimage,screen);
						currentRect = (SDL_Rect*) rectList.data;
						rectList.size = 0;
						currentRect->w = 0;
						currentRect->h = 0;
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

	delete[] baseimage.data;

}
