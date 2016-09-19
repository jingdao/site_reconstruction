#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <vector>
#include <math.h>
#include <SDL/SDL.h>
#define K_PARAM 5
#define EXPAND 2.0
#define OUTLIER_RATIO 0.05
#define BOX_COLOR blue

SDL_Surface *screen;
bool mouseDrag = false;
bool mouseRotate = false;
bool mouseTranslate = false;
int previousX,previousY;
int targetIndex = 1;

typedef struct {
	int width,height;
	unsigned char* data;
} Image;

typedef struct {
	unsigned char r,g,b;
} Color;

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

void rotateRect(float* box, float centerX, float centerY, int x1, int y1, int x2, int y2) {
	float ux = x1 - centerX;
	float uy = y1 - centerY;
	float mu = sqrt(ux*ux + uy*uy);
	ux /= mu;
	uy /= mu;
	float vx = x2 - centerX;
	float vy = y2 - centerY;
	float mv = sqrt(vx*vx + vy*vy);
	vx /= mv;
	vy /= mv;
	float cos_t = ux*vx + uy*vy;
	float sin_t = ux*vy - uy*vx;
	float scale = mv / mu;
	for (int i=0;i<4;i++) {
		float x = box[i*2] - centerX;
		float y = box[i*2+1] - centerY;
		box[i*2] = scale * (cos_t * x - sin_t * y) + centerX;
		box[i*2+1] = scale * (sin_t * x + cos_t * y) + centerY;
	}
}

void translateRect(float* box, float* centerX, float* centerY, int x1, int y1, int x2, int y2) {
	float dx = x2 - x1;
	float dy = y2 - y1;
	for (int i=0;i<4;i++) {
		box[i*2] += dx;
		box[i*2+1] += dy;
	}
	*centerX += dx;
	*centerY += dy;
}

void imgcpy(Image image,SDL_Surface* surf) {
	unsigned char* src = image.data, *dst = (unsigned char*)surf->pixels;
	for (int i=0;i<image.height;i++) {
		memcpy(dst,src,image.width*3);
		src += image.width*3;
		dst += surf->pitch;
	}
}

void writeBoxToFile(std::vector< std::vector<float> > *box, FILE* f) {
	fprintf(f,"%lu",box->size() * 4);
	for (size_t i=0;i<box->size();i++) {
		for (size_t j=0;j<8;j++) {
			fprintf(f," %f",box->at(i)[j]);
			if (j%2 == 1)
				fprintf(f," 1");
		}
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

int main(int argc, char* argv[]) {

	if (argc < 3) {
		printf("%s label_point.txt maxID [1.ppm ..]\n",argv[0]);
		return 1;
	}

	SDL_Init(SDL_INIT_VIDEO);
	SDL_WM_SetCaption("imlabel",NULL);
	Image baseimage = {0,0,NULL};
	Image modimage;
	if (!loadImageByIndex(targetIndex,&baseimage,true))
		return 1;
	screen = SDL_SetVideoMode(baseimage.width,baseimage.height,24,SDL_SWSURFACE);

	Color yellow = {255,255,0};
	Color red = {255,0,0};
	Color blue = {0,0,255};
	std::vector<SDL_Rect> rectList;
	SDL_Rect currentRect;
	std::vector< std::vector<float> > boxList;
	std::vector<float> centerX;
	std::vector<float> centerY;
	int boxID = 0;
	FILE* label_point = fopen(argv[1],"w");
	int maxID = atoi(argv[2]);

	imgcpy(baseimage,screen);
	SDL_Flip(screen);

	SDL_Event event;
	while (true) {
		while (SDL_PollEvent(&event)) {
			switch(event.type){
				case SDL_KEYDOWN:
					switch( event.key.keysym.sym ){
						case 'n':
						if (targetIndex >= maxID) {
							printf("maxID reached\n");
							break;
						}
						writeBoxToFile(&boxList,label_point);
						printf("%d.ppm: wrote %lu boxes\n",targetIndex,boxList.size());
						do {
							targetIndex++;
							char buffer[128];
							sprintf(buffer,"%d.key",targetIndex);
							FILE* key = fopen(buffer,"r");
							if (key) {
								fclose(key);
								break;
							} else 
								fprintf(label_point,"0\n");
						} while (targetIndex < maxID);
						if (loadImageByIndex(targetIndex,&baseimage,true)) {
							imgcpy(baseimage,screen);
							for (size_t i=0;i<boxList.size();i++)
								drawBox(screen,boxList[i].data(),yellow);
							SDL_Flip(screen);
						}
						break;
						default:
						break;
					}
					break;
				case SDL_MOUSEBUTTONDOWN:
					if (event.button.button == SDL_BUTTON_LEFT) {
						for (size_t j=0;j<boxList.size();j++) {
							float dx = centerX[j] - event.button.x;
							float dy = centerY[j] - event.button.y;
							if (fabs(dx)<20 && fabs(dy)<20) {
								mouseTranslate = true;
								boxID = j;
								break;
							}
							std::vector<float> box = boxList[j];
							for (size_t i=0;i<box.size()/2;i++) {
								dx = box[i*2] - event.button.x;
								dy = box[i*2+1] - event.button.y;
								if (fabs(dx)<20 && fabs(dy)<20) {
									mouseRotate = true;
									boxID = j;
									break;
								}
							}
							if (mouseRotate)
								break;
						}
						if (!mouseRotate && !mouseTranslate)
							mouseDrag = true;
						previousX = event.button.x;
						previousY = event.button.y;
					}
					break;
				case SDL_MOUSEMOTION:
					if (mouseDrag) {
						imgcpy(baseimage,screen);
						currentRect = getVarRect(previousX,previousY,event.motion.x,event.motion.y);
						for (size_t i=0;i<boxList.size();i++)
							drawBox(screen,boxList[i].data(),yellow);
						drawRect(screen,currentRect,yellow);
						SDL_Flip(screen);
					} else if (mouseRotate) {
						imgcpy(baseimage,screen);
						rotateRect(boxList[boxID].data(),centerX[boxID],centerY[boxID],previousX,previousY,event.motion.x,event.motion.y);
						previousX = event.button.x;
						previousY = event.button.y;
						for (size_t i=0;i<boxList.size();i++)
							drawBox(screen,boxList[i].data(),yellow);
						SDL_Flip(screen);
					} else if (mouseTranslate) {
						imgcpy(baseimage,screen);
						translateRect(boxList[boxID].data(),centerX.data()+boxID,centerY.data()+boxID,previousX,previousY,event.motion.x,event.motion.y);
						previousX = event.button.x;
						previousY = event.button.y;
						for (size_t i=0;i<boxList.size();i++)
							drawBox(screen,boxList[i].data(),yellow);
						SDL_Flip(screen);
					}
					break;
				case SDL_MOUSEBUTTONUP:
					if (mouseDrag) {
						mouseDrag = false;
						rectList.push_back(currentRect);
						std::vector<float> box;
                        box.push_back(currentRect.x);
                        box.push_back(currentRect.y);
                        box.push_back(currentRect.x + currentRect.w);
                        box.push_back(currentRect.y);
                        box.push_back(currentRect.x);
                        box.push_back(currentRect.y + currentRect.h);
                        box.push_back(currentRect.x + currentRect.w);
                        box.push_back(currentRect.y + currentRect.h);
						boxList.push_back(box);
						centerX.push_back(currentRect.x + currentRect.w / 2);
						centerY.push_back(currentRect.y + currentRect.h / 2);
					} else if (mouseRotate) {
						mouseRotate = false;
					} else if (mouseTranslate) {
						mouseTranslate = false;
					} else {
						imgcpy(baseimage,screen);
						rectList.clear();
						boxList.clear();
						centerX.clear();
						centerY.clear();
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

	fclose(label_point);
	delete[] baseimage.data;

}
