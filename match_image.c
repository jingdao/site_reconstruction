#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <SDL/SDL.h>
#define POINT_SIZE 5

SDL_Surface *screen, *leftscreen, *rightscreen;
bool mouseDrag = false;
int previousX,previousY;
SDL_Rect leftrect,rightrect;
int targetIndex = 0;
int width,height;
float cx,cy;
FILE* target_point;

typedef struct {
	int width,height;
	unsigned char* data;
} Image;

typedef struct {
	unsigned char r,g,b;
} Color;

typedef struct {
	float x1,y1,x2,y2;
} Match;

typedef struct {
	int size,capacity;
	void* data;
} List;

int compare_float (const void * a, const void * b) {
	float c = *(float*) a - *(int*) b;
	if (c > 0) return 1;
	else if (c < 0) return -1;
	else return 0;
}

void drawRect(SDL_Surface *surf,SDL_Rect rect,Color color) {
	int top = rect.y < 0 ? 0 : rect.y ;
	int bottom = rect.y + rect.h >= surf->h ? surf->h-1 : rect.y + rect.h;
	int left = rect.x < 0 ? 0 : rect.x;
	int right = rect.x + rect.w >= surf->w ? surf->w-1 : rect.x + rect.w;
	//TOP
	int i=left;
	unsigned char* c = surf->pixels + top * surf->pitch + left*3;
	while (i++ < right) {
		*c++ = color.b;
		*c++ = color.g;
		*c++ = color.r;
	}	
	//BOTTOM
	i=left;
	c = surf->pixels + bottom * surf->pitch + left*3;
	while (i++ < right) {
		*c++ = color.b;
		*c++ = color.g;
		*c++ = color.r;
	}	
	//LEFT
	i=top;
	c = surf->pixels + top * surf->pitch + left*3;
	while (i++ < bottom) {
		c[0] = color.b;
		c[1] = color.g;
		c[2] = color.r;
		c += surf->pitch;
	}	
	//RIGHT
	i=top;
	c = surf->pixels + top * surf->pitch + right*3;
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
	SDL_Rect r = {left,top,right-left,bottom-top};
	return r;
}

void drawKeyPoint(SDL_Surface *surf,int x,int y,Color color) {
	SDL_Rect r = {x-POINT_SIZE,y-POINT_SIZE,POINT_SIZE*2,POINT_SIZE*2};
	drawRect(surf,r,color);
}

void imgcpy(Image image,SDL_Surface* surf) {
	unsigned char* src = image.data, *dst = surf->pixels;
	for (int i=0;i<image.height;i++) {
		memcpy(dst,src,image.width*3);
		src += image.width*3;
		dst += surf->pitch;
	}
}

void showRightImage(int index,List* matches,Image rightimage) {
	char buffer[128];
	char filename[128];
	sprintf(filename,"%d.match",index);
	FILE* key_match = fopen(filename,"r");
	if (!key_match) {
		printf("Cannot open %s\n",filename);
		return;
	}
	matches->size=0;
	while (fgets(buffer,128,key_match)) {
		Match m;
		int id1,id2;
		if (sscanf(buffer,"%d %f %f %d %f %f",&id1,&m.x1,&m.y1,&id2,&m.x2,&m.y2)==6) {
			if (matches->size == matches->capacity) {
				matches->data = realloc(matches->data,matches->capacity*2*sizeof(Match));
				matches->capacity *= 2;
			}
			((Match*)matches->data)[matches->size++] = m;
		}
	}
	printf("Loaded %d matches from %s\n",matches->size,filename);
	fclose(key_match);
	sprintf(filename,"%d.pgm",index);
	FILE* ppm = fopen(filename,"r");
	fgets(buffer,128,ppm); //P5 or P6
	fgets(buffer,128,ppm);
	fgets(buffer,128,ppm); //255
	fread(rightimage.data,1,rightimage.width*rightimage.height,ppm);
	//populate gray level
	for (int i=rightimage.width*rightimage.height;i>=0;i--) {
		rightimage.data[i*3+2] = rightimage.data[i*3+1] = rightimage.data[i*3] = rightimage.data[i];
	}
	imgcpy(rightimage,rightscreen);
	fclose(ppm);
}

void findCorrespondence(List matches,SDL_Rect rect) {
	if (rect.w <= 0 || rect.h <= 0)
		return;
	Match* m = (Match*)matches.data;
	Color red = {255,0,0};
	List xl = {0,8,NULL};
	List yl = {0,8,NULL};
	xl.data = malloc(8*sizeof(float));
	yl.data = malloc(8*sizeof(float));
	for (int i=0;i<matches.size;i++) {
		if (m[i].x1 > rect.x && m[i].x1 < rect.x + rect.w &&
			m[i].y1 > rect.y && m[i].y1 < rect.y + rect.h) {
			drawKeyPoint(rightscreen,m[i].x2,m[i].y2,red);
			if (xl.size == xl.capacity) {
				xl.data = realloc(xl.data,xl.capacity*2*sizeof(float));
				xl.capacity *= 2;
				yl.data = realloc(yl.data,yl.capacity*2*sizeof(float));
				yl.capacity *= 2;
			}
			((float*)xl.data)[xl.size++] = m[i].x2;
			((float*)yl.data)[yl.size++] = m[i].y2;
		}
	}
	qsort(xl.data,xl.size,sizeof(float),compare_float);
	qsort(yl.data,yl.size,sizeof(float),compare_float);
	float x_med = (((float*)xl.data)[xl.size/2] - cx) / width;
	float y_med = (((float*)yl.data)[yl.size/2] - cy) / height;
	fprintf(target_point,"%f %f\n",x_med,y_med);
	free(xl.data);
	free(yl.data);
}

int main(int argc, char* argv[]) {

	if (argc < 2) {
		printf("./match_image target_point.txt target.pgm\n");
		return 1;
	}

	SDL_Init(SDL_INIT_VIDEO);
	SDL_WM_SetCaption("match_image",NULL);
	char buffer[128];
	target_point = fopen(argv[1],"w");
	FILE* ppm = fopen(argv[2],"r");
	if (!ppm) {
		printf("%s not found\n",argv[2]);
		return 1;
	}
	fgets(buffer,128,ppm); //P5 or P6
	fgets(buffer,128,ppm);
	char *c = buffer;
	width = strtol(c,&c,10);
	height = strtol(c,&c,10);
	fgets(buffer,128,ppm); //255
	screen = SDL_SetVideoMode(width*2,height,24,SDL_SWSURFACE);
	leftscreen = SDL_CreateRGBSurface(0,width,height,24,0xFF0000,0xFF00,0xFF,0);
	rightscreen = SDL_CreateRGBSurface(0,width,height,24,0xFF0000,0xFF00,0xFF,0);
	Image leftimage = {width,height,NULL};
	Image rightimage = {width,height,NULL};
	leftimage.data = malloc(width*height*3);
	rightimage.data = malloc(width*height*3);
	leftrect.x = 0; leftrect.y = 0; leftrect.w = width; leftrect.h = height;
	rightrect.x = width; rightrect.y = 0; rightrect.w = width; rightrect.h = height;
	cx = 0.5 * (width-1);
	cy = 0.5 * (height-1);

	fread(leftimage.data,1,width*height,ppm);
	//populate gray level
	for (int i=width*height;i>=0;i--) {
		leftimage.data[i*3+2] = leftimage.data[i*3+1] = leftimage.data[i*3] = leftimage.data[i];
	}
	fclose(ppm);
	Color yellow = {255,255,0};
	List matches = {0,8,NULL};
	SDL_Rect currentRect = {0,0,0,0};
	matches.data = malloc(8*sizeof(Match));

	imgcpy(leftimage,leftscreen);
	SDL_BlitSurface(leftscreen,NULL,screen,&leftrect);
	showRightImage(targetIndex++,&matches,rightimage);
	SDL_BlitSurface(rightscreen,NULL,screen,&rightrect);
	SDL_Flip(screen);

	SDL_Event event;
	while (true) {
		while (SDL_PollEvent(&event)) {
			switch(event.type){
				case SDL_KEYDOWN:
					switch( event.key.keysym.sym ){
						case 'b':
						showRightImage(targetIndex--,&matches,rightimage);
						findCorrespondence(matches,currentRect);
						SDL_BlitSurface(rightscreen,NULL,screen,&rightrect);
						SDL_Flip(screen);
						break;
						case 'n':
						showRightImage(targetIndex++,&matches,rightimage);
						findCorrespondence(matches,currentRect);
						SDL_BlitSurface(rightscreen,NULL,screen,&rightrect);
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
					}
					break;
				case SDL_MOUSEMOTION:
					if (mouseDrag) {
						imgcpy(leftimage,leftscreen);
						currentRect = getVarRect(previousX,previousY,event.motion.x,event.motion.y);
						drawRect(leftscreen,currentRect,yellow);
						SDL_BlitSurface(leftscreen,NULL,screen,&leftrect);
						SDL_Flip(screen);
					}
					break;
				case SDL_MOUSEBUTTONUP:
					imgcpy(rightimage,rightscreen);
					SDL_BlitSurface(rightscreen,NULL,screen,&rightrect);
					if (mouseDrag) {
						mouseDrag = false;
						findCorrespondence(matches,currentRect);
						SDL_BlitSurface(rightscreen,NULL,screen,&rightrect);
					} else {
						imgcpy(leftimage,leftscreen);
						SDL_BlitSurface(leftscreen,NULL,screen,&leftrect);
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
	free(leftimage.data);
	free(rightimage.data);

}
