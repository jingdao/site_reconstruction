#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <SDL/SDL.h>
#define POINT_SIZE 5

SDL_Surface *screen;
unsigned char* palette;
bool useColor;
bool mouseDrag = false;
int previousX,previousY;

typedef struct {
	int width,height;
	unsigned char* data;
} Image;

typedef struct {
	int x,y,width,height;
} Rect;

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

void drawRect(Image image,Rect rect,Color color) {
	int top = rect.y < 0 ? 0 : rect.y ;
	int bottom = rect.y + rect.height >= image.height ? image.height-1 : rect.y + rect.height;
	int left = rect.x < 0 ? 0 : rect.x;
	int right = rect.x + rect.width >= image.width ? image.width-1 : rect.x + rect.width;
	//TOP
	int i=left;
	unsigned char* c = image.data + (top*image.width+left)*3;
	while (i++ < right) {
		*c++ = color.b;
		*c++ = color.g;
		*c++ = color.r;
	}	
	//BOTTOM
	i=left;
	c = image.data + (bottom*image.width+left)*3;
	while (i++ < right) {
		*c++ = color.b;
		*c++ = color.g;
		*c++ = color.r;
	}	
	//LEFT
	i=top;
	c = image.data + (top*image.width+left)*3;
	while (i++ < bottom) {
		c[0] = color.b;
		c[1] = color.g;
		c[2] = color.r;
		c += image.width * 3;
	}	
	//RIGHT
	i=top;
	c = image.data + (top*image.width+right)*3;
	while (i++ < bottom) {
		c[0] = color.b;
		c[1] = color.g;
		c[2] = color.r;
		c += image.width * 3;
	}
}

Rect getVarRect(int x1,int y1,int x2,int y2) {
	int left = x1 < x2 ? x1 : x2;
	int right = x1 > x2 ? x1 : x2;
	int top = y1 < y2 ? y1 : y2;
	int bottom = y1 > y2 ? y1 : y2;
	Rect r = {left,top,right-left,bottom-top};
	return r;
}

void drawKeyPoint(Image image,int x,int y,Color color) {
	Rect r = {x-POINT_SIZE,y-POINT_SIZE,POINT_SIZE*2,POINT_SIZE*2};
	drawRect(image,r,color);
}

void findCorrespondence(List matches,Rect rect) {
	Match* m = (Match*)matches.data;
	int numMatches = 0;
	for (int i=0;i<matches.size;i++) {
		if (m[i].x1 > rect.x && m[i].x1 < rect.x + rect.width &&
			m[i].y1 > rect.y && m[i].y1 < rect.y + rect.height) {
			printf("%d %f %f\n",numMatches++,m[i].x2,m[i].y2);
		}
	}
}

int main(int argc, char* argv[]) {

	if (argc < 3) {
		printf("./mark_image in.ppm/in.pgm key.match\n");
		return 1;
	}

	SDL_Init(SDL_INIT_VIDEO);
	SDL_WM_SetCaption("mark_image",NULL);
	char buffer[128];
	FILE* ppm = fopen(argv[1],"r");
	if (!ppm) {
		printf("%s not found\n",argv[1]);
		return 1;
	}
	fgets(buffer,128,ppm); //P5 or P6
	useColor = strncmp(buffer,"P6",2) == 0;
	fgets(buffer,128,ppm);
	char *c = buffer;
	int width = strtol(c,&c,10);
	int height = strtol(c,&c,10);
	fgets(buffer,128,ppm); //255
	screen = SDL_SetVideoMode(width,height,24,SDL_SWSURFACE);
	palette = malloc(width*height*3);
	if (useColor) {
		fread(palette,1,width*height*3,ppm);
		//swap RGB order
		for (int i=0;i<width*height*3;i+=3) {
			unsigned char tmp = palette[i];
			palette[i] = palette[i+2];
			palette[i+2] = tmp;
		}
	} else {
		fread(palette,1,width*height,ppm);
		//populate gray level
		for (int i=width*height;i>=0;i--) {
			palette[i*3+2] = palette[i*3+1] = palette[i*3] = palette[i];
		}
	}
	memcpy(screen->pixels,palette,width*height*3);
	SDL_Flip(screen);
	fclose(ppm);
	Image image = {width,height,palette};
	Image screenImage = {width,height,screen->pixels};
	Color red = {255,0,0};
	Color yellow = {255,255,0};

	FILE* key_match = fopen(argv[2],"r");
	List matches = {0,8,NULL};
	Rect currentRect;
	matches.data = malloc(8*sizeof(Match));
	while (fgets(buffer,128,key_match)) {
		Match m;
		int id1,id2;
		if (sscanf(buffer,"%d %f %f %d %f %f",&id1,&m.x1,&m.y1,&id2,&m.x2,&m.y2)==6) {
			if (matches.size == matches.capacity) {
				matches.data = realloc(matches.data,matches.capacity*2*sizeof(Match));
				matches.capacity *= 2;
			}
			((Match*)matches.data)[matches.size++] = m;
		}
	}
	printf("Loaded %d matches from %s\n",matches.size,argv[2]);
	fclose(key_match);

	SDL_Event event;
	while (true) {
		while (SDL_PollEvent(&event)) {
			switch(event.type){
				case SDL_MOUSEBUTTONDOWN:
					if (event.button.button == SDL_BUTTON_LEFT) {
						mouseDrag = true;
						previousX = event.button.x;
						previousY = event.button.y;
					} else {
						drawKeyPoint(image,event.button.x,event.button.y,red);
						memcpy(screen->pixels,palette,width*height*3);
						SDL_Flip(screen);
					}
					break;
				case SDL_MOUSEMOTION:
					if (mouseDrag) {
						memcpy(screen->pixels,palette,width*height*3);
						currentRect = getVarRect(previousX,previousY,event.motion.x,event.motion.y);
						drawRect(screenImage,currentRect,yellow);
						SDL_Flip(screen);
					}
					break;
				case SDL_MOUSEBUTTONUP:
					if (mouseDrag) {
						mouseDrag = false;
						findCorrespondence(matches,currentRect);
//						memcpy(screen->pixels,palette,width*height*3);
//						SDL_Flip(screen);
					}
					break;
				case SDL_QUIT:
					exit(0);
					break;
			}
		}
		usleep(1000);
	}

	free(palette);

}
