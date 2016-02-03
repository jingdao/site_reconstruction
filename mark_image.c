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

void drawKeyPoint(Image image,int x,int y,Color color) {
	Rect r = {x-POINT_SIZE,y-POINT_SIZE,POINT_SIZE*2,POINT_SIZE*2};
	drawRect(image,r,color);
}

void imgcpy(Image image,SDL_Surface* surf) {
	unsigned char* src = image.data, *dst = surf->pixels;
	for (int i=0;i<image.height;i++) {
		memcpy(dst,src,image.width*3);
		src += image.width*3;
		dst += surf->pitch;
	}
}

int main(int argc, char* argv[]) {

	if (argc < 2) {
		printf("./mark_image in.ppm/in.pgm\n");
		return 1;
	}

	SDL_Init(SDL_INIT_VIDEO);
	SDL_WM_SetCaption(argv[1],NULL);
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
	fclose(ppm);
	Image image = {width,height,palette};
	Image screenImage = {width,height,screen->pixels};
	Color red = {255,0,0};
	if (useColor)
		strcpy(strstr(argv[1],".ppm"),".rpt");
	else
		strcpy(strstr(argv[1],".pgm"),".rpt");
	FILE* point_file = fopen(argv[1],"w");
	if (!point_file) {
		printf("%s not found\n",argv[1]);
		return 1;
	}
	double cx = 0.5 * (width-1);
	double cy = 0.5 * (height-1);
	imgcpy(image,screen);
	SDL_Flip(screen);

	SDL_Event event;
	while (true) {
		while (SDL_PollEvent(&event)) {
			switch(event.type){
				case SDL_MOUSEBUTTONDOWN:
					if (event.button.button == SDL_BUTTON_LEFT) {
						fprintf(point_file,"%f %f\n",(event.button.x-cx)/width,(event.button.y-cy)/height);
						drawKeyPoint(image,event.button.x,event.button.y,red);
						imgcpy(image,screen);
						SDL_Flip(screen);
					}
					break;
				case SDL_MOUSEMOTION:
					break;
				case SDL_MOUSEBUTTONUP:
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
