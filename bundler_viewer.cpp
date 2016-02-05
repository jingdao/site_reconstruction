#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <vector>
#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>

struct Point {
	float x,y,z;
	unsigned char r,g,b;
};

struct CamModel {
	Point center,ul,ur,bl,br;
};

double cameraX=0,cameraY=-5,cameraZ=5;
double centerX=0,centerY=0,centerZ=0;
double upX=0,upY=0,upZ=1;
int screenWidth = 1100, screenHeight = 750;
float cameraSize = 0.1;
std::vector<CamModel> camera_location;
std::vector<Point> pointcloud;
int mouseIndex = 0;
int previousX,previousY;
double scrollSpeed = 1.1;
float Twc[] = {0,-0.25,-0.18};

void invertTransformation(float* R,float* T) {
	float Tcw[3] = {T[0],T[1],T[2]};
	float tmp;
	tmp = R[1]; R[1] = R[3]; R[3] = tmp;
	tmp = R[2]; R[2] = R[6]; R[6] = tmp;
	tmp = R[5]; R[5] = R[7]; R[7] = tmp;
	T[0] = - (R[0] * Tcw[0] + R[1] * Tcw[1] + R[2] * Tcw[2]);
	T[1] = - (R[3] * Tcw[0] + R[4] * Tcw[1] + R[5] * Tcw[2]);
	T[2] = - (R[6] * Tcw[0] + R[7] * Tcw[1] + R[8] * Tcw[2]);
}

Point transformPoint(Point p,float* R,float* T) {
	Point res;
	res.x = R[0] * p.x  + R[1] * p.y + R[2] * p.z + T[0];
	res.y = R[3] * p.x  + R[4] * p.y + R[5] * p.z + T[1];
	res.z = R[6] * p.x  + R[7] * p.y + R[8] * p.z + T[2];
	return res;
}

void drawLine(Point p1,Point p2) {
	glBegin(GL_LINES);
	glVertex3d(p1.x,p1.y,p1.z);
	glVertex3d(p2.x,p2.y,p2.z);
	glEnd();
}

void draw() {
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glPushMatrix();
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(cameraX,cameraY,cameraZ,centerX,centerY,centerZ,upX,upY,upZ);

	glLineWidth(1.0);
	for (size_t n = 0; n < camera_location.size(); n++){
		glColor3ub(0,255,0);
		drawLine(camera_location[n].center,camera_location[n].ul);
		drawLine(camera_location[n].center,camera_location[n].ur);
		drawLine(camera_location[n].center,camera_location[n].bl);
		drawLine(camera_location[n].center,camera_location[n].br);
		drawLine(camera_location[n].ul,camera_location[n].bl);
		drawLine(camera_location[n].ur,camera_location[n].br);
		drawLine(camera_location[n].ul,camera_location[n].ur);
		drawLine(camera_location[n].bl,camera_location[n].br);
	}

	glPointSize(2.0);
	glBegin(GL_POINTS);
	for (size_t n = 0; n < pointcloud.size(); n++){
		glColor3ub(pointcloud[n].r,pointcloud[n].g,pointcloud[n].b);
		glVertex3d(pointcloud[n].x,pointcloud[n].y,pointcloud[n].z);
	}
	glEnd();

	glFlush();
	SDL_GL_SwapBuffers();

	glPopMatrix();
	glPopAttrib();
}

int main(int argc,char* argv[]) {
	if (argc < 2) {
		printf("./bundler_viewer bundle.out\n");
		return 1;
	}

	FILE* f = fopen(argv[1],"r");
	if (!f) {
		printf("Cannot open %s\n",argv[1]);
		return 1;
	}
	int num_cameras,num_points;
	char buffer[512];
	float R[9];
	float T[3];
	float focal_length,k1,k2,t;
	Point center = {0,0,0};
	Point ul = {-cameraSize,cameraSize,-cameraSize};
	Point ur = {cameraSize,cameraSize,-cameraSize};
	Point bl = {-cameraSize,-cameraSize,-cameraSize};
	Point br = {cameraSize,-cameraSize,-cameraSize};
	fgets(buffer,512,f);
	fgets(buffer,512,f);
	sscanf(buffer,"%d %d",&num_cameras,&num_points);
	for (int i=0;i<num_cameras;i++) {
		fgets(buffer,512,f);
		sscanf(buffer,"%f %f %f\n",&focal_length,&k1,&k2);
		fgets(buffer,512,f);
		sscanf(buffer,"%f %f %f\n",R,R+1,R+2);
		fgets(buffer,512,f);
		sscanf(buffer,"%f %f %f\n",R+3,R+4,R+5);
		fgets(buffer,512,f);
		sscanf(buffer,"%f %f %f\n",R+6,R+7,R+8);
		fgets(buffer,512,f);
		sscanf(buffer,"%f %f %f\n",T,T+1,T+2);
		if (focal_length > 0) {
			invertTransformation(R,T);
			CamModel cam;
			cam.center = transformPoint(center,R,T);
			cam.ul = transformPoint(ul,R,T);
			cam.ur = transformPoint(ur,R,T);
			cam.bl = transformPoint(bl,R,T);
			cam.br = transformPoint(br,R,T);
			camera_location.push_back(cam);
		}
	}
	for (int i=0;i<num_points;i++) {
		fgets(buffer,512,f);
		Point p;
		sscanf(buffer,"%f %f %f",&p.x,&p.y,&p.z);
		fgets(buffer,512,f);
		sscanf(buffer,"%hhu %hhu %hhu",&p.r,&p.g,&p.b);
		pointcloud.push_back(p);
		fgets(buffer,512,f);
	}
	printf("Parsed %lu cameras %lu points\n",camera_location.size(),pointcloud.size());
	fclose(f);

	SDL_Init(SDL_INIT_VIDEO);
	SDL_WM_SetCaption("Viz Cam", NULL);
	SDL_SetVideoMode(screenWidth,screenHeight, 32, SDL_OPENGL);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(70,(double)screenWidth/screenHeight,1,1000);

	int interval = 10000;
	SDL_Event event;
	while (SDL_PollEvent(&event)); //clear event buffer
	draw();
	while (true) {
		while (SDL_PollEvent(&event)) {
			switch(event.type){
				case SDL_KEYDOWN:
					switch( event.key.keysym.sym ){
						case SDLK_LEFT:
						cameraX -= 1;
						break;
						case SDLK_RIGHT:
						cameraX += 1;
						break;
						case SDLK_UP:
						cameraZ += 1;
						break;
						case SDLK_DOWN:
						cameraZ -= 1;
						break;
						default:
						break;
					}
					draw();
					break;
				case SDL_MOUSEBUTTONDOWN:
					if (event.button.button == SDL_BUTTON_WHEELUP) {
						cameraX /= scrollSpeed;
						cameraY /= scrollSpeed;
						cameraZ /= scrollSpeed;
						draw();
					} else if (event.button.button == SDL_BUTTON_WHEELDOWN) {
						cameraX *= scrollSpeed;
						cameraY *= scrollSpeed;
						cameraZ *= scrollSpeed;
						draw();
					} else {
						mouseIndex = event.button.button == SDL_BUTTON_LEFT ? 1 : 2;
						previousX = event.button.x;
						previousY = event.button.y;
					}
					break;
				case SDL_MOUSEBUTTONUP:
					mouseIndex = 0;
					break;
				case SDL_MOUSEMOTION:
					if (mouseIndex == 1) {
						double rho = sqrt(cameraX*cameraX+cameraY*cameraY);
						double xstep = cameraY / rho;
						double ystep = -cameraX / rho;
						cameraX += 0.05 * (event.motion.x-previousX) * xstep;
						cameraY += 0.05 * (event.motion.x-previousX) * ystep;
						cameraZ += 0.05 * (event.motion.y-previousY);
						previousX = event.motion.x;
						previousY = event.motion.y;
						draw();
					}
					break;
				case SDL_QUIT:
					exit(0);
					break;
			}
		}
		usleep(interval);
	}

//	atexit(SQL_Quit);

	return 0;
}
