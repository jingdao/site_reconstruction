#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <vector>
#define DEBUG 1

struct PCD {
	int numPoints;
	float* float_data;
	unsigned char* color_data;
};
enum PCD_data_storage {
	ASCII,
	BINARY,
	NONE
};
struct Box {
	float minX,minY,minZ,maxX,maxY,maxZ;
};
struct Point {
	float x,y,z;
};
struct CamModel {
	Point center,ul,ur,bl,br;
};


double cameraX=30,cameraY=50,cameraZ=40;
double centerX=0,centerY=0,centerZ=0;
double upX=0,upY=0,upZ=1;
int mouseIndex = 0;
int previousX,previousY;
double scrollSpeed = 1.01;
float cameraSize = 0.5;
PCD* cloud;
PCD* object;
std::vector<Point> object_location;
int location_index = 0;
std::vector<CamModel> camera_location;
std::vector<float> location_error;
std::vector<int> match_target;
std::vector<Point> scan_line;

float getDistance(float x1,float y1,float x2,float y2) {
	float S = 0;
	S += (x1-x2) * (x1-x2);
	S += (y1-y2) * (y1-y2);
	return sqrt(S);
}

float getAngle(float x1,float y1,float x2,float y2,float x3,float y3) {
	float a = getDistance(x1,y1,x2,y2);
	float b = getDistance(x2,y2,x3,y3);
	float c = getDistance(x1,y1,x3,y3);
	return acos( (a*a+b*b-c*c) / (2*a*b));
}

Point getLocation(Point p1,Point p2,float theta,float distance) {
	float phi = atan2(p2.y-p1.y,p2.x-p1.x);
	Point q;
	q.x = p1.x + distance * cos(phi + theta);
	q.y = p1.y + distance * sin(phi + theta);
	q.z = 5;
	return q;
}

void rotationFromAngle(float* M, float rx,float ry, float rz) {
	float r1 = rx;
	float r2 = ry;
	float r3 = rz;
	M[0] = cos(r2)*cos(r3);
	M[1] = - cos(r1) * sin(r3) + cos(r3) * sin(r1) * sin(r2);
	M[2] = sin(r1) * sin(r3) + cos(r1) * cos(r3) * sin(r2);
	M[3] = cos(r2) * sin(r3);
	M[4] = cos(r1) * cos(r3) + sin(r1) * sin(r2) * sin(r3);
	M[5] = - cos(r3) * sin(r1) + cos(r1) * sin(r2) * sin(r3);
	M[6] = -sin(r2);
	M[7] = cos(r2) * sin(r1);
	M[8] = cos(r1) * cos(r2);
}

Point transformPoint(Point p,float* R,float* T) {
	Point res;
	res.x = R[0] * p.x  + R[1] * p.y + R[2] * p.z + T[0];
	res.y = R[3] * p.x  + R[4] * p.y + R[5] * p.z + T[1];
	res.z = R[6] * p.x  + R[7] * p.y + R[8] * p.z + T[2];
	return res;
}

bool isValid(int index) {
	return location_error[index]>=0 && location_error[index]<=10 && match_target[index] > 0;
}

PCD* NewPCD(const char* fileName, Box* box) {
	PCD* pcd = new PCD();
	PCD_data_storage data_storage = NONE;
	FILE* f = fopen(fileName, "r");
	if (!f) {
		printf("File not found: %s\n", fileName);
		return NULL;
	}
	char buf[256];
	int pointsParsed = 0;
	while (fgets(buf, 256, f)) {
		if (sscanf(buf, "POINTS %d", &pcd->numPoints) == 1) {
			pcd->float_data = (float*)malloc(4 * pcd->numPoints * sizeof(float));
			pcd->color_data = (unsigned char*)malloc(3 * pcd->numPoints);
		} else if (strncmp(buf,"DATA ascii",10)==0) {
			data_storage = ASCII;
		} else if (strncmp(buf,"DATA binary_compressed",23)==0) {
			data_storage = BINARY;
			fread(pcd->float_data,sizeof(float),pcd->numPoints*4,f);
			break;
		}
		else if (data_storage == ASCII) {
			float x,y,z,c;
			int numParsed = sscanf(buf,"%f %f %f %f",&x,&y,&z,&c);
			if (numParsed >= 3 && !(box && x>box->minX && x<box->maxX && y>box->minY && y<box->maxY && z>box->minZ && z<box->maxZ)) {
				pcd->float_data[pointsParsed * 4] = x;
				pcd->float_data[pointsParsed * 4 + 1] = y;
				pcd->float_data[pointsParsed * 4 + 2] = z;
				if (numParsed >= 4) {
					int rgb = (int) c;
					pcd->color_data[pointsParsed*3] = (rgb >> 16) & 0xFF;
					pcd->color_data[pointsParsed*3 + 1] = (rgb >> 8) & 0xFF;
					pcd->color_data[pointsParsed*3 + 2] = rgb & 0xFF;
				}
				pointsParsed++;
			}
		}
	}
	fclose(f);
	return pcd;
}

Box getBoundingBox(PCD* pcd) {
	Box b;
	b.minX = b.maxX = pcd->float_data[0];
	b.minY = b.maxY = pcd->float_data[1];
	b.minZ = b.maxZ = pcd->float_data[2];
	for (int i=0;i<pcd->numPoints;i++) {
		if (pcd->float_data[i*4] < b.minX) b.minX = pcd->float_data[i*4];
		else if (pcd->float_data[i*4] > b.maxX) b.maxX = pcd->float_data[i*4];
		else if (pcd->float_data[i*4+1] < b.minY) b.minY = pcd->float_data[i*4+1];
		else if (pcd->float_data[i*4+1] > b.maxY) b.maxY = pcd->float_data[i*4+1];
		else if (pcd->float_data[i*4+2] < b.minZ) b.minZ = pcd->float_data[i*4+2];
		else if (pcd->float_data[i*4+2] > b.maxZ) b.maxZ = pcd->float_data[i*4+2];
	}
	return b;
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

//	glMatrixMode(GL_PROJECTION);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	gluLookAt(cameraX,cameraY,cameraZ,centerX,centerY,centerZ,upX,upY,upZ);

	glPointSize(1.0);
	glBegin(GL_POINTS);
	for (int n = 0; n < cloud->numPoints; n++){
		glColor3ub(cloud->color_data[n*3],cloud->color_data[n*3+1],cloud->color_data[n*3+2]);
		glVertex3d(cloud->float_data[n*4],cloud->float_data[n*4+1],cloud->float_data[n*4+2]);
	}
	glEnd();

	glPointSize(12.0);
	glBegin(GL_POINTS);
	if (location_index == 0)
		glColor3ub(255,255,0);
	else if (location_error[location_index-1] < 0 || match_target[location_index-1] <=2)
		glColor3ub(255,0,0);
	else {
		double intensity = 255.0 / ((int)location_error[location_index-1] + 1);
		glColor3ub(0,(unsigned char)intensity,0);
	}

	for (int n = 0; n < object->numPoints; n++){
		glVertex3d(object->float_data[n*4],object->float_data[n*4+1],object->float_data[n*4+2]);
	}
	glEnd();

	glLineWidth(1.0);
	glColor3ub(0,255,0);
	if (location_index == 0) {
		for (size_t n = 0; n < camera_location.size(); n++){
			drawLine(camera_location[n].center,camera_location[n].ul);
			drawLine(camera_location[n].center,camera_location[n].ur);
			drawLine(camera_location[n].center,camera_location[n].bl);
			drawLine(camera_location[n].center,camera_location[n].br);
			drawLine(camera_location[n].ul,camera_location[n].bl);
			drawLine(camera_location[n].ur,camera_location[n].br);
			drawLine(camera_location[n].ul,camera_location[n].ur);
			drawLine(camera_location[n].bl,camera_location[n].br);
		}
	} else {
		drawLine(camera_location[location_index-1].center,camera_location[location_index-1].ul);
		drawLine(camera_location[location_index-1].center,camera_location[location_index-1].ur);
		drawLine(camera_location[location_index-1].center,camera_location[location_index-1].bl);
		drawLine(camera_location[location_index-1].center,camera_location[location_index-1].br);
		drawLine(camera_location[location_index-1].ul,camera_location[location_index-1].bl);
		drawLine(camera_location[location_index-1].ur,camera_location[location_index-1].br);
		drawLine(camera_location[location_index-1].ul,camera_location[location_index-1].ur);
		drawLine(camera_location[location_index-1].bl,camera_location[location_index-1].br);
#if DEBUG
		if (isValid(location_index - 1)) {
			glColor3ub(255,0,0);
			drawLine(scan_line[2*(location_index-1)],scan_line[2*(location_index-1)+1]);
		}
#endif
	}


	glFlush();
	SDL_GL_SwapBuffers();

	glPopMatrix();
	glPopAttrib();
}

void translatePCD(PCD* pcd,float dx,float dy,float dz) {
	for (int i=0;i<pcd->numPoints;i++) {
		pcd->float_data[i*4] += dx;
		pcd->float_data[i*4+1] += dy;
		pcd->float_data[i*4+2] += dz;
	}
}

int main(int argc,char* argv[]) {
	if (argc < 4) {
		printf("./site_viewer scenario.pcd object.pcd camera_location.txt\n");
		return 1;
	}

	char buffer[128];
	object = NewPCD(argv[2],NULL);
	if (!object)
		return 1;
	Box object_box = getBoundingBox(object);
	Point previous_location = {
		(object_box.maxX + object_box.minX)/2,
		(object_box.maxY + object_box.minY)/2,
		(object_box.maxZ + object_box.minZ)/2,
	};
	Point current_location;
	printf("Loaded %s (%f %f %f %f %f %f)\n",argv[2],object_box.minX,object_box.maxX,object_box.minY,object_box.maxY,object_box.minZ,object_box.maxZ);
	cloud = NewPCD(argv[1],&object_box);
	if (!cloud)
		return 1;

	Point center = {0,0,0};
	Point ul = {-cameraSize,-cameraSize,cameraSize};
	Point ur = {cameraSize,-cameraSize,cameraSize};
	Point bl = {-cameraSize,cameraSize,cameraSize};
	Point br = {cameraSize,cameraSize,cameraSize};
	FILE* ref_camera_loc = fopen(argv[3],"r");
	if (!ref_camera_loc) {
		printf("Cannot open %s\n",argv[3]);
		return 1;
	}
	float T[3];
	float R[9];
	float rx,ry,rz;
	Point viewVec;
	float err;
	int numMatch;
	while (fgets(buffer,128,ref_camera_loc)) {
		if (sscanf(buffer,"%f %f %f %f %f %f %f %f %f %f %d",T,T+1,T+2,&rx,&ry,&rz,&viewVec.x,&viewVec.y,&viewVec.z,&err,&numMatch)==11) {
			rotationFromAngle(R,rx,ry,rz);
			CamModel cam;
			cam.center = transformPoint(center,R,T);
			cam.ul = transformPoint(ul,R,T);
			cam.ur = transformPoint(ur,R,T);
			cam.bl = transformPoint(bl,R,T);
			cam.br = transformPoint(br,R,T);
			camera_location.push_back(cam);
			float scale = (previous_location.z - T[2]) / viewVec.z;
			Point projection = {T[0]+scale*viewVec.x,T[1]+scale*viewVec.y,previous_location.z};
			object_location.push_back(projection);
			location_error.push_back(err);
			match_target.push_back(numMatch);
			scan_line.push_back(cam.center);
			scan_line.push_back(projection);
		}
	}
	printf("Loaded %lu reference cameras\n",camera_location.size());
	fclose(ref_camera_loc);

	SDL_Init(SDL_INIT_VIDEO);
	SDL_WM_SetCaption("Point Cloud", NULL);
	SDL_SetVideoMode(1800,1000, 32, SDL_OPENGL);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(70,(double)640/480,1,1000);

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
						case 'n':
						if (location_index >= object_location.size())
							location_index = 0;
						printf("location_index: %d error: %f\n",location_index,location_error[location_index]);
						if (isValid(location_index)) {
							current_location = object_location[location_index];
							translatePCD(object,
								current_location.x-previous_location.x,
								current_location.y-previous_location.y,
								current_location.z-previous_location.z);
							previous_location = current_location;
						}
						location_index++;
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

}
