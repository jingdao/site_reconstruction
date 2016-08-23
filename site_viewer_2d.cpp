#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <vector>
#define OBJ_HEIGHT 2

struct PCD {
	int numPoints;
	float* float_data;
	unsigned char* color_data;
};
struct Box {
	float minX,minY,minZ,maxX,maxY,maxZ;
};
enum PCD_data_storage {
	ASCII,
	BINARY,
	NONE
};
struct Point {
	float x,y,z;
};
struct Quaternion{
    float r;
    float i;
    float j;
    float k;
};
struct Camera{
    Quaternion position;
    Quaternion rotation;
    float focal_length,yaw,pitch,roll;
};

double cameraX,cameraY,cameraZ;
double centerX=0,centerY=0,centerZ=0;
double upX=0,upY=0,upZ=1;
int mouseIndex = 0;
int previousX,previousY;
double scrollSpeed = 1.01;
float fov = 70;
PCD* cloud;
std::vector< std::vector<Point> > object;
int location_index = 0;

PCD* NewPCD(const char* fileName) {
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
			if (numParsed >= 3) {
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

void centerPCD(PCD* pcd, Box box) {
	float cx = (box.minX + box.maxX)/2;
	float cy = (box.minY + box.maxY)/2;
	float cz = box.minZ;
	for (int n=0;n<pcd->numPoints;n++) {
		pcd->float_data[n*4] -= cx;
		pcd->float_data[n*4+1] -= cy;
		pcd->float_data[n*4+2] -= cz;
	}
}

void drawLine(Point p1,Point p2) {
	glBegin(GL_LINES);
	glVertex3d(p1.x,p1.y,p1.z);
	glVertex3d(p2.x,p2.y,p2.z);
	glEnd();
}

void drawBox(Point* p) {
	glLineWidth(3.0);
	glBegin(GL_LINES);

	glVertex3d(p[0].x,p[0].y,p[0].z);
	glVertex3d(p[1].x,p[1].y,p[1].z);
	glVertex3d(p[0].x,p[0].y,p[0].z);
	glVertex3d(p[2].x,p[2].y,p[2].z);
	glVertex3d(p[1].x,p[1].y,p[1].z);
	glVertex3d(p[3].x,p[3].y,p[3].z);
	glVertex3d(p[2].x,p[2].y,p[2].z);
	glVertex3d(p[3].x,p[3].y,p[3].z);

	glVertex3d(p[0].x,p[0].y,p[0].z);
	glVertex3d(p[4].x,p[4].y,p[4].z);
	glVertex3d(p[1].x,p[1].y,p[1].z);
	glVertex3d(p[5].x,p[5].y,p[5].z);
	glVertex3d(p[2].x,p[2].y,p[2].z);
	glVertex3d(p[6].x,p[6].y,p[6].z);
	glVertex3d(p[3].x,p[3].y,p[3].z);
	glVertex3d(p[7].x,p[7].y,p[7].z);

	glVertex3d(p[4].x,p[4].y,p[4].z);
	glVertex3d(p[5].x,p[5].y,p[5].z);
	glVertex3d(p[4].x,p[4].y,p[4].z);
	glVertex3d(p[6].x,p[6].y,p[6].z);
	glVertex3d(p[5].x,p[5].y,p[5].z);
	glVertex3d(p[7].x,p[7].y,p[7].z);
	glVertex3d(p[6].x,p[6].y,p[6].z);
	glVertex3d(p[7].x,p[7].y,p[7].z);

	glEnd();
}

Point project3D(float u,float v,float z,Camera cam,float cx, float cy) {
	Point p;
	p.x = cam.position.i + (u-cx)*(cam.position.k-z)/cam.focal_length;
	p.y = cam.position.j - (v-cy)*(cam.position.k-z)/cam.focal_length;
	p.z = z;
	return p;
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

	glColor3ub(0,255,0);
	for (size_t i=0;i<object[location_index].size()/8;i++) {
		drawBox(object[location_index].data() + i*8);
	}

	glFlush();
	SDL_GL_SwapBuffers();

	glPopMatrix();
	glPopAttrib();
}

int main(int argc,char* argv[]) {
	if (argc < 3) {
		printf("./site_viewer scenario.pcd camera_location_2d.txt\n");
		return 1;
	}

	char buffer[128];
	cloud = NewPCD(argv[1]);
	if (!cloud)
		return 1;
	Box extent = getBoundingBox(cloud);
	centerPCD(cloud,extent);
	cameraX = (extent.maxX - extent.minX)/2;
	cameraY = (extent.maxY - extent.minY)/2;
	cameraZ = extent.maxZ - extent.minZ;
	Camera cam = {
		{0,2,-1,60},
		{0,0,0,0},
		600,180,0,0
	};
	cam.position.i -= (extent.maxX + extent.minX) / 2;
	cam.position.j -= (extent.maxY + extent.minY) / 2;
	cam.position.k -= extent.minZ;
	FILE* ppm = fopen("site_blur.ppm","r");
	fgets(buffer,128,ppm);
	fgets(buffer,128,ppm);
	char *c = buffer;
	int bmp_width = strtol(c,&c,10);
	int bmp_height = strtol(c,&c,10);
	fclose(ppm);

	FILE* target_point = fopen(argv[2],"r");
	if (!target_point) {
		printf("Cannot open %s\n",argv[2]);
		return 1;
	}
	while (fgets(buffer,128,target_point)) {
		float x,y;
		int numMatch;
		char* c = buffer;
		int n = strtol(c,&c,10);
		std::vector<float> x_target;
		std::vector<float> y_target;
		std::vector<Point> v;
		for (int i=0;i<n;i++) {
			x = strtod(c,&c);
			y = strtod(c,&c);
			numMatch = strtol(c,&c,10);
			x_target.push_back(x);
			y_target.push_back(y);
		}
		for (int i=0;i<n/4;i++) {
			Point p1 = project3D(x_target[i*4],y_target[i*4],0,cam,bmp_width/2,bmp_height/2);
			Point p2 = project3D(x_target[i*4+1],y_target[i*4+1],0,cam,bmp_width/2,bmp_height/2);
			Point p3 = project3D(x_target[i*4+2],y_target[i*4+2],0,cam,bmp_width/2,bmp_height/2);
			Point p4 = project3D(x_target[i*4+3],y_target[i*4+3],0,cam,bmp_width/2,bmp_height/2);
			Point q1 = {p1.x,p1.y,OBJ_HEIGHT};
			Point q2 = {p2.x,p2.y,OBJ_HEIGHT};
			Point q3 = {p3.x,p3.y,OBJ_HEIGHT};
			Point q4 = {p4.x,p4.y,OBJ_HEIGHT};
			v.push_back(p1);
			v.push_back(p2);
			v.push_back(p3);
			v.push_back(p4);
			v.push_back(q1);
			v.push_back(q2);
			v.push_back(q3);
			v.push_back(q4);
		}
		object.push_back(v);
	}
	fclose(target_point);

	SDL_Init(SDL_INIT_VIDEO);
	SDL_WM_SetCaption("Point Cloud", NULL);
	SDL_SetVideoMode(800,600, 32, SDL_OPENGL);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fov,(double)640/480,1,1000);

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
						if (location_index >= object.size())
							location_index = 0;
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
