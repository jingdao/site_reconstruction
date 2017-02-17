#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <vector>
#define TOP_VIEW 1

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
struct Triangle {
	size_t id1, id2, id3;
};

double cameraX,cameraY,cameraZ;
double centerX=0,centerY=0,centerZ=0;
#if TOP_VIEW
	double upX=0,upY=1,upZ=0;
#else
	double upX=0,upY=0,upZ=1;
#endif
int mouseIndex = 0;
int previousX,previousY;
double scrollSpeed = 1.01;
float fov = 70;
PCD* cloud;
std::vector<Point> modelVertices;
std::vector<Triangle> modelFaces;
SDL_Surface *screen;

bool loadPLY(const char* filename) {
	modelVertices.clear();
	modelFaces.clear();
	FILE* f = fopen(filename, "r");
	if (!f) {
		printf("File not found: %s\n", filename);
		return false;
	}
	char buf[256];
	int numVertex, numFace;
	while (fgets(buf, 256, f)) {
		if (sscanf(buf, "element vertex %d", &numVertex) == 1) {
		}
		else if (sscanf(buf, "element face %d", &numFace) == 1) {
		}
		else if (strncmp(buf, "end_header", 10) == 0) {
			for (int i = 0; i<numVertex; i++) {
				fgets(buf, 256, f);
				Point p = {0,0,0};
				if (sscanf(buf, "%f %f %f", &(p.x), &(p.y), &(p.z)) == 3) {
					p.x *= 0.0254;
					p.y *= 0.0254;
					p.z *= 0.0254;
					modelVertices.push_back(p);
				}
				else {
					printf("Error parsing %s\n", filename);
					printf("Line %d: %s\n", i, buf);
					break;
				}
			}
			for (int i = 0; i<numFace; i++) {
				fgets(buf, 256, f);
				Triangle t;
				if (sscanf(buf, "3 %lu %lu %lu", &(t.id1), &(t.id2), &(t.id3)) == 3) {
					modelFaces.push_back(t);
				}
				else {
					printf("Error parsing %s\n", filename);
					printf("Line %d: %s\n", i, buf);
					break;
				}
			}
			break;
		}
	}
	fclose(f);
	printf("Loaded %lu vertices %lu triangles from %s\n", modelVertices.size(), modelFaces.size(), filename);
	float centerX = 0;
	float centerY = 0;
	float bottomZ = 0;
	for (size_t i = 0; i < modelVertices.size(); i++) {
		centerX += modelVertices[i].x;
		centerY += modelVertices[i].y;
		if (i==0 || modelVertices[i].z < bottomZ)
			bottomZ = modelVertices[i].z;
	}
	centerX /= modelVertices.size();
	centerY /= modelVertices.size();
	centerX -= 10; 
	centerY += 5;
	for (size_t i = 0; i < modelVertices.size(); i++) {
		modelVertices[i].x -= centerX;
		modelVertices[i].y -= centerY;
		modelVertices[i].z -= bottomZ;
	}
	return true;
}

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

void getPCA(PCD *cloud, std::vector<Point> *box, Box extent) {
	float cx = (extent.minX + extent.maxX)/2;
	float cy = (extent.minY + extent.maxY)/2;
	float cz = extent.minZ;
	double cov[4] = {}; //column major
	Point center = {};
	for (int i = 0; i < cloud->numPoints; i++) {
		center.x += cloud->float_data[i*4];
		center.y += cloud->float_data[i*4+1];
		center.z += cloud->float_data[i*4+2];
	}
	center.x /= cloud->numPoints; 
	center.y /= cloud->numPoints;
	center.z /= cloud->numPoints;
	for (int j = 0; j<cloud->numPoints; j++) {
		float deltaP[2] = {
			cloud->float_data[j*4]- center.x,
			cloud->float_data[j*4+1]- center.y,
		};
		cov[0] += deltaP[0] * deltaP[0];
		cov[1] += deltaP[1] * deltaP[0];
		cov[2] += deltaP[0] * deltaP[1];
		cov[3] += deltaP[1] * deltaP[1];
	}
	cov[0] /= cloud->numPoints * cloud->numPoints;
	cov[1] /= cloud->numPoints * cloud->numPoints;
	cov[2] /= cloud->numPoints * cloud->numPoints;
	cov[3] /= cloud->numPoints * cloud->numPoints;
	float trace = cov[0] + cov[3];
	float det = cov[0] * cov[3] - cov[1] * cov[2];
	float L1 = trace / 2 + sqrt(trace*trace / 4 - det);
	float L2 = trace / 2 - sqrt(trace*trace / 4 - det);
	float minScale[3], maxScale[3];
	float v[9] = {
		0,0,0,
		0,0,0,
		0,0,1
	};
	if (cov[2] != 0) {
		v[0] = L1 - cov[3];
		v[1] = L2 - cov[3];
		v[3] = v[4] = cov[2];
	}
	else if (cov[1] != 0) {
		v[0] = v[1] = cov[1];
		v[3] = L1 - cov[0];
		v[4] = L2 - cov[0];
	}
	else {
		v[0] = v[4] = 1;
	}
	float m1 = sqrt(v[0] * v[0] + v[3] * v[3]);
	float m2 = sqrt(v[1] * v[1] + v[4] * v[4]);
	v[0] /= m1;
	v[3] /= m1;
	v[1] /= m2;
	v[4] /= m2;
	for (int j = 0; j<cloud->numPoints; j++) {
		for (int i = 0; i<3; i++) {
			float dotProduct =
				cloud->float_data[j*4] * v[i * 3] +
				cloud->float_data[j*4+1] * v[i * 3 + 1] +
				cloud->float_data[j*4+2] * v[i * 3 + 2];
			if (j == 0 || dotProduct < minScale[i])
				minScale[i] = dotProduct;
			if (j == 0 || dotProduct > maxScale[i])
				maxScale[i] = dotProduct;
		}
	}
	float bbCenter[3] = {0,0,0};
	for (int i = 0; i<3; i++) {
		bbCenter[0] += (minScale[i] + maxScale[i]) / 2 * v[i * 3];
		bbCenter[1] += (minScale[i] + maxScale[i]) / 2 * v[i * 3 + 1];
		bbCenter[2] += (minScale[i] + maxScale[i]) / 2 * v[i * 3 + 2];
	}
	for (int i = 0; i<8; i++) {
		float coords[3];
		for (int j = 0; j<3; j++) {
			coords[j] = bbCenter[j];
			for (int axis = 0; axis<3; axis++) {
				float sign = (i & 1 << axis) ? 1 : -1;
				coords[j] += sign * (maxScale[axis]-minScale[axis]) / 2 * v[axis * 3 + j];
			}
		}
		Point p = {coords[0]-cx,coords[1]-cy,coords[2]-cz};
		box->push_back(p);
	}
}

void drawLine(Point p1,Point p2) {
	glBegin(GL_LINES);
	glVertex3d(p1.x,p1.y,p1.z);
	glVertex3d(p2.x,p2.y,p2.z);
	glEnd();
}

void drawBox(Point* p) {
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

	glLineWidth(1.0);
	glColor3ub(200, 200, 200);
	glBegin(GL_LINES);
	for (size_t i = 0; i < modelFaces.size(); i++) {
		Point p1 = modelVertices[modelFaces[i].id1];
		Point p2 = modelVertices[modelFaces[i].id2];
		Point p3 = modelVertices[modelFaces[i].id3];
		glVertex3f(p1.x, p1.y, p1.z);
		glVertex3f(p2.x, p2.y, p2.z);
		glVertex3f(p2.x, p2.y, p2.z);
		glVertex3f(p3.x, p3.y, p3.z);
		glVertex3f(p3.x, p3.y, p3.z);
		glVertex3f(p1.x, p1.y, p1.z);
	}
	glEnd();

	glBegin(GL_TRIANGLES);
	glColor3ub(100, 100, 100);
	for (size_t i = 0; i < modelFaces.size(); i++) {
		Point p1 = modelVertices[modelFaces[i].id1];
		Point p2 = modelVertices[modelFaces[i].id2];
		Point p3 = modelVertices[modelFaces[i].id3];
		glVertex3f(p1.x, p1.y, p1.z);
		glVertex3f(p2.x, p2.y, p2.z);
		glVertex3f(p3.x, p3.y, p3.z);
	}
	glEnd();

	glFlush();
	SDL_GL_SwapBuffers();

	glPopMatrix();
	glPopAttrib();
}

void writeImageByIndex(int index,SDL_Surface *surf) {
	char buffer[128];
	sprintf(buffer,"%d-pcd.ppm",index);
	FILE* f = fopen(buffer,"w");
	fprintf(f,"P6\n%d %d\n255\n",surf->w,surf->h);
	unsigned char* pixels = new unsigned char[surf->w*surf->h*3]();
	glReadPixels(0,0,surf->w,surf->h,GL_RGB,GL_UNSIGNED_BYTE,pixels);
	unsigned char* src = pixels + (surf->h-1)*surf->w*3;
	for (int i=0;i<surf->h;i++) {
		fwrite(src,1,surf->w*3,f);
		src -= surf->w*3;
	}
	delete[] pixels;
	fclose(f);
}

int main(int argc,char* argv[]) {
	if (argc < 2) {
		printf("./%s scenario.pcd\n",argv[0]);
		return 1;
	}

	loadPLY("backhoe2.ply");
	char buffer[1024];
	cloud = NewPCD(argv[1]);
	if (!cloud)
		return 1;
	Box extent = getBoundingBox(cloud);
	centerPCD(cloud,extent);
#if TOP_VIEW
	cameraX = centerX = 12;
	cameraY = centerY = -5;
	cameraZ = 20;
#else
	cameraX = (extent.maxX - extent.minX)/2;
	cameraY = (extent.maxY - extent.minY)/2;
	cameraZ = (extent.maxZ - extent.minZ)*2;
#endif

	SDL_Init(SDL_INIT_VIDEO);
	SDL_WM_SetCaption("Point Cloud", NULL);
	screen = SDL_SetVideoMode(800,600, 24, SDL_OPENGL);
    glEnable(GL_DEPTH_TEST);
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
						centerX -= 1;
						draw();
						break;
						case SDLK_RIGHT:
						cameraX += 1;
						centerX += 1;
						draw();
						break;
						case SDLK_UP:
						cameraY += 1;
						centerY += 1;
						draw();
						break;
						case SDLK_DOWN:
						cameraY -= 1;
						centerY -= 1;
						draw();
						break;
						case 'n':
						break;
						case 'm':
						break;
						default:
						break;
					}
					break;
				case SDL_MOUSEBUTTONDOWN:
					if (event.button.button == SDL_BUTTON_WHEELUP) {
#if !TOP_VIEW
						cameraX /= scrollSpeed;
						cameraY /= scrollSpeed;
#endif
						cameraZ /= scrollSpeed;
						draw();
					} else if (event.button.button == SDL_BUTTON_WHEELDOWN) {
#if !TOP_VIEW
						cameraX *= scrollSpeed;
						cameraY *= scrollSpeed;
#endif
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
