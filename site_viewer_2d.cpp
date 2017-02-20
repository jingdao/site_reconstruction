#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <vector>
#include <ft2build.h>
#include FT_FREETYPE_H
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

const int labelWidth=80, labelHeight=30,fontpixels=20,grayLevel=50;
FT_Library ft;
FT_Face face;
FT_GlyphSlot glyph;
GLuint textures;
unsigned char raster[labelWidth * labelHeight * 3];
double cameraX,cameraY,cameraZ;
double centerX=0,centerY=0,centerZ=0;
double upX=0,upY=0,upZ=1;
int mouseIndex = 0;
int previousX,previousY;
double scrollSpeed = 1.01;
float fov = 70;
PCD* cloud;
std::vector< std::vector<Point> > object;
std::vector< std::vector<Point> > background;
int location_index = -1;
SDL_Surface *screen;

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

Point project3D(float u,float v,float z,Camera cam,float cx, float cy) {
	Point p;
	p.x = cam.position.i + (u-cx)*(cam.position.k-z)/cam.focal_length;
	p.y = cam.position.j - (v-cy)*(cam.position.k-z)/cam.focal_length;
	p.z = z;
	return p;
}

void render_text(const char *text, unsigned char *data) {
	memset(data,0,labelWidth*labelHeight*3);
	const char *p;
	unsigned char color[] = {255,255,255};
	int x = 0;
	for(p = text; *p; p++) {
		if(FT_Load_Char(face, *p, FT_LOAD_RENDER))
			continue;
		unsigned char *src = glyph->bitmap.buffer;
		int k=0;
		for (int i=0;i<glyph->bitmap.rows;i++) {
			unsigned char *dest = data + ((glyph->bitmap_top - i + fontpixels/2) * labelWidth + x + glyph->bitmap_left)* 3;
			for (int j=0;j<glyph->bitmap.width;j++) {
				memset(dest,*src,3); // draw in grayscale
				//*dest = *src; //draw in red
				src++;
				dest+=3;
			}
		}
		x += glyph->bitmap_left + glyph->bitmap.width;
//		x += glyph->bitmap.width;
		if (x >= labelWidth)
			break;
	}
}

void drawText(Point* p, const char* description) {
	glRasterPos3f((p[4].x+p[5].x+p[6].x+p[7].x)/4, (p[4].y+p[5].y+p[6].y+p[7].y)/4, p[4].z);
	render_text(description,raster);
	glDrawPixels(labelWidth,labelHeight,GL_RGB,GL_UNSIGNED_BYTE,raster);
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

	glColor3ub(255,255,255);
	glLineWidth(2.0);
	for (size_t i=0;i<background.size();i++) {
		drawBox(background[i].data());
	}

	if (location_index >= 0) {
		glLineWidth(5.0);
		for (size_t i=0;i<object[location_index].size()/8;i++) {
			if (i==0) {
				glColor3ub(0,255,0);
				drawBox(object[location_index].data() + i*8);
				drawText(object[location_index].data() + i*8, "car");
			} else {
				glColor3ub(0,0,255);
				drawBox(object[location_index].data() + i*8);
				drawText(object[location_index].data() + i*8, "worker");
			}
		}
	}
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
	if (argc < 3) {
		printf("./site_viewer scenario.pcd camera_location_2d.txt\n");
		return 1;
	}

	char buffer[1024];
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
	FILE* ppm = fopen("site.ppm","r");
	fgets(buffer,128,ppm);
	fgets(buffer,128,ppm);
	char *c = buffer;
	int bmp_width = strtol(c,&c,10);
	int bmp_height = strtol(c,&c,10);
	fclose(ppm);

	int background_id=0;
	while (true) {
		sprintf(buffer,"../clusters/%d-cloud.pcd",background_id);
		PCD* cloud = NewPCD(buffer);
		if (!cloud)
			break;
		std::vector<Point> box;
		getPCA(cloud,&box,extent);
		background.push_back(box);
		background_id++;
	}

	FILE* target_point = fopen(argv[2],"r");
	FILE* target_point_3d = fopen("target_point_3d.txt","w");
	if (!target_point) {
		printf("Cannot open %s\n",argv[2]);
		return 1;
	}
	while (fgets(buffer,1024,target_point)) {
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
			fprintf(target_point_3d,"4 %.3f %.3f 1 %.3f %.3f 1 %.3f %.3f 1 %.3f %.3f 1\n",p1.x,p1.y,p2.x,p2.y,p3.x,p3.y,p4.x,p4.y);
		}
		object.push_back(v);
	}
	fclose(target_point);
	fclose(target_point_3d);

	FT_Init_FreeType(&ft);
	FT_New_Face(ft,"/usr/share/fonts/truetype/freefont/FreeSans.ttf",0,&face);
	FT_Set_Pixel_Sizes(face,fontpixels,fontpixels);
	glyph = face->glyph;

	SDL_Init(SDL_INIT_VIDEO);
	SDL_WM_SetCaption("Point Cloud", NULL);
	screen = SDL_SetVideoMode(800,600, 24, SDL_OPENGL);
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
						draw();
						break;
						case SDLK_RIGHT:
						cameraX += 1;
						draw();
						break;
						case SDLK_UP:
						cameraZ += 1;
						draw();
						break;
						case SDLK_DOWN:
						cameraZ -= 1;
						draw();
						break;
						case 'n':
						do {
							location_index++;
							if (location_index >= object.size()) {
								location_index = -1;
								break;
							}
						} while (object[location_index].size()==0);
						draw();
						if (location_index >= 0)
							writeImageByIndex(location_index+1,screen);
						break;
						case 'm':
						while (true) {
							do {
								location_index++;
								if (location_index >= object.size()) {
									location_index = -1;
									break;
								}
							} while (object[location_index].size()==0);
							draw();
							if (location_index >= 0)
								writeImageByIndex(location_index+1,screen);
							else break;
							double rho = sqrt(cameraX*cameraX+cameraY*cameraY);
							double xstep = cameraY / rho * 50;
							double ystep = -cameraX / rho * 50;
							cameraX += 0.05 * xstep;
							cameraY += 0.05 * ystep;
						}
						break;
						default:
						break;
					}
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
