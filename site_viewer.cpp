#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <vector>
#define DEBUG 1
#define DRAW_TRAJECTORY 0

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
float fov = 70;
PCD* cloud;
std::vector<Point> object_location;
int location_index = 0;
std::vector<CamModel> camera_location;
std::vector<float> location_error;
std::vector<int> match_target;
std::vector<Point> scan_line;
std::vector<Point> trajectory;
std::vector<PCD*> object;
std::vector<Box> object_box;

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

PCD* NewPCD(const char* fileName,int num_intercept) {
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
				bool intercept = false;
				for (size_t i=0;i<num_intercept;i++) {
					if (x>object_box[i].minX && x<object_box[i].maxX &&
						y>object_box[i].minY && y<object_box[i].maxY &&
						z>object_box[i].minZ && z<object_box[i].maxZ)
						intercept=true;
				}
				if (!intercept) {
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

void drawBox(Box b) {
	glLineWidth(3.0);
	glColor3ub(0,255,0);
	glBegin(GL_LINES);

	glVertex3d(b.minX,b.minY,b.minZ);
	glVertex3d(b.maxX,b.minY,b.minZ);
	glVertex3d(b.minX,b.minY,b.minZ);
	glVertex3d(b.minX,b.maxY,b.minZ);
	glVertex3d(b.minX,b.maxY,b.minZ);
	glVertex3d(b.maxX,b.maxY,b.minZ);
	glVertex3d(b.maxX,b.minY,b.minZ);
	glVertex3d(b.maxX,b.maxY,b.minZ);

	glVertex3d(b.minX,b.minY,b.minZ);
	glVertex3d(b.minX,b.minY,b.maxZ);
	glVertex3d(b.maxX,b.minY,b.minZ);
	glVertex3d(b.maxX,b.minY,b.maxZ);
	glVertex3d(b.minX,b.maxY,b.minZ);
	glVertex3d(b.minX,b.maxY,b.maxZ);
	glVertex3d(b.maxX,b.maxY,b.minZ);
	glVertex3d(b.maxX,b.maxY,b.maxZ);

	glVertex3d(b.minX,b.minY,b.maxZ);
	glVertex3d(b.maxX,b.minY,b.maxZ);
	glVertex3d(b.minX,b.minY,b.maxZ);
	glVertex3d(b.minX,b.maxY,b.maxZ);
	glVertex3d(b.minX,b.maxY,b.maxZ);
	glVertex3d(b.maxX,b.maxY,b.maxZ);
	glVertex3d(b.maxX,b.minY,b.maxZ);
	glVertex3d(b.maxX,b.maxY,b.maxZ);

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

	glPointSize(1.0);
	glBegin(GL_POINTS);
	if (location_index == 0)
		glColor3ub(255,255,0);
	else if (location_error[location_index-1] < 0 || match_target[location_index-1] <=2)
		glColor3ub(255,0,0);
	else {
		double intensity = 255.0 / ((int)location_error[location_index-1] + 1);
		glColor3ub(0,(unsigned char)intensity,0);
	}

	for (size_t i=0;i<object.size();i++) {
		for (int n = 0; n < object[i]->numPoints; n++){
			glColor3ub(object[i]->color_data[n*3],object[i]->color_data[n*3+1],object[i]->color_data[n*3+2]);
			glVertex3d(object[i]->float_data[n*4],object[i]->float_data[n*4+1],object[i]->float_data[n*4+2]);
		}
	}
	glEnd();
	for (size_t i=0;i<object_box.size();i++)
		drawBox(object_box[i]);

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
		int li = location_index-object.size();
		drawLine(camera_location[li].center,camera_location[li].ul);
		drawLine(camera_location[li].center,camera_location[li].ur);
		drawLine(camera_location[li].center,camera_location[li].bl);
		drawLine(camera_location[li].center,camera_location[li].br);
		drawLine(camera_location[li].ul,camera_location[li].bl);
		drawLine(camera_location[li].ur,camera_location[li].br);
		drawLine(camera_location[li].ul,camera_location[li].ur);
		drawLine(camera_location[li].bl,camera_location[li].br);
#if DEBUG
		glColor3ub(255,0,0);
		for (size_t i=0;i<object_box.size();i++) {
			int li = location_index-object.size();
			if (isValid(li+i)) {
				drawLine(scan_line[2*(li+i)],scan_line[2*(li+i)+1]);
			}
		}
#endif
	}

#if DRAW_TRAJECTORY
	glColor3ub(0,0,255);
	for (int i=object.size();i<trajectory.size();i++) {
		drawLine(trajectory[i],trajectory[i - object.size()]);
	}
#endif

	glFlush();
	SDL_GL_SwapBuffers();

	glPopMatrix();
	glPopAttrib();
}

Box translateBox(Box b,float dx,float dy,float dz) {
	b.minX += dx;
	b.maxX += dx;
	b.minY += dy;
	b.maxY += dy;
	b.minZ += dz;
	b.maxZ += dz;
	return b;
}

void transformPCD(PCD* pcd,float* R,float* T) {
	for (int i=0;i<pcd->numPoints;i++) {
		Point p = {pcd->float_data[i*4],pcd->float_data[i*4+1],pcd->float_data[i*4+2]};
		pcd->float_data[i*4] = R[0]*p.x + R[1]*p.y + R[2]*p.z + T[0]; 
		pcd->float_data[i*4+1] = R[3]*p.x + R[4]*p.y + R[5]*p.z + T[1]; 
		pcd->float_data[i*4+2] = R[6]*p.x + R[7]*p.y + R[8]*p.z + T[2]; 
	}
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
		printf("./site_viewer scenario.pcd object.txt camera_location.txt [cam_settings.txt]\n");
		return 1;
	}

	char buffer[512];
	FILE* object_data = fopen(argv[2],"r");
	std::vector<Point> previous_location;
	std::vector<Point> current_location;
	std::vector<float> object_rotation;
	while (fgets(buffer,512,object_data)) {
		char filename[512];
		float rx,ry,rz;
		if (sscanf(buffer,"%s %f %f %f",filename,&rx,&ry,&rz)==4) {
			PCD* obj = NewPCD(filename,0);
			if (!obj)
				return 1;
			Box ob = getBoundingBox(obj);
			printf("Loaded %s (%f %f %f %f %f %f)\n",filename,ob.minX,ob.maxX,ob.minY,ob.maxY,ob.minZ,ob.maxZ);
			Point pl = {
				(ob.maxX + ob.minX)/2,
				(ob.maxY + ob.minY)/2,
				(ob.maxZ + ob.minZ)/2,
			};
			object.push_back(obj);
			object_box.push_back(ob);
			object_rotation.push_back(1.0 * rx / 180 * M_PI);
			object_rotation.push_back(1.0 * ry / 180 * M_PI);
			object_rotation.push_back(1.0 * rz / 180 * M_PI);
			previous_location.push_back(pl);
			current_location.push_back(pl);
		}
	}
	fclose(object_data);
	cloud = NewPCD(argv[1],object_box.size());
	if (!cloud)
		return 1;

	for (size_t i=0;i<object.size();i++) {
		float R[9];
		float T[3] = {previous_location[i].x,previous_location[i].y,previous_location[i].z};
		Point centroid = previous_location[i];
		translatePCD(object[i],-T[0],-T[1],-T[2]);
		float rx = object_rotation[3*i];
		float ry = object_rotation[3*i+1];
		float rz = object_rotation[3*i+2];
		rotationFromAngle(R,rx,ry,rz);
		transformPCD(object[i],R,T);
		object_box[i] = getBoundingBox(object[i]);
	}

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
	fgets(buffer,512,ref_camera_loc);
	int numScenes,numObjects;
	sscanf(buffer,"%d %d",&numScenes,&numObjects);
	for (int i=0;i<numScenes;i++) {
		fgets(buffer,512,ref_camera_loc);
		char *c = buffer;
		T[0] = strtod(c,&c);
		T[1] = strtod(c,&c);
		T[2] = strtod(c,&c);
		rx = strtod(c,&c);
		ry = strtod(c,&c);
		rz = strtod(c,&c);
		err = strtod(c,&c);
		rotationFromAngle(R,rx,ry,rz);
		CamModel cam;
		cam.center = transformPoint(center,R,T);
		cam.ul = transformPoint(ul,R,T);
		cam.ur = transformPoint(ur,R,T);
		cam.bl = transformPoint(bl,R,T);
		cam.br = transformPoint(br,R,T);
		for (int j=0;j<numObjects;j++) {
			viewVec.x = strtod(c,&c);
			viewVec.y = strtod(c,&c);
			viewVec.z = strtod(c,&c);
			numMatch = strtol(c,&c,10);
			float scale = (previous_location[j].z - T[2]) / viewVec.z;
			Point projection = {T[0]+scale*viewVec.x,T[1]+scale*viewVec.y,previous_location[j].z};
			camera_location.push_back(cam);
			object_location.push_back(projection);
			location_error.push_back(err);
			match_target.push_back(numMatch);
			scan_line.push_back(cam.center);
			scan_line.push_back(projection);
		}
	}
	printf("Loaded %lu reference cameras\n",camera_location.size());
	fclose(ref_camera_loc);

	if (argc > 4) {
		FILE* cam_settings = fopen(argv[4],"r");
		if (cam_settings) {
			float yaw,pitch,roll;
			float R[9];
			float T[3] = {0,0,0};
			Point upVector = {0,0,1};
			Point centerVector = {0,1,0};
			fscanf(cam_settings,"%lf %lf %lf %f %f %f",&cameraX,&cameraY,&cameraZ,&yaw,&pitch,&roll);
			rotationFromAngle(R,0,0,yaw);
			upVector = transformPoint(upVector,R,T);
			centerVector = transformPoint(centerVector,R,T);
			upX = upVector.x;
			upY = upVector.y;
			upZ = upVector.z;
			centerX = cameraX + centerVector.x;
			centerY = cameraY + centerVector.y;
			centerZ = cameraZ + centerVector.z;
			fclose(cam_settings);
		}
	} else {
		Box cloud_box = getBoundingBox(cloud);
		upX = 0; upY = -1; upZ = 0;
		centerX = (cloud_box.minX + cloud_box.maxX) / 2;
		centerY = (cloud_box.minY + cloud_box.maxY) / 2;
		centerZ = (cloud_box.minZ + cloud_box.maxZ) / 2;
		cameraX = centerX;
		cameraY = centerY;
		cameraZ = centerZ + (cloud_box.maxY - cloud_box.minY) / 2 / tan(fov/360*M_PI);
	}

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
						for (int i=0;i<numObjects;i++) {
							if (location_index >= object_location.size())
								location_index = 0;
							printf("location_index: %d error: %f\n",location_index,location_error[location_index]);
							if (isValid(location_index)) {
								current_location[i] = object_location[location_index];
								trajectory.push_back(current_location[i]);
								translatePCD(object[i],
									current_location[i].x-previous_location[i].x,
									current_location[i].y-previous_location[i].y,
									current_location[i].z-previous_location[i].z);
								object_box[i] = translateBox(object_box[i],
									current_location[i].x-previous_location[i].x,
									current_location[i].y-previous_location[i].y,
									current_location[i].z-previous_location[i].z);
								previous_location[i] = current_location[i];
							}
							location_index++;
						}
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
