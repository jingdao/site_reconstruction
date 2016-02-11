#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <SDL/SDL.h>
#include <unistd.h>

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
struct Point {
	float x,y,z;
};
struct Box {
	float minX,minY,minZ,maxX,maxY,maxZ;
};
struct Image {
	int width,height;
	unsigned char* data;
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

Quaternion quaternionMult(Quaternion qa, Quaternion qb) {
	Quaternion qc;
	qc.r = qa.r*qb.r - qa.i*qb.i - qa.j*qb.j - qa.k*qb.k;
	qc.i = qa.r*qb.i + qa.i*qb.r + qa.j*qb.k - qa.k*qb.j;
	qc.j = qa.r*qb.j + qa.j*qb.r + qa.k*qb.i - qa.i*qb.k;
	qc.k = qa.r*qb.k + qa.k*qb.r + qa.i*qb.j - qa.j*qb.i;
	return qc;
}

Quaternion quaternionInv(Quaternion q) {
	Quaternion qc;
	qc.r = q.r;
	qc.i = -q.i;
	qc.j = -q.j;
	qc.k = -q.k;
	return qc;
}

Quaternion quaternionRot(Quaternion q, Quaternion vec) {
	return quaternionMult(quaternionMult(q, vec), quaternionInv(q));
}

Quaternion quaternionFromAngle(float rx,float ry, float rz) {
	float r,i,j,k;
	float r1 = rz * M_PI / 180;
	float r2 = ry * M_PI / 180;
	float r3 = rx * M_PI / 180;
	float M[9];
	M[0] = cos(r2)*cos(r3);
	M[1] = -cos(r2) * sin(r3);
	M[2] = sin(r2);
	M[3] = cos(r1) * sin(r3) + cos(r3) * sin(r1) * sin(r2);
	M[4] = cos(r1) * cos(r3) - sin(r1) * sin(r2) * sin(r3);
	M[5] = -cos(r2) * sin(r1);
	M[6] = sin(r1) * sin(r3) - cos(r1) * cos(r3) * sin(r2);
	M[7] = cos(r3) * sin(r1) + cos(r1) * sin(r2) * sin(r3);
	M[8] = cos(r1) * cos(r2);

	double tr = M[0] + M[4] + M[8];

	if (tr > 0) { 
		float S = sqrt(tr+1.0) * 2;
		r = (float)(0.25 * S);
		i = (float)((M[7] - M[5]) / S);
		j = (float)((M[2] - M[6]) / S);
		k = (float)((M[3] - M[1]) / S);
	} else if ((M[0] > M[4])&&(M[0]> M[8])) {
		float S = sqrt(1.0 + M[0]- M[4] - M[8]) * 2;
		i = (float)(0.25 * S);
		r = (float)((M[7] - M[5]) / S);
		k = (float)((M[2] + M[6]) / S);
		j = (float)((M[3] + M[1]) / S);
	} else if (M[4]>M[8]) {
		float S = sqrt(1.0 - M[0] + M[4] - M[8]) * 2; 
		j = (float)(0.25 * S);
		k = (float)((M[7] + M[5]) / S);
		r = (float)((M[2] - M[6]) / S);
		i = (float)((M[3] + M[1]) / S); 
	} else {
		float S = sqrt(1.0 - M[0] - M[4] + M[8]) * 2; 
		k = (float)(0.25 * S);
		j = (float)((M[7] + M[5]) / S);
		i = (float)((M[2] + M[6]) / S);
		r = (float)((M[3] - M[1]) / S); 
	}
	Quaternion res = {r,i,j,k};
	return res;
}

int SetRGB(Image* bmp, int x, int y, unsigned char r, unsigned char g, unsigned char b) {
	if (!bmp || x < 0 || y < 0 || x >= bmp->width || y >= bmp->height) return 0;
	int offset = (x + y*bmp->width) * 3;
	bmp->data[offset] = r;
	bmp->data[offset + 1] = g;
	bmp->data[offset + 2] = b;
	return 1;
}

void Show3DProjection(Image* bmp, PCD* pcd, Camera camera, Point* depth_buffer) {
	memset(bmp->data,0,bmp->width*bmp->height*3);
	memset(depth_buffer,0,bmp->width*bmp->height*sizeof(Point));
	int i;
	Quaternion q, newPoint;
	std::vector<float> xl,yl,zl;
	std::vector<int> index; 
	q.r = 0;
	for (i = 0; i < pcd->numPoints; i++) {
		//first transform point cloud to camera coordinates
		q.i = pcd->float_data[i * 4] - camera.position.i;
		q.j = pcd->float_data[i * 4 + 1] - camera.position.j;
		q.k = pcd->float_data[i * 4 + 2] - camera.position.k;
		newPoint = quaternionRot(camera.rotation, q);
		//next map points to image plane based on camera focal length
		if (newPoint.k > 0) { //if not behind camera
			xl.push_back(newPoint.i);
			yl.push_back(newPoint.j);
			zl.push_back(newPoint.k);
			index.push_back(i);
		}
	}
	//lastly plot these points on xy plane of bitmap
	unsigned char r, g, b;
	int x,y;
	int pixelsDrawn = 0;
	float maxX = -1;
	float maxY = -1;
	float scaleX=1, scaleY=1;
	for (unsigned int j=0;j<xl.size();j++) {
		x = (int) ((xl[j]*camera.focal_length/zl[j]) * scaleX + bmp->width / 2);
		y = (int) ((yl[j]*camera.focal_length/zl[j]) * scaleX + bmp->height / 2);
		if (x < 0 || y < 0 || x >= bmp->width || y >= bmp->height)
			continue;
		Point* p = depth_buffer + y*bmp->width+x;
		if (p->z == 0 || p->z > zl[j]) {
			p->x = pcd->float_data[index[j] * 4];
			p->y = pcd->float_data[index[j] * 4 + 1];
			p->z = pcd->float_data[index[j] * 4 + 2];
			r = pcd->color_data[index[j]*3];
			g = pcd->color_data[index[j]*3+1];
			b = pcd->color_data[index[j]*3+2];
			pixelsDrawn += SetRGB(bmp, x, y, r, g, b);
		}
	}
	printf("(%.0f %.0f %.0f %f %f %f) %d pixels drawn\n",
		camera.position.i,camera.position.j,camera.position.k,
		camera.yaw,camera.pitch,camera.roll,
		pixelsDrawn);
}

void fillImage(Image* bmp,Point* depth_buffer) {
	bool updated = true;
	bool* filled = new bool[bmp->height * bmp->width];
	unsigned char* c = bmp->data;
	for (int i=0;i<bmp->height;i++) {
		for (int j=0;j<bmp->width;j++) {
			filled[j + i * bmp->width] = (c[0] || c[1] || c[2]);
			c += 3;
		}
	}
	//average neighbors
	while (updated) {
		updated = false;
		c = bmp->data + (bmp->width+1)*3;
		for (int i=1;i<bmp->height-1;i++) {
			for (int j=1;j<bmp->width-1;j++) {
				if (!filled[j + i * bmp->width]) {
					int r=0,g=0,b=0,count=0;
					std::vector<Point> neighbors;
					if (filled[j - 1 + i * bmp->width]) {
						r+=c[-3]; g+=c[-2]; b+=c[-1]; count++;
						neighbors.push_back(depth_buffer[j - 1 + i * bmp->width]);
					}
					if (filled[j + 1 + i * bmp->width]) {
						r+=c[3]; g+=c[4]; b+=c[5]; count++;
						neighbors.push_back(depth_buffer[j + 1 + i * bmp->width]);
					}
					if (filled[j + (i-1) * bmp->width]) {
						r+=c[-bmp->width*3]; g+=c[-bmp->width*3+1]; b+=c[-bmp->width*3+2]; count++;
						neighbors.push_back(depth_buffer[j + (i-1) * bmp->width]);
					}
					if (filled[j + (i+1) * bmp->width]) {
						r+=c[bmp->width*3]; g+=c[bmp->width*3+1]; b+=c[bmp->width*3+2]; count++;
						neighbors.push_back(depth_buffer[j + (i+1) * bmp->width]);
					}
					if (count) {
						filled[j + i * bmp->width] = true;
						c[0] = (unsigned char) (r/count);
						c[1] = (unsigned char) (g/count);
						c[2] = (unsigned char) (b/count);
						for (int k=0;k<count;k++) {
							depth_buffer[j + i * bmp->width].x += neighbors[k].x;
							depth_buffer[j + i * bmp->width].y += neighbors[k].y;
							depth_buffer[j + i * bmp->width].z += neighbors[k].z;
						}
						depth_buffer[j + i * bmp->width].x /= count;
						depth_buffer[j + i * bmp->width].y /= count;
						depth_buffer[j + i * bmp->width].z /= count;
					} else {
						updated = true;
					}
				}
				c += 3;
			}
			c += 6;
		}
	}
	//fill edges
	for (int i=1;i<bmp->width-1;i++) {
		bmp->data[i*3] = bmp->data[(bmp->width+i)*3];
		bmp->data[i*3+1] = bmp->data[(bmp->width+i)*3+1];
		bmp->data[i*3+2] = bmp->data[(bmp->width+i)*3+2];
		bmp->data[((bmp->height-1)*bmp->width+i)*3] = bmp->data[((bmp->height-2)*bmp->width+i)*3];
		bmp->data[((bmp->height-1)*bmp->width+i)*3+1] = bmp->data[((bmp->height-2)*bmp->width+i)*3+1];
		bmp->data[((bmp->height-1)*bmp->width+i)*3+2] = bmp->data[((bmp->height-2)*bmp->width+i)*3+2];
		depth_buffer[i] = depth_buffer[bmp->width+i];
		depth_buffer[(bmp->height-1)*bmp->width+i] = depth_buffer[(bmp->height-2)*bmp->width+i];
	}
	for (int i=0;i<bmp->height;i++) {
		bmp->data[i*bmp->width*3] = bmp->data[(i*bmp->width+1)*3];
		bmp->data[i*bmp->width*3+1] = bmp->data[(i*bmp->width+1)*3+1];
		bmp->data[i*bmp->width*3+2] = bmp->data[(i*bmp->width+1)*3+2];
		bmp->data[((i+1)*bmp->width-1)*3] = bmp->data[((i+1)*bmp->width-2)*3];
		bmp->data[((i+1)*bmp->width-1)*3+1] = bmp->data[((i+1)*bmp->width-2)*3+1];
		bmp->data[((i+1)*bmp->width-1)*3+2] = bmp->data[((i+1)*bmp->width-2)*3+2];
		depth_buffer[i*bmp->width] = depth_buffer[i*bmp->width+1];
		depth_buffer[(i+1)*bmp->width-1] = depth_buffer[(i+1)*bmp->width-2];
	}
	free(filled);
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
				pcd->float_data[pointsParsed * 4 + 3] = c;
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

void writePPM(char* filename,unsigned char* imageBuffer,int width,int height) {
	FILE* f = fopen(filename,"w");
	if (!f)
		return;
	fprintf(f,"P6\n%d %d\n255\n",width,height);
	fwrite(imageBuffer,1,width*height*3,f);
	printf("Saved to %s\n",filename);
	fclose(f);
}

void imgcpy(Image image,SDL_Surface* surf) {
	unsigned char* src = image.data, *dst = (unsigned char*) surf->pixels;
	for (int i=0;i<image.height;i++) {
		memcpy(dst,src,image.width*3);
		src += image.width*3;
		dst += surf->pitch;
	}
}

int main(int argc,char* argv[]) {

	if (argc < 6) {
		printf("./im_synth scenario.pcd width height site.ppm depth_buffer.txt\n");
		return 1;
	}

	PCD* site = NewPCD(argv[1],NULL);
	if (!site)
		return 1;
	
	int width=atoi(argv[2]);
	int height=atoi(argv[3]);
	float speed=5; 
	float max_depth,min_depth,depth_scale;
	Image ppm = {width,height,NULL};
	ppm.data = new unsigned char[width*height*3];
	SDL_Surface* screen = SDL_SetVideoMode(width,height,24,SDL_SWSURFACE);
	Point* depth_buffer = new Point[width*height]();
	FILE* depth_buffer_txt, *depth_buffer_pgm;

	Camera cam = {
		{0,30,40,30},
		{0,0,0,0},
		600,130,205,0
	};
	cam.rotation = quaternionFromAngle(cam.pitch,cam.roll,cam.yaw);

	SDL_Init(SDL_INIT_VIDEO);
	SDL_WM_SetCaption(argv[1],NULL);
	Show3DProjection(&ppm,site,cam,depth_buffer);
	imgcpy(ppm,screen);
	SDL_Flip(screen);

	SDL_Event event;
	while (true) {
		while (SDL_PollEvent(&event)) {
			switch(event.type){
				case SDL_KEYDOWN:
					switch( event.key.keysym.sym ){
						case SDLK_LEFT:
						cam.position.i -= speed;
						break;
						case SDLK_RIGHT:
						cam.position.i += speed;
						break;
						case SDLK_UP:
						cam.position.j -= speed;
						break;
						case SDLK_DOWN:
						cam.position.j += speed;
						break;
						case 'm':
						cam.position.k -= speed;
						break;
						case ',':
						cam.position.k += speed;
						break;
						case 'z':
						cam.yaw -= speed;
						cam.rotation = quaternionFromAngle(cam.pitch,cam.roll,cam.yaw);
						break;
						case 'x':
						cam.yaw += speed;
						cam.rotation = quaternionFromAngle(cam.pitch,cam.roll,cam.yaw);
						break;
						case 'c':
						cam.pitch -= speed;
						cam.rotation = quaternionFromAngle(cam.pitch,cam.roll,cam.yaw);
						break;
						case 'v':
						cam.pitch += speed;
						cam.rotation = quaternionFromAngle(cam.pitch,cam.roll,cam.yaw);
						break;
						case 'b':
						cam.roll -= speed;
						cam.rotation = quaternionFromAngle(cam.pitch,cam.roll,cam.yaw);
						break;
						case 'n':
						cam.roll += speed;
						cam.rotation = quaternionFromAngle(cam.pitch,cam.roll,cam.yaw);
						break;
						case 's':
						writePPM(argv[4],ppm.data,ppm.width,ppm.height);
						depth_buffer_txt = fopen(argv[5],"w");
						for (int i=0;i<width*height;i++) {
							fprintf(depth_buffer_txt,"%f %f %f\n",
								depth_buffer[i].x,
								depth_buffer[i].y,
								depth_buffer[i].z);
						}
						fclose(depth_buffer_txt);
						strcpy(strstr(argv[5],".txt"),".pgm");
						depth_buffer_pgm = fopen(argv[5],"w");
						fprintf(depth_buffer_pgm,"P5\n%d %d\n255\n",width,height);
						max_depth = depth_buffer[0].z;
						min_depth = depth_buffer[0].z;
						for (int i=1;i<width*height;i++) {
							if (depth_buffer[i].z > max_depth)
								max_depth = depth_buffer[i].z;
							if (depth_buffer[i].z < min_depth)
								min_depth = depth_buffer[i].z;
						}
						depth_scale = 255.0 / (max_depth - min_depth); 
						for (int i=0;i<width*height;i++) {
							unsigned char c = (unsigned char) ((depth_buffer[i].z-min_depth) * depth_scale);
							fwrite(&c,1,1,depth_buffer_pgm);
						}
						fclose(depth_buffer_pgm);
						break;
						default:
						break;
					}
					Show3DProjection(&ppm,site,cam,depth_buffer);
					imgcpy(ppm,screen);
					SDL_Flip(screen);
					break;
				case SDL_MOUSEBUTTONDOWN:
					fillImage(&ppm,depth_buffer);
					imgcpy(ppm,screen);
					SDL_Flip(screen);
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

	delete[] ppm.data;
	delete[] depth_buffer;
	free(site->float_data);
	free(site->color_data);

}
