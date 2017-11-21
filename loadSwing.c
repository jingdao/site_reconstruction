#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <getopt.h>             /* getopt_long() */

#include <fcntl.h>              /* low-level i/o */
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/ioctl.h>

#include <linux/videodev2.h>
#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <math.h>
#include <stdbool.h>
#define USE_3D 0

#define CLEAR(x) memset(&(x), 0, sizeof(x))

#ifndef V4L2_PIX_FMT_H264
#define V4L2_PIX_FMT_H264     v4l2_fourcc('H', '2', '6', '4') /* H264 with start codes */
#endif

enum io_method {
        IO_METHOD_READ,
        IO_METHOD_MMAP,
        IO_METHOD_USERPTR,
};

struct buffer {
        void   *start;
        size_t  length;
};

typedef struct {
	int width,height;
	unsigned char* data;
} Image;

typedef struct {
	unsigned char r,g,b;
} Color;

static char            *dev_name;
static enum io_method   io = IO_METHOD_MMAP;
static int              fd = -1;
struct buffer          *buffers;
static unsigned int     n_buffers;
static int              out_buf;
static int              force_format;
static int              frame_count = 200;
static int              frame_number = 0;

SDL_Surface *screen;
unsigned char* rgb_data;
unsigned char* tmp_data = NULL;
bool* mask = NULL;

double cameraX=2,cameraY=-5,cameraZ=1.5;
double centerX=0,centerY=0,centerZ=0;
double upX=0,upY=0,upZ=1;
float fov = 70;
double base_x=0;
double base_y=0;

static void errno_exit(const char *s)
{
        fprintf(stderr, "%s error %d, %s\n", s, errno, strerror(errno));
        exit(EXIT_FAILURE);
}

static int xioctl(int fh, int request, void *arg)
{
        int r;

        do {
                r = ioctl(fh, request, arg);
        } while (-1 == r && EINTR == errno);

        return r;
}

void drawPendulum(double center_x, double center_y) {
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glPushMatrix();
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(cameraX,cameraY,cameraZ,centerX,centerY,centerZ,upX,upY,upZ);

	//project to world space
	base_x = 1.0*(frame_number-1)/frame_number*base_x + 1.0 / frame_number * center_x;
	base_y = 1.0*(frame_number-1)/frame_number*base_y + 1.0 / frame_number * center_y;
	center_x -= base_x;
	center_y -= base_y;
	double focal_length = 500;
	double cable_length = 1.32;
	double base_dim = 0.86;
	double base_height = 1.98;
	double k = center_x*center_x/focal_length/focal_length + center_y*center_y/focal_length/focal_length + 1;
	double world_z = sqrt(cable_length*cable_length / k);
	double world_x = center_x / focal_length * world_z;
	double world_y = center_y / focal_length * world_z;
	world_z = base_height - world_z;

	printf("Frame %d (%6.2f %6.2f,%6.2f)\r",frame_number, world_x,world_y,world_z);
	fflush(stdout);

	//draw grid
	double grid_dim = 2.0;
	glLineWidth(2.0);
	glBegin(GL_LINES);
	glColor3ub(0,255,0);

	glVertex3d(-grid_dim, -grid_dim, 0);
	glVertex3d(-grid_dim, grid_dim, 0);
	glVertex3d(-grid_dim, -grid_dim, 0);
	glVertex3d(grid_dim, -grid_dim, 0);
	glVertex3d(-grid_dim, grid_dim, 0);
	glVertex3d(grid_dim, grid_dim, 0);
	glVertex3d(grid_dim, -grid_dim, 0);
	glVertex3d(grid_dim, grid_dim, 0);

	glVertex3d(-grid_dim, -grid_dim, 0);
	glVertex3d(-grid_dim, -grid_dim, grid_dim);
	glVertex3d(-grid_dim, grid_dim, 0);
	glVertex3d(-grid_dim, grid_dim, grid_dim);
	glVertex3d(grid_dim, -grid_dim, 0);
	glVertex3d(grid_dim, -grid_dim, grid_dim);
	glVertex3d(grid_dim, grid_dim, 0);
	glVertex3d(grid_dim, grid_dim, grid_dim);

	glVertex3d(-grid_dim, -grid_dim, grid_dim);
	glVertex3d(-grid_dim, grid_dim, grid_dim);
	glVertex3d(-grid_dim, -grid_dim, grid_dim);
	glVertex3d(grid_dim, -grid_dim, grid_dim);
	glVertex3d(-grid_dim, grid_dim, grid_dim);
	glVertex3d(grid_dim, grid_dim, grid_dim);
	glVertex3d(grid_dim, -grid_dim, grid_dim);
	glVertex3d(grid_dim, grid_dim, grid_dim);


//	glColor3ub(255,0,0);
//	glVertex3d(0,0,0);
//	glVertex3d(grid_dim,0,0);
//	glColor3ub(0,255,0);
//	glVertex3d(0,0,0);
//	glVertex3d(0,grid_dim,0);
//	glColor3ub(0,0,255);
//	glVertex3d(0,0,0);
//	glVertex3d(0,0,grid_dim);

	//draw pendulum
	glColor3ub(255,255,255);
	glVertex3d(-base_dim,-base_dim,0);
	glVertex3d(0,0,base_height);
	glVertex3d(0,0,base_height);
	glVertex3d(world_x, world_y, world_z);
	glEnd();

	glPointSize(10.0);
	glBegin(GL_POINTS);
	glVertex3d(world_x,world_y,world_z);
	glEnd();

	glFlush();
	SDL_GL_SwapBuffers();
	glPopMatrix();
	glPopAttrib();
}

inline void yuv2rgb(unsigned char y,unsigned char u,unsigned char v,unsigned char *r,unsigned char *g,unsigned char *b) {
	int c = y - 16;
	int d = u - 128;
	int e = v - 128;

	int _r = (298 * c + 409 * e + 128) >> 8;
	int _g = (298 * c - 100 * d - 208 * e + 128) >> 8;
	int _b = (298 * c + 516 * d + 128) >> 8;

//	int cb = u - 128;
//	int cr = v - 128;
//
//	int _r = y + cr + (cr>>2) + (cr>>3) + (cr>>5);
//	int _g = y - (cb>>2) - (cb>>4) - (cb>>5) + (cr>>1) + (cr>>3) + (cr>>4) + (cr>>5);
//	int _b = y + cb + (cb>>1) + (cb>>2) + (cb>>6);

	_r = (_r>=0 ? _r : 0);
	_g = (_g>=0 ? _g : 0);
	_b = (_b>=0 ? _b : 0);

	*r = (_r<256 ? _r : 255);
	*g = (_g<256 ? _g : 255);
	*b = (_b<256 ? _b : 255);
}

void convertLAB(Image image) {
	unsigned char *src = image.data;
	for (int i=0;i<image.height;i++) {
		for (int j=0;j<image.width;j++) {
			float r = (*src++) / 255.0;
			float g = (*src++) / 255.0;
			float b = (*src++) / 255.0;
			r = (r > 0.04045 ? pow((r + 0.055) / 1.055, 2.4) : r / 12.92) * 100.0;
			g = (g > 0.04045 ? pow((g + 0.055) / 1.055, 2.4) : g / 12.92) * 100.0;
			b = (b > 0.04045 ? pow((b + 0.055) / 1.055, 2.4) : b / 12.92) * 100.0;
			float X = r * 0.4124 + g * 0.3576 + b * 0.1805;
			float Y = r * 0.2126 + g * 0.7152 + b * 0.0722;
			float Z = r * 0.0193 + g * 0.1192 + b * 0.9505;

			float var_X = X / 95.047;
			float var_Y = Y / 100;
			float var_Z = Z / 108.883;

			if ( var_X > 0.008856 ) var_X = pow(var_X,( 1.0/3 ));
			else var_X = (903.3*var_X + 16) / 116;
			if ( var_Y > 0.008856 ) var_Y = pow(var_Y,( 1.0/3 ));
			else var_Y = (903.3*var_Y + 16) / 116;
			if ( var_Z > 0.008856 ) var_Z = pow(var_Z,( 1.0/3 ));
			else var_Z = (903.3*var_Z + 16) / 116;

			float L = 2.56 * (( 116 * var_Y ) - 16);
			float A = 1.388 * 500 * ( var_X - var_Y ) + 119.624;
			float B = 1.26 * 200 * ( var_Y - var_Z ) + 135.932;
//			printf("%f %f %f\n",L,A,B);
			src[-3] = L < 0 ? 0 : L > 255 ? 255 : (unsigned char) L;
			src[-2] = A < 0 ? 0 : A > 255 ? 255 : (unsigned char) A;
			src[-1] = B < 0 ? 0 : B > 255 ? 255 : (unsigned char) B;
		}
	}
}

int getDiff(Color c1, Color c2) {
	int d=0;
	d += (c1.r - c2.r) * (c1.r - c2.r);
	d += (c1.g - c2.g) * (c1.g - c2.g);
	d += (c1.b - c2.b) * (c1.b - c2.b);
	return d;
}

void getKMeans(Image image, Color *palette, int* match, int k) {
	Color* colors = malloc(image.height * image.width * sizeof(Color));
	for (int i=0;i<image.height;i++) {
		unsigned char* src = image.data + (i * image.width)*3;
		for (int j=0;j<image.width;j++) {
			Color* c = colors + (i * image.width) + j;
			c->r = *src++;
			c->g = *src++;
			c->b = *src++;
		}
	}
	for (int i=0;i<k;i++) {
		palette[i] = colors[rand() % (image.width * image.height)];
	}
	int* newR = malloc(k * sizeof(int));
	int* newG = malloc(k * sizeof(int));
	int* newB = malloc(k * sizeof(int));
	int* count = malloc(k * sizeof(int));
	bool updated = true;
	while (updated) {
		updated = false;
		memset(newR,0,k*sizeof(int));
		memset(newG,0,k*sizeof(int));
		memset(newB,0,k*sizeof(int));
		memset(count,0,k*sizeof(int));
		for (size_t i=0;i<image.height*image.width;i++) {
			int minD=255*255*3,minID=0;
			for (int j=0;j<k;j++) {
				int d = getDiff(colors[i],palette[j]);
				if (d < minD) {
					minD = d;
					minID = j;
				}
			}
			newR[minID] += colors[i].r;
			newG[minID] += colors[i].g;
			newB[minID] += colors[i].b;
			count[minID]++;
			if (minID != match[i]) {
				match[i] = minID;
				updated = true;
			}
		}
		for (int j=0;j<k;j++) {
			if (count[j] == 0) {
				palette[j] = colors[rand() % (image.height*image.width)];
			} else {
				palette[j].r = newR[j] / count[j];
				palette[j].g = newG[j] / count[j];
				palette[j].b = newB[j] / count[j];
			}
		}
	}
	unsigned char* dst = image.data;
	for (int i=0;i<image.height * image.width;i++) {
		*dst++ = palette[match[i]].r;
		*dst++ = palette[match[i]].g;
		*dst++ = palette[match[i]].b;
	}
	free(newR);
	free(newG);
	free(newB);
	free(count);
	free(colors);
}

void recursive_fill(Image image, int x, int y, double* center_x, double* center_y, int* cluster_size) {
	if (x < 0 || x >= image.width || y < 0 || y >= image.height)
		return;
	if (!mask[y*image.width + x])
		return;
	mask[y*image.width + x] = false;
	*center_x += x;
	*center_y += y;
	*cluster_size += 1;
	recursive_fill(image, x-1, y, center_x, center_y, cluster_size);
	recursive_fill(image, x+1, y, center_x, center_y, cluster_size);
	recursive_fill(image, x, y-1, center_x, center_y, cluster_size);
	recursive_fill(image, x, y+1, center_x, center_y, cluster_size);
}

void filterBW(Image image, Image rgb_image, unsigned char l_threshold, unsigned char ab_threshold, int cluster_threshold) {
	unsigned char* src = image.data;
	for (int i=0;i<image.height*image.width;i++) {
		if (src[0] > l_threshold && src[1] > 120-ab_threshold && src[1] < 120+ab_threshold && src[2] > 136-ab_threshold && src[2] < 136+ab_threshold) {
			src[0] = 0xFF;
			src[1] = 0xFF;
			src[2] = 0xFF;
			mask[i] = true;
		} else {
			src[0] = 0;
			src[1] = 0;
			src[2] = 0;
			mask[i] = false;
		}
		src += 3;
	}

	double center_x;
	double center_y;
	int cluster_size=0;
	for (int i=0;i<image.height;i++) {
		for (int j=0;j<image.width;j++) {
			if (!mask[i*image.width + j])
				continue;
			center_x = 0;
			center_y = 0;
			cluster_size = 0;
			recursive_fill(image, j,i, &center_x, &center_y, &cluster_size);
//			printf("Frame %d cluster_size %d\n", frame_number, cluster_size);
			if (cluster_size > cluster_threshold)
				break;
		}
		if (cluster_size > cluster_threshold)
			break;
	}

	center_x /= cluster_size;
	center_y /= cluster_size;

#if USE_3D
	drawPendulum(center_x, center_y);
#else
	printf("Frame %d (%d: %6.2f,%6.2f)\r",frame_number, cluster_size, center_x, center_y);
	fflush(stdout);
	int box_size = 10;
	int box_x1 = center_x-box_size < 0 ? 0 : center_x-box_size >= image.width ? image.width-1 : center_x-box_size;
	int box_x2 = center_x+box_size < 0 ? 0 : center_x+box_size >= image.width ? image.width-1 : center_x+box_size;
	int box_y1 = center_y-box_size < 0 ? 0 : center_y-box_size >= image.height ? image.height-1 : center_y-box_size;
	int box_y2 = center_y+box_size < 0 ? 0 : center_y+box_size >= image.height ? image.height-1 : center_y+box_size;
	for (int i=box_y1;i<=box_y2;i++) {
		unsigned char* dst = rgb_image.data + (i * image.width + box_x1) * 3;
		for (int j=box_x1;j<=box_x2;j++) {
			dst[0] = 0;
			dst[1] = 0;
			dst[2] = 255;
			dst += 3;
		}
	}
	SDL_Flip(screen);
#endif
}

static void process_image(const void *p, int size)
{
        frame_number++;
		const unsigned char* yuyv_data = p;
		if (tmp_data == NULL) {
			tmp_data = malloc(640*480*3);
			mask = malloc(640*480*sizeof(bool));
		}
#if USE_3D
		unsigned char* rgb_data = tmp_data;
		Image rgb_img = {640,480,rgb_data};
		Image img = {640,480,tmp_data};
#else
		unsigned char* rgb_data = screen->pixels;
		Image rgb_img = {640,480,screen->pixels};
		Image img = {640,480,tmp_data};
//		Image img = {640,480,screen->pixels};
#endif
		for (int i=0;i<640*480/2;i++) {
			yuv2rgb(yuyv_data[0],yuyv_data[1],yuyv_data[3],rgb_data+2,rgb_data+1,rgb_data+0);
			yuv2rgb(yuyv_data[2],yuyv_data[1],yuyv_data[3],rgb_data+5,rgb_data+4,rgb_data+3);
			yuyv_data += 4;
			rgb_data += 6;
		}

#if !USE_3D
		memcpy(tmp_data, screen->pixels, img.width*img.height*3);
#endif

//		int K = 4;
//		int* match = malloc(img.width*img.height*sizeof(int));
//		Color* palette = malloc(K * sizeof(Color));
		convertLAB(img);
//		getKMeans(img,palette,match,K);
		filterBW(img, rgb_img, 250, 2, 2000);
//		free(match);
//		free(palette);

//		printf("Frame %d\r",frame_number);
//		fflush(stdout);
//        char filename[15];
//        sprintf(filename, "frame-%d.raw", frame_number);
//        FILE *fp=fopen(filename,"wb");
//        
//        if (out_buf)
//                fwrite(p, size, 1, fp);
//
//        fflush(fp);
//        fclose(fp);
}

static int read_frame(void)
{
        struct v4l2_buffer buf;
        unsigned int i;

        switch (io) {
        case IO_METHOD_READ:
                if (-1 == read(fd, buffers[0].start, buffers[0].length)) {
                        switch (errno) {
                        case EAGAIN:
                                return 0;

                        case EIO:
                                /* Could ignore EIO, see spec. */

                                /* fall through */

                        default:
                                errno_exit("read");
                        }
                }

                process_image(buffers[0].start, buffers[0].length);
                break;

        case IO_METHOD_MMAP:
                CLEAR(buf);

                buf.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
                buf.memory = V4L2_MEMORY_MMAP;

                if (-1 == xioctl(fd, VIDIOC_DQBUF, &buf)) {
                        switch (errno) {
                        case EAGAIN:
                                return 0;

                        case EIO:
                                /* Could ignore EIO, see spec. */

                                /* fall through */

                        default:
                                errno_exit("VIDIOC_DQBUF");
                        }
                }

                assert(buf.index < n_buffers);

                process_image(buffers[buf.index].start, buf.bytesused);

                if (-1 == xioctl(fd, VIDIOC_QBUF, &buf))
                        errno_exit("VIDIOC_QBUF");
                break;

        case IO_METHOD_USERPTR:
                CLEAR(buf);

                buf.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
                buf.memory = V4L2_MEMORY_USERPTR;

                if (-1 == xioctl(fd, VIDIOC_DQBUF, &buf)) {
                        switch (errno) {
                        case EAGAIN:
                                return 0;

                        case EIO:
                                /* Could ignore EIO, see spec. */

                                /* fall through */

                        default:
                                errno_exit("VIDIOC_DQBUF");
                        }
                }

                for (i = 0; i < n_buffers; ++i)
                        if (buf.m.userptr == (unsigned long)buffers[i].start
                            && buf.length == buffers[i].length)
                                break;

                assert(i < n_buffers);

                process_image((void *)buf.m.userptr, buf.bytesused);

                if (-1 == xioctl(fd, VIDIOC_QBUF, &buf))
                        errno_exit("VIDIOC_QBUF");
                break;
        }

        return 1;
}

static void mainloop(void)
{
//        unsigned int count;
//
//        count = frame_count;

//        while (count-- > 0) {
        while (1) {
                for (;;) {
                        fd_set fds;
                        struct timeval tv;
                        int r;

                        FD_ZERO(&fds);
                        FD_SET(fd, &fds);

                        /* Timeout. */
                        tv.tv_sec = 2;
                        tv.tv_usec = 0;

                        r = select(fd + 1, &fds, NULL, NULL, &tv);

                        if (-1 == r) {
                                if (EINTR == errno)
                                        continue;
                                errno_exit("select");
                        }

                        if (0 == r) {
                                fprintf(stderr, "select timeout\n");
                                exit(EXIT_FAILURE);
                        }

                        if (read_frame())
                                break;
                        /* EAGAIN - continue select loop. */
                }
        }
}

static void stop_capturing(void)
{
        enum v4l2_buf_type type;

        switch (io) {
        case IO_METHOD_READ:
                /* Nothing to do. */
                break;

        case IO_METHOD_MMAP:
        case IO_METHOD_USERPTR:
                type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
                if (-1 == xioctl(fd, VIDIOC_STREAMOFF, &type))
                        errno_exit("VIDIOC_STREAMOFF");
                break;
        }
}

static void start_capturing(void)
{
        unsigned int i;
        enum v4l2_buf_type type;

        switch (io) {
        case IO_METHOD_READ:
                /* Nothing to do. */
                break;

        case IO_METHOD_MMAP:
                for (i = 0; i < n_buffers; ++i) {
                        struct v4l2_buffer buf;

                        CLEAR(buf);
                        buf.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
                        buf.memory = V4L2_MEMORY_MMAP;
                        buf.index = i;

                        if (-1 == xioctl(fd, VIDIOC_QBUF, &buf))
                                errno_exit("VIDIOC_QBUF");
                }
                type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
                if (-1 == xioctl(fd, VIDIOC_STREAMON, &type))
                        errno_exit("VIDIOC_STREAMON");
                break;

        case IO_METHOD_USERPTR:
                for (i = 0; i < n_buffers; ++i) {
                        struct v4l2_buffer buf;

                        CLEAR(buf);
                        buf.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
                        buf.memory = V4L2_MEMORY_USERPTR;
                        buf.index = i;
                        buf.m.userptr = (unsigned long)buffers[i].start;
                        buf.length = buffers[i].length;

                        if (-1 == xioctl(fd, VIDIOC_QBUF, &buf))
                                errno_exit("VIDIOC_QBUF");
                }
                type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
                if (-1 == xioctl(fd, VIDIOC_STREAMON, &type))
                        errno_exit("VIDIOC_STREAMON");
                break;
        }
}

static void uninit_device(void)
{
        unsigned int i;

        switch (io) {
        case IO_METHOD_READ:
                free(buffers[0].start);
                break;

        case IO_METHOD_MMAP:
                for (i = 0; i < n_buffers; ++i)
                        if (-1 == munmap(buffers[i].start, buffers[i].length))
                                errno_exit("munmap");
                break;

        case IO_METHOD_USERPTR:
                for (i = 0; i < n_buffers; ++i)
                        free(buffers[i].start);
                break;
        }

        free(buffers);
}

static void init_read(unsigned int buffer_size)
{
        buffers = calloc(1, sizeof(*buffers));

        if (!buffers) {
                fprintf(stderr, "Out of memory\n");
                exit(EXIT_FAILURE);
        }

        buffers[0].length = buffer_size;
        buffers[0].start = malloc(buffer_size);

        if (!buffers[0].start) {
                fprintf(stderr, "Out of memory\n");
                exit(EXIT_FAILURE);
        }
}

static void init_mmap(void)
{
        struct v4l2_requestbuffers req;

        CLEAR(req);

        req.count = 4;
        req.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
        req.memory = V4L2_MEMORY_MMAP;

        if (-1 == xioctl(fd, VIDIOC_REQBUFS, &req)) {
                if (EINVAL == errno) {
                        fprintf(stderr, "%s does not support "
                                 "memory mapping\n", dev_name);
                        exit(EXIT_FAILURE);
                } else {
                        errno_exit("VIDIOC_REQBUFS");
                }
        }

        if (req.count < 2) {
                fprintf(stderr, "Insufficient buffer memory on %s\n",
                         dev_name);
                exit(EXIT_FAILURE);
        }

        buffers = calloc(req.count, sizeof(*buffers));

        if (!buffers) {
                fprintf(stderr, "Out of memory\n");
                exit(EXIT_FAILURE);
        }

        for (n_buffers = 0; n_buffers < req.count; ++n_buffers) {
                struct v4l2_buffer buf;

                CLEAR(buf);

                buf.type        = V4L2_BUF_TYPE_VIDEO_CAPTURE;
                buf.memory      = V4L2_MEMORY_MMAP;
                buf.index       = n_buffers;

                if (-1 == xioctl(fd, VIDIOC_QUERYBUF, &buf))
                        errno_exit("VIDIOC_QUERYBUF");

                buffers[n_buffers].length = buf.length;
                buffers[n_buffers].start =
                        mmap(NULL /* start anywhere */,
                              buf.length,
                              PROT_READ | PROT_WRITE /* required */,
                              MAP_SHARED /* recommended */,
                              fd, buf.m.offset);

                if (MAP_FAILED == buffers[n_buffers].start)
                        errno_exit("mmap");
        }
}

static void init_userp(unsigned int buffer_size)
{
        struct v4l2_requestbuffers req;

        CLEAR(req);

        req.count  = 4;
        req.type   = V4L2_BUF_TYPE_VIDEO_CAPTURE;
        req.memory = V4L2_MEMORY_USERPTR;

        if (-1 == xioctl(fd, VIDIOC_REQBUFS, &req)) {
                if (EINVAL == errno) {
                        fprintf(stderr, "%s does not support "
                                 "user pointer i/o\n", dev_name);
                        exit(EXIT_FAILURE);
                } else {
                        errno_exit("VIDIOC_REQBUFS");
                }
        }

        buffers = calloc(4, sizeof(*buffers));

        if (!buffers) {
                fprintf(stderr, "Out of memory\n");
                exit(EXIT_FAILURE);
        }

        for (n_buffers = 0; n_buffers < 4; ++n_buffers) {
                buffers[n_buffers].length = buffer_size;
                buffers[n_buffers].start = malloc(buffer_size);

                if (!buffers[n_buffers].start) {
                        fprintf(stderr, "Out of memory\n");
                        exit(EXIT_FAILURE);
                }
        }
}

static void init_device(void)
{
        struct v4l2_capability cap;
        struct v4l2_cropcap cropcap;
        struct v4l2_crop crop;
        struct v4l2_format fmt;
        unsigned int min;

        if (-1 == xioctl(fd, VIDIOC_QUERYCAP, &cap)) {
                if (EINVAL == errno) {
                        fprintf(stderr, "%s is no V4L2 device\n",
                                 dev_name);
                        exit(EXIT_FAILURE);
                } else {
                        errno_exit("VIDIOC_QUERYCAP");
                }
        }

        if (!(cap.capabilities & V4L2_CAP_VIDEO_CAPTURE)) {
                fprintf(stderr, "%s is no video capture device\n",
                         dev_name);
                exit(EXIT_FAILURE);
        }

        switch (io) {
        case IO_METHOD_READ:
                if (!(cap.capabilities & V4L2_CAP_READWRITE)) {
                        fprintf(stderr, "%s does not support read i/o\n",
                                 dev_name);
                        exit(EXIT_FAILURE);
                }
                break;

        case IO_METHOD_MMAP:
        case IO_METHOD_USERPTR:
                if (!(cap.capabilities & V4L2_CAP_STREAMING)) {
                        fprintf(stderr, "%s does not support streaming i/o\n",
                                 dev_name);
                        exit(EXIT_FAILURE);
                }
                break;
        }


        /* Select video input, video standard and tune here. */


        CLEAR(cropcap);

        cropcap.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;

        if (0 == xioctl(fd, VIDIOC_CROPCAP, &cropcap)) {
                crop.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
                crop.c = cropcap.defrect; /* reset to default */

                if (-1 == xioctl(fd, VIDIOC_S_CROP, &crop)) {
                        switch (errno) {
                        case EINVAL:
                                /* Cropping not supported. */
                                break;
                        default:
                                /* Errors ignored. */
                                break;
                        }
                }
        } else {
                /* Errors ignored. */
        }


        CLEAR(fmt);

        fmt.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
        if (force_format) {
//				fprintf(stderr, "Set H264\r\n");
				fprintf(stderr, "Set YUYV\r\n");
                fmt.fmt.pix.width       = 640; //replace
                fmt.fmt.pix.height      = 480; //replace
//                fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_H264; //replace
                fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_YUYV; //replace
                fmt.fmt.pix.field       = V4L2_FIELD_ANY;

                if (-1 == xioctl(fd, VIDIOC_S_FMT, &fmt))
                        errno_exit("VIDIOC_S_FMT");

                /* Note VIDIOC_S_FMT may change width and height. */
        } else {
                /* Preserve original settings as set by v4l2-ctl for example */
                if (-1 == xioctl(fd, VIDIOC_G_FMT, &fmt))
                        errno_exit("VIDIOC_G_FMT");
        }

        /* Buggy driver paranoia. */
        min = fmt.fmt.pix.width * 2;
        if (fmt.fmt.pix.bytesperline < min)
                fmt.fmt.pix.bytesperline = min;
        min = fmt.fmt.pix.bytesperline * fmt.fmt.pix.height;
        if (fmt.fmt.pix.sizeimage < min)
                fmt.fmt.pix.sizeimage = min;

        switch (io) {
        case IO_METHOD_READ:
                init_read(fmt.fmt.pix.sizeimage);
                break;

        case IO_METHOD_MMAP:
                init_mmap();
                break;

        case IO_METHOD_USERPTR:
                init_userp(fmt.fmt.pix.sizeimage);
                break;
        }
}

static void close_device(void)
{
        if (-1 == close(fd))
                errno_exit("close");

        fd = -1;
}

static void open_device(void)
{
        struct stat st;

        if (-1 == stat(dev_name, &st)) {
                fprintf(stderr, "Cannot identify '%s': %d, %s\n",
                         dev_name, errno, strerror(errno));
                exit(EXIT_FAILURE);
        }

        if (!S_ISCHR(st.st_mode)) {
                fprintf(stderr, "%s is no device\n", dev_name);
                exit(EXIT_FAILURE);
        }

        fd = open(dev_name, O_RDWR /* required */ | O_NONBLOCK, 0);

        if (-1 == fd) {
                fprintf(stderr, "Cannot open '%s': %d, %s\n",
                         dev_name, errno, strerror(errno));
                exit(EXIT_FAILURE);
        }
}

static void usage(FILE *fp, int argc, char **argv)
{
        fprintf(fp,
                 "Usage: %s [options]\n\n"
                 "Version 1.3\n"
                 "Options:\n"
                 "-d | --device name   Video device name [%s]\n"
                 "-h | --help          Print this message\n"
                 "-m | --mmap          Use memory mapped buffers [default]\n"
                 "-r | --read          Use read() calls\n"
                 "-u | --userp         Use application allocated buffers\n"
                 "-o | --output        Outputs stream to stdout\n"
                 "-f | --format        Force format to 640x480 YUYV\n"
                 "-c | --count         Number of frames to grab [%i]\n"
                 "",
                 argv[0], dev_name, frame_count);
}

static const char short_options[] = "d:hmruofc:";

static const struct option
long_options[] = {
        { "device", required_argument, NULL, 'd' },
        { "help",   no_argument,       NULL, 'h' },
        { "mmap",   no_argument,       NULL, 'm' },
        { "read",   no_argument,       NULL, 'r' },
        { "userp",  no_argument,       NULL, 'u' },
        { "output", no_argument,       NULL, 'o' },
        { "format", no_argument,       NULL, 'f' },
        { "count",  required_argument, NULL, 'c' },
        { 0, 0, 0, 0 }
};

int main(int argc, char **argv)
{
        dev_name = "/dev/video0";

        for (;;) {
                int idx;
                int c;

                c = getopt_long(argc, argv,
                                short_options, long_options, &idx);

                if (-1 == c)
                        break;

                switch (c) {
                case 0: /* getopt_long() flag */
                        break;

                case 'd':
                        dev_name = optarg;
                        break;

                case 'h':
                        usage(stdout, argc, argv);
                        exit(EXIT_SUCCESS);

                case 'm':
                        io = IO_METHOD_MMAP;
                        break;

                case 'r':
                        io = IO_METHOD_READ;
                        break;

                case 'u':
                        io = IO_METHOD_USERPTR;
                        break;

                case 'o':
                        out_buf++;
                        break;

                case 'f':
                        force_format++;
                        break;

                case 'c':
                        errno = 0;
                        frame_count = strtol(optarg, NULL, 0);
                        if (errno)
                                errno_exit(optarg);
                        break;

                default:
                        usage(stderr, argc, argv);
                        exit(EXIT_FAILURE);
                }
        }

		SDL_Init(SDL_INIT_VIDEO);
		SDL_WM_SetCaption("loadSwing", NULL);
#if USE_3D
		screen = SDL_SetVideoMode(800,600, 24, SDL_OPENGL);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(fov,(double)800/600,1,1000);
#else
		screen = SDL_SetVideoMode(640,480,24,SDL_SWSURFACE);
#endif

        open_device();
        init_device();
        start_capturing();
        mainloop();
        stop_capturing();
        uninit_device();
        close_device();
        fprintf(stderr, "\n");
        return 0;
}
