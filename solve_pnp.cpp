#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <opencv/cv.h>
#define DEBUG 0
#define EXTRINSIC_GUESS 0

struct Point {
	float x,y,z;
};

int main(int argc, char* argv[]) {
	if (argc<6) {
		printf("./solve_pnp width height target_point.txt depth_buffer.txt camera_location.txt [0.match ..]\n");
		return 1;
	}

	srand(time(NULL));
	int width = atoi(argv[1]);
	int height = atoi(argv[2]);
	char buffer[128];
	double fx=width,fy=width,cx=(width-1)/2,cy=(height-1)/2;
	std::vector<float> x_target,y_target;
	std::vector<int> match_target;
	std::vector<Point> depth_buffer;
	FILE* target_point = fopen(argv[3],"r");
	while (fgets(buffer,128,target_point)) {
		float x,y;
		int numMatch;
		if (sscanf(buffer,"%f %f %d",&x,&y,&numMatch)==3) {
			x_target.push_back(x);
			y_target.push_back(y);
			match_target.push_back(numMatch);
		}
	}
	fclose(target_point);

	FILE* depth_buffer_txt = fopen(argv[4],"r");
	while (fgets(buffer,128,depth_buffer_txt)) {
		Point p;
		if (sscanf(buffer,"%f %f %f",&p.x,&p.y,&p.z)==3)
			depth_buffer.push_back(p);
	}
	fclose(depth_buffer_txt);

	FILE* camera_location = fopen(argv[5],"w");
	int index = 0;
	cv::Mat rvec,tvec;
	while (true) {
		std::vector<cv::Point3f> objectPoints;
		std::vector<cv::Point2f> imagePoints;
		sprintf(buffer,"%d.site.match",index);
		FILE* key_match = fopen(buffer,"r");
		if (!key_match)
			break;
		while (fgets(buffer,128,key_match)) {
			int id1,id2;
			float x1,y1,x2,y2;
			if (sscanf(buffer,"%d %f %f %d %f %f",&id1,&x1,&y1,&id2,&x2,&y2)==6) {
				cv::Point2f p(x1,y1);
				Point q = depth_buffer[(int)y2 * width + (int)x2];
				cv::Point3f r(q.x,q.y,q.z);
				objectPoints.push_back(r);
				imagePoints.push_back(p);
			}
		}
		fclose(key_match);

		if (objectPoints.size() >= 4) {
			cv::Mat cameraMatrix = (cv::Mat_<double>(3,3) << fx, 0, cx, 0, fy, cy, 0, 0, 1);
			cv::Mat distCoeffs;
#if EXTRINSIC_GUESS
			cv::solvePnP(objectPoints,imagePoints,cameraMatrix,distCoeffs,rvec,tvec,index!=0);
#else
			cv::solvePnP(objectPoints,imagePoints,cameraMatrix,distCoeffs,rvec,tvec);
#endif
		}

		double theta = cv::norm(rvec);
		cv::Mat rvec_n = rvec / theta;
		double rx = rvec_n.at<double>(0,0);
		double ry = rvec_n.at<double>(1,0);
		double rz = rvec_n.at<double>(2,0);
		cv::Mat R = (cv::Mat_<double>(3,3) << 
			cos(theta) + (1-cos(theta))*rx*rx,
			(1-cos(theta))*rx*ry - sin(theta)*rz,
			(1-cos(theta))*rx*rz + sin(theta)*ry,
			(1-cos(theta))*rx*ry + sin(theta)*rz,
			cos(theta) + (1-cos(theta))*ry*ry,
			(1-cos(theta))*ry*rz - sin(theta)*rx,
			(1-cos(theta))*rx*rz - sin(theta)*ry,
			(1-cos(theta))*ry*rz + sin(theta)*rx,
			cos(theta) + (1-cos(theta))*rz*rz
			);
		cv::Mat R_inv = R.t();
		cv::Mat cameraPos = -R_inv * tvec;
		Point camera_angle;
		camera_angle.x = atan2(R_inv.at<double>(2,1),R_inv.at<double>(2,2));
		//camera_angle.y = atan2(-R_inv.at<double>(2,0),sqrt(R_inv.at<double>(2,1)*R_inv.at<double>(2,1)+R_inv.at<double>(2,2)*R_inv.at<double>(2,2)));
		camera_angle.y = -asin(R_inv.at<double>(2,0));
		camera_angle.z = atan2(R_inv.at<double>(1,0),R_inv.at<double>(0,0));
		cv::Mat viewVec = (cv::Mat_<double>(3,1) << (x_target[index]-cx)*fy, (y_target[index]-cy)*fx, fx*fy);
		viewVec = R_inv * viewVec;
		viewVec = viewVec / cv::norm(viewVec);
#if DEBUG
		printf("rvec: %8.2f %8.2f %8.2f\n",rvec.at<double>(0,0),rvec.at<double>(1,0),rvec.at<double>(2,0));
		printf("tvec: %8.2f %8.2f %8.2f\n",tvec.at<double>(0,0),tvec.at<double>(1,0),tvec.at<double>(2,0));
		printf("R: %8.2f %8.2f %8.2f\n",R.at<double>(0,0),R.at<double>(0,1),R.at<double>(0,2));
		printf("   %8.2f %8.2f %8.2f\n",R.at<double>(1,0),R.at<double>(1,1),R.at<double>(1,2));
		printf("   %8.2f %8.2f %8.2f\n",R.at<double>(2,0),R.at<double>(2,1),R.at<double>(2,2));
		printf("cameraPos: %8.2f %8.2f %8.2f\n",cameraPos.at<double>(0,0),cameraPos.at<double>(1,0),cameraPos.at<double>(2,0));
#endif
		if (objectPoints.size() >= 4 && cameraPos.at<double>(2,0) > 0) {
			double reproj_err = 0;
			for (size_t i=0;i<objectPoints.size();i++) {
				cv::Mat O = (cv::Mat_<double>(3,1) << objectPoints[i].x,objectPoints[i].y,objectPoints[i].z);
				cv::Mat m = R * O + tvec;
				cv::Point3f p(m.at<double>(0,0),m.at<double>(1,0),m.at<double>(2,0));
				cv::Point2f q = imagePoints[i];
				cv::Point2f r(p.x/p.z*fx+cx,p.y/p.z*fy+cy);
#if DEBUG
				printf("(%f %f %f) reprojection (%f %f) (%f %f)\n",
					objectPoints[i].x,objectPoints[i].y,objectPoints[i].z,
					r.x,r.y,
					q.x,q.y);
#endif
				reproj_err += (r.x - q.x) * (r.x - q.x);
				reproj_err += (r.y - q.y) * (r.y - q.y);
			}
			reproj_err = sqrt(reproj_err / objectPoints.size());
			printf("%d.match PnP: optimized %lu points (RMSE = %f)\n",index,objectPoints.size(),reproj_err);
			fprintf(camera_location,"%f %f %f %f %f %f %f %f %f %f %d\n",
				cameraPos.at<double>(0,0),cameraPos.at<double>(1,0),cameraPos.at<double>(2,0),
				camera_angle.x,camera_angle.y,camera_angle.z,
				viewVec.at<double>(0,0),viewVec.at<double>(1,0),viewVec.at<double>(2,0),
				reproj_err,match_target[index]);
		} else {
			printf("%d.match PnP: unable to optimize %lu points\n",index,objectPoints.size());
			fprintf(camera_location,"%f %f %f %f %f %f %f %f %f %f %d\n",
				cameraPos.at<double>(0,0),cameraPos.at<double>(1,0),cameraPos.at<double>(2,0),
				camera_angle.x,camera_angle.y,camera_angle.z,
				viewVec.at<double>(0,0),viewVec.at<double>(1,0),viewVec.at<double>(2,0),
				-1.0,match_target[index]);
		}
		index++;
	}
	fclose(camera_location);
	return 1;

}
