#include <iostream>
#include <ctime>
#include <algorithm>
#include <random>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;
//rectangle(targetImage, Point2f(width/2-1,height-3*2), Point2f(width / 2-1, height), Scalar(0, 0, 255),-1);

double getdisTheta(Point2f uav, Point2f target);
double getDistance(Point2f p1, Point2f p2);

class uavData
{
public:
	Point2f position = Point2f(-1, -1); // 我方飞行器位置
	Point2f uavSpeed = Point2f(-1, -1); // 我方飞行器速度Point2f(速度大小, 速度方向);
	vector<double> distance; // 目标雷达距离
	vector<double> disTheta; // 目标雷达方位角 0-360
	vector<Point2f> targetPosition; // 目标转化到同一坐标系后的坐标
	vector<double> speed; // 目标径向速度
						  /* 其他数据 */
	vector<vector<double>> HRRP; // 目标一维距离像
								 /* ... */
	uavData(Point2f position, Point2f uavSpeed)
	{
		this->position = position;
		this->uavSpeed = uavSpeed;
	}
};

/* 
函数功能: 以主要目标为中心，在笛卡尔坐标系中，生成目标的坐标和位置信息
vector<Point2f> targetLocation: 存放目标坐标位置
vector<Point2f> targetSpeed: 存放坐标速度
double colsMidIndex: 主要目标的横坐标
double rowsIndex: 主要目标的纵坐标
*/
void targetInfoInit(vector<Point2f>& targetLocation, vector<Point2f>& targetSpeed, double colsMidIndex, double rowsIndex)
{
	// 位置信息
	targetLocation.push_back(Point2f(colsMidIndex, rowsIndex)); //主要目标
	int angle1 = 60, dist1 = 500;
	targetLocation.push_back(Point2f(colsMidIndex - dist1, rowsIndex)); // 距主要目标500米(60度等分排布)  1
	targetLocation.push_back(Point2f(colsMidIndex - dist1 * cos(1.0 * angle1 * CV_PI / 180), rowsIndex - dist1 * sin(1.0 * angle1 * CV_PI / 180))); // 2
	targetLocation.push_back(Point2f(colsMidIndex + dist1 * cos(1.0 * angle1 * CV_PI / 180), rowsIndex - dist1 * sin(1.0 * angle1 * CV_PI / 180))); // 3
	targetLocation.push_back(Point2f(colsMidIndex + dist1, rowsIndex)); // 4
	//int angle2 = 45, dist2 = 1000; 
	//targetLocation.push_back(Point2f(colsMidIndex - dist2, rowsIndex)); // 相距1千米(45度等分排布) 1
	//targetLocation.push_back(Point2f(colsMidIndex - dist2 * cos(1.0 * angle2 * CV_PI / 180), rowsIndex + dist2 * sin(1.0 * angle2 * CV_PI / 180))); // 2
	//targetLocation.push_back(Point2f(colsMidIndex - dist2 * cos(1.0 * angle2 * 2 * CV_PI / 180), rowsIndex + dist2 * sin(1.0 * angle2 * 2 * CV_PI / 180))); // 3
	//targetLocation.push_back(Point2f(colsMidIndex + dist2 * cos(1.0 * angle2 * CV_PI / 180), rowsIndex + dist2 * sin(1.0 * angle2 * CV_PI / 180))); // 4
	//targetLocation.push_back(Point2f(colsMidIndex + dist2, rowsIndex)); // 5
	//int angle3 = 36, dist3 = 2000;
	//targetLocation.push_back(Point2f(colsMidIndex - dist3, rowsIndex)); // 相距2千米(36度等分排布) 1
	//targetLocation.push_back(Point2f(colsMidIndex - dist3 * cos(1.0 * angle3 * CV_PI / 180), rowsIndex + dist3 * sin(1.0 * angle3 * CV_PI / 180))); // 2
	//targetLocation.push_back(Point2f(colsMidIndex - dist3 * cos(1.0 * angle3 * 2 * CV_PI / 180), rowsIndex + dist3 * sin(1.0 * angle3 * 2 * CV_PI / 180))); // 3
	//targetLocation.push_back(Point2f(colsMidIndex + dist3 * cos(1.0 * angle3 * 2 * CV_PI / 180), rowsIndex + dist3 * sin(1.0 * angle3 * 2 * CV_PI / 180))); // 4
	//targetLocation.push_back(Point2f(colsMidIndex + dist3 * cos(1.0 * angle3 * CV_PI / 180), rowsIndex + dist3 * sin(1.0 * angle3 * CV_PI / 180))); // 5
	//targetLocation.push_back(Point2f(colsMidIndex + dist3, rowsIndex)); // 6			

	// 速度信息
	targetSpeed.push_back(Point(15, 270 * CV_PI / 180)); targetSpeed.push_back(Point(15, 270 * CV_PI / 180));
	targetSpeed.push_back(Point(15, 270 * CV_PI / 180)); targetSpeed.push_back(Point(15, 270 * CV_PI / 180));
	targetSpeed.push_back(Point(15, 270 * CV_PI / 180));
}

/*
函数功能: 以主要目标为中心，在笛卡尔坐标系中，生成我方无人机

double colsMidIndex: 主要目标的横坐标
double rowsIndex: 主要目标的纵坐标
*/
void uavInfoInit(vector<uavData>& uavDatas, double colsMidIndex, double rowsIndex)
{
	uavData uav1(Point2f(colsMidIndex - 5000, rowsIndex - 10000), Point2f(15, 90 * CV_PI / 180)), 
		uav2(Point2f(colsMidIndex + 5000, rowsIndex - 10000), Point2f(15, 90 * CV_PI / 180)),
		uav3(Point2f(colsMidIndex - 5000, rowsIndex - 20000), Point2f(15, 90 * CV_PI / 180)),
		uav4(Point2f(colsMidIndex + 5000, rowsIndex - 20000), Point2f(15, 90 * CV_PI / 180));

	uavDatas.push_back(uav1); uavDatas.push_back(uav2); uavDatas.push_back(uav3); uavDatas.push_back(uav4);
}

/*
函数功能: 根据生成的目标坐标、速度等信息，计算实际雷达获取的含有测量误差的信息
vector<uavData>& uavDatas: 我方无人机的数据
vector<Point2f> targetLocation: 存放目标坐标位置
vector<Point2f> targetSpeed: 存放坐标速度
*/
void calculateRadarInfo(vector<uavData>& uavDatas, vector<Point2f> targetLocation, vector<Point2f> targetSpeed)
{
	for (int m = 0; m < uavDatas.size(); ++m)
	{
		random_device rd;
		mt19937 gen(rd());
		shuffle(targetLocation.begin(), targetLocation.end(), gen);
		int radarDisError = 30; // 正负radarDisError米的测距误差
		int radarAngError = 0.2; // 左右radarAngError度的测角误差
		for (int i = 0; i < targetLocation.size(); ++i)
		{
			// 雷达误差
			double ethe = ((static_cast<double>(rand()) / RAND_MAX) * radarAngError * 2 - radarAngError) * CV_PI / 180;
			double edis = (rand() % radarDisError * 2 - radarDisError) + rand() / (RAND_MAX * 1.f);
			/*double ethe = 0;
			double edis = 0;*/
			uavDatas[m].distance.push_back(getDistance(uavDatas[m].position, targetLocation[i]) + edis);// 无人机与目标距离
			uavDatas[m].disTheta.push_back(getdisTheta(uavDatas[m].position, targetLocation[i]) + ethe);// 无人机与目标角度(笛卡尔坐标系0-180)
		}
	}

}



// 把实际笛卡尔坐标系的点，在opencv的Mat格式图像上显示
int drawPoint(Mat& input, vector<Point2f> vPoints)
{
	



	return 0;
}


int main()
{
	// 5000x5000像素， 1像素 = 1米
	Mat RealLocationImage = Mat::zeros(Size(5000, 5000), CV_8UC3);
	int colsMidIndex = RealLocationImage.cols / 2; // 主目标x方向坐标
	int rowsIndex = RealLocationImage.rows - 100; // 主目标y方向坐标
	int textSize = 2; // 绘图文字大小
	int pointSize = 10; // 绘图目标大小
	// 生成目标信息 位置 速度
	vector<Point2f> targetLocation;
	vector<Point2f> targetSpeed; // 速度都定为15m/s,沿y正方向
	targetInfoInit(targetLocation, targetSpeed, colsMidIndex, rowsIndex); // 生成目标信息
	// 我方无人机坐标初始化
	vector<uavData> uavDatas;
	uavInfoInit(uavDatas, colsMidIndex, rowsIndex);
	

	// 显示在图上
	//int kk = 0;
	//for (int i = 0; i < targetLocation.size(); ++i)
	//{
	//	circle(RealLocationImage, targetLocation[i], 10, Scalar(0, 255, 255), -1);
	//	putText(RealLocationImage, std::to_string(kk++), targetLocation[i] + Point2f(-5, 30),
	//		cv::FONT_HERSHEY_SIMPLEX, textSize, CV_RGB(255, 230, 255), 1.8);
	//}
		
	
	
	// 计算雷达获取的信息
	calculateRadarInfo(uavDatas, targetLocation, targetSpeed);
	//for (int m = 0; m < uavDatas.size(); ++m)
	//{
	//	random_device rd;
	//	mt19937 gen(rd());
	//	shuffle(targetLocation.begin(), targetLocation.end(), gen);
	//	int radarDisError = 30; // 正负radarDisError米的测距误差
	//	int radarAngError = 0.2; // 左右radarAngError度的测角误差
	//	for (int i = 0; i < targetLocation.size(); ++i)
	//	{
	//		// 雷达误差
	//		double ethe = ((static_cast<double>(rand()) / RAND_MAX) * radarAngError * 2 - radarAngError) * CV_PI / 180;
	//		double edis = (rand() % radarDisError * 2 - radarDisError) + rand() / (RAND_MAX * 1.f);
	//		/*double ethe = 0;
	//		double edis = 0;*/
	//		uavDatas[m].distance.push_back(getDistance(uavDatas[m].position, targetLocation[i]) + edis);// 无人机与目标距离
	//		uavDatas[m].disTheta.push_back(getdisTheta(uavDatas[m].position, targetLocation[i]) + ethe);// 无人机与目标角度(笛卡尔坐标系0-180)
	//	}
	//}
	
	// 以飞行器1和2初始位置中点为原点，进行坐标系转换
	Point2f center = (uavDatas[0].position + uavDatas[1].position) / 2;
	Mat sysLocationImage = Mat::zeros(Size(5000, 5000), CV_8UC3);
	int cc = 0;
	Point2f drawCenter(sysLocationImage.cols / 2, sysLocationImage.rows / 2);
	for (int i = 0; i < targetLocation.size(); ++i) // 1
	{
		for (int j = 0; j < uavDatas.size(); ++j)
		{
			double dx = uavDatas[j].distance[i] * cos(uavDatas[j].disTheta[i]);
			double dy = uavDatas[j].distance[i] * sin(uavDatas[j].disTheta[i]);
			// 转化为同一坐标系
			if (uavDatas[j].position.x - center.x <= 0)
				dx += center.x - uavDatas[j].position.x;
			else
				dx -= uavDatas[j].position.x - center.x;
			if (uavDatas[j].position.y - center.y <= 0)
				dy += center.y - uavDatas[j].position.y;
			else
				dy -= uavDatas[j].position.y - center.y;
			Point2f dpt = Point2f(dx, -dy);
			Point2f targ = center + dpt;
			uavDatas[j].targetPosition.push_back(targ);

			// 绘图
			switch (j)
			{
			case 0:
				circle(sysLocationImage, targ, 10, Scalar(255, 123, 255), -1);
				putText(sysLocationImage, std::to_string(cc), targ + Point2f(0, 30),
					cv::FONT_HERSHEY_SIMPLEX, textSize, CV_RGB(255, 123, 255), 1.8);
				break;
			case 1:
				circle(sysLocationImage, targ, 10, Scalar(0, 255, 0), -1);
				putText(sysLocationImage, std::to_string(cc), targ + Point2f(0, -15),
					cv::FONT_HERSHEY_SIMPLEX, textSize, CV_RGB(0, 255, 0), 1.8);
				break;
			case 2:
				circle(sysLocationImage, targ, 10, Scalar(0, 0, 255), -1);
				putText(sysLocationImage, std::to_string(cc), targ + Point2f(-30, 0),
					cv::FONT_HERSHEY_SIMPLEX, textSize, CV_RGB(255, 0, 0), 1.8);
				break;
			case 3:
				circle(sysLocationImage, targ, 10, Scalar(0, 255, 255), -1);
				putText(sysLocationImage, std::to_string(cc), targ + Point2f(-30, 0),
					cv::FONT_HERSHEY_SIMPLEX, textSize, CV_RGB(255, 0, 0), 1.8);
				break;
			}
		}
		cc++;
	}

	waitKey();
	return 0;
}

// 获取两点间的距离
double getDistance(Point2f p1, Point2f p2)
{
	double dx = (p1.x - p2.x)*(p1.x - p2.x);
	double dy = (p1.y - p2.y)*(p1.y - p2.y);
	return sqrt(dx + dy);
}

// 获取两点之间的角度
double getdisTheta(Point2f uav, Point2f target)
{
	double dx = (uav.x - target.x);
	double dy = (uav.y - target.y);
	double angle = atan2(dy, dx); // -180到+180
	return angle; // 弧度
				  //return angle * 180 / CV_PI; // 角度
}