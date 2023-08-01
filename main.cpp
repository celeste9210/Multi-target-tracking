#include <iostream>
#include <ctime>
#include <algorithm>
#include <random>
#include <map>
// #include <opencv2/opencv.hpp>
#include "D:/Softwares/OpenCV320/build/include/opencv2/opencv.hpp"

using namespace std;
using namespace cv;
//rectangle(targetImage, Point2f(width/2-1,height-3*2), Point2f(width / 2-1, height), Scalar(0, 0, 255),-1);

double getdisTheta(Point2f uav, Point2f target);
double getDistance(Point2f p1, Point2f p2);
double getRSpeed(Point2f targetSpeed, double theta);

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

Scalar colors[] = { Scalar(0, 0, 255),    // 红色 0
Scalar(0, 255, 0),    // 绿色 1
Scalar(0, 128, 255),  // 橙色 2
Scalar(0, 255, 255),  // 黄色 3
Scalar(255, 255, 0),  // 青色 4
Scalar(255, 0, 255),  // 洋红色 5
Scalar(0, 0, 128),    // 深红色 6
Scalar(0, 128, 0),    // 深绿色 7
Scalar(255, 0, 128),  // 紫色 8
Scalar(128, 128, 0),  // 深黄色 9
Scalar(128, 0, 128),  // 深紫色 10
Scalar(128, 128, 0),  // 深青色 11
Scalar(128, 255, 0),  // 青绿色 12
Scalar(255, 0, 0),    // 深蓝色 13
Scalar(128, 128, 255), // 浅红色 14
Scalar(255, 0, 0),    // 蓝色 15
Scalar(128, 255, 128), // 浅绿色 16
Scalar(255, 128, 128), // 浅蓝色 17
Scalar(128, 255, 255), // 浅黄色 18
Scalar(255, 128, 255),  // 浅洋红色 19
Scalar(255, 255, 255)  // 白色 20

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
	int angle2 = 45, dist2 = 1000; 
	targetLocation.push_back(Point2f(colsMidIndex - dist2, rowsIndex)); // 相距1千米(45度等分排布) 1
	targetLocation.push_back(Point2f(colsMidIndex - dist2 * cos(1.0 * angle2 * CV_PI / 180), rowsIndex - dist2 * sin(1.0 * angle2 * CV_PI / 180))); // 2
	targetLocation.push_back(Point2f(colsMidIndex - dist2 * cos(1.0 * angle2 * 2 * CV_PI / 180), rowsIndex - dist2 * sin(1.0 * angle2 * 2 * CV_PI / 180))); // 3
	targetLocation.push_back(Point2f(colsMidIndex + dist2 * cos(1.0 * angle2 * CV_PI / 180), rowsIndex - dist2 * sin(1.0 * angle2 * CV_PI / 180))); // 4
	targetLocation.push_back(Point2f(colsMidIndex + dist2, rowsIndex)); // 5
	int angle3 = 36, dist3 = 2000;
	targetLocation.push_back(Point2f(colsMidIndex - dist3, rowsIndex)); // 相距2千米(36度等分排布) 1
	targetLocation.push_back(Point2f(colsMidIndex - dist3 * cos(1.0 * angle3 * CV_PI / 180), rowsIndex - dist3 * sin(1.0 * angle3 * CV_PI / 180))); // 2
	targetLocation.push_back(Point2f(colsMidIndex - dist3 * cos(1.0 * angle3 * 2 * CV_PI / 180), rowsIndex - dist3 * sin(1.0 * angle3 * 2 * CV_PI / 180))); // 3
	targetLocation.push_back(Point2f(colsMidIndex + dist3 * cos(1.0 * angle3 * 2 * CV_PI / 180), rowsIndex - dist3 * sin(1.0 * angle3 * 2 * CV_PI / 180))); // 4
	targetLocation.push_back(Point2f(colsMidIndex + dist3 * cos(1.0 * angle3 * CV_PI / 180), rowsIndex - dist3 * sin(1.0 * angle3 * CV_PI / 180))); // 5
	targetLocation.push_back(Point2f(colsMidIndex + dist3, rowsIndex)); // 6			

	// 速度信息
	targetSpeed.push_back(Point(15, 270 * CV_PI / 180)); targetSpeed.push_back(Point(15, 270 * CV_PI / 180));
	targetSpeed.push_back(Point(15, 270 * CV_PI / 180)); targetSpeed.push_back(Point(15, 270 * CV_PI / 180));
	targetSpeed.push_back(Point(15, 270 * CV_PI / 180));
}

/*
函数功能: 以主要目标为中心，在笛卡尔坐标系中，生成我方无人机坐标和速度
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
void calculateRadarInfo(vector<uavData>& uavDatas, vector<Point2f>& targetLocation, vector<Point2f>& targetSpeed)
{
	for (int m = 0; m < uavDatas.size(); ++m)
	{
		random_device rd;
		mt19937 gen(rd());
		shuffle(targetLocation.begin(), targetLocation.end(), gen);
		int radarDisError = 99; // 正负radarDisError米的测距误差
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
			//uavDatas[m].speed.push_back(getRSpeed(targetSpeed[i], getdisTheta(uavDatas[m].position, targetLocation[i]) + ethe));//目标径向速度
		}
	}

}

/*
函数功能: 根据飞行雷达信息，目标信息，将目标位置转化到同一坐标系下
Mat& sysLocationImage: 系统坐标系
vector<uavData>& uavDatas: 我方无人机数据
Point2f center: 统一坐标系原点在大地坐标系下的坐标
*/
void coordinatCvt(Mat& sysLocationImage, vector<uavData>& uavDatas, Point2f center)
{
	// 
	center = (uavDatas[0].position + uavDatas[1].position) / 2;
	//Mat sysLocationImage = Mat::zeros(Size(5000, 5000), CV_8UC3);
	int cc = 0;
	Point2f drawCenter(sysLocationImage.cols / 2, sysLocationImage.rows / 2);
	for (int i = 0; i < uavDatas[0].distance.size(); ++i) // 1
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
			/*switch (j)
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
			}*/
		}
		cc++;
	}
 	waitKey();
}

/*
函数功能: 把实际笛卡尔坐标系的点，在opencv的Mat格式图像上显示
Mat& input: 系统坐标系
vector<Point2f> vPoints: 要绘制的点的笛卡尔坐标系下的坐标
Scalar pointColor: 目标颜色
int textIndex：目标顺序
int pointSize: 绘制目标大小
int textSize: 文字大小
*/
int drawPoint(Mat& input, vector<Point2f> vPoints, Scalar pointColor, int textIndex, int pointSize = 10, int textSize = 2)
{
	for (int i = 0; i < vPoints.size(); ++i)
	{
		vPoints[i].y = input.rows - vPoints[i].y;
		circle(input, vPoints[i], pointSize, pointColor, -1);
		/*putText(input, std::to_string(textIndex), vPoints[i] + Point2f(-30, 0),
			cv::FONT_HERSHEY_SIMPLEX, textSize, CV_RGB(255, 255, 255), 2);
		textIndex++;*/
	}
	return 0;
}



/*计算欧式距离*/
float calcuDistance(double* ptr, double* ptrCen, int cols) {
	float d = 0.0;
	for (int j = 0; j < cols; j++)
	{
		d += (double)(ptr[j] - ptrCen[j])*(ptr[j] - ptrCen[j]);
	}
	d = sqrt(d);
	return d;
}


/** @brief   最大最小距离聚类
@param data  输入样本数据，每一行为一个样本，每个样本可以存在多个特征数据
@param Theta 阈值，一般设置为0.5，阈值越小聚类中心越多
@param stt 第一个聚类中心是哪个数据
@param centerIndex 聚类中心的下标
@return 返回每个样本的类别，类别从1开始，0表示未分类或者分类失败
*/
Mat  MaxMinDisFun(Mat data1, float Theta, int stt, vector<int>& centerIndex)
{
	centerIndex.clear();
	Mat data;
	// 归一化
	normalize(data1, data);
	double maxDistance = 0;
	int start = stt;    //初始选一个中心点
	int index = start; //相当于指针指示新中心点的位置
	int k = 0;        //中心点计数，也即是类别
	int dataNum = data.rows; //输入的样本数
	//vector<int>	centerIndex;// 记录聚类中心是第几个样本
	Mat distance = cv::Mat::zeros(cv::Size(1, dataNum), CV_32FC1); //表示所有样本到当前聚类中心的距离
	Mat minDistance = cv::Mat::zeros(cv::Size(1, dataNum), CV_32FC1); //取较小距离
	Mat classes = cv::Mat::zeros(cv::Size(1, dataNum), CV_32SC1);     //表示每个点属于哪个类别
	centerIndex.push_back(index); //保存第一个聚类中心

	for (int i = 0; i < dataNum; i++) 
	{
		double* ptr1 = data.ptr<double>(i);
		double* ptrCen = data.ptr<double>(centerIndex[0]);
		//int ta = *ptr1;
		//int ba = *ptrCen;
		float d = calcuDistance(ptr1, ptrCen, data.cols);
		distance.at<float>(i, 0) = d;
		classes.at<int>(i, 0) = k + 1; // 初始时全划为第1类
		if (maxDistance < d)
		{
			maxDistance = d;
			index = i; //与第一个聚类中心距离最大的样本
		}// index:第二个聚类中心，maxDistance:第1-2个聚类中心之间的距离 distance:每个点和第一个聚类中心的距离
	}
	

	minDistance = distance.clone(); // minDistance:每个点和第一个(上一个)聚类中心的距离
	double minVal; double maxVal; Point minLoc; Point maxLoc;
	maxVal = maxDistance; // Z1-Z2距离
	while (maxVal > (maxDistance*Theta)) 
	{
		k = k + 1;
		centerIndex.push_back(index); //新的聚类中心
		for (int i = 0; i < dataNum; i++)
		{
			double* ptr1 = data.ptr<double>(i);
			double* ptrCen = data.ptr<double>(centerIndex[k]);
			float d = calcuDistance(ptr1, ptrCen, data.cols);
			distance.at<float>(i, 0) = d; //和当前聚类中心的距离
			//按照当前最近临方式分类，哪个近就分哪个类别
			// 到当前聚类中心的距离 < 到上一个聚类中心的距离，则划分为当前聚类中心
			// 但应该考虑同一目标位置相近且不超过一定误差范围
			if (distance.at<float>(i, 0) < minDistance.at<float>(i, 0)/* && distance.at<float>(i, 0) < 200*/)
			{
				minDistance.at<float>(i, 0) = distance.at<float>(i, 0);
				classes.at<int>(i, 0) = k + 1;
			}
		} 
		//查找minDistance中最大值
		cv::minMaxLoc(minDistance, &minVal, &maxVal, &minLoc, &maxLoc);
		index = maxLoc.y;
	}
	return classes;
}
int findKmax(vector<uavData> uavDatas)
{
	int Kmax = uavDatas[0].targetPosition.size();
	for (int i = 1; i < uavDatas.size(); ++i)
	{
		if (uavDatas[i].targetPosition.size()>Kmax)
		{
			Kmax = uavDatas[i].targetPosition.size();
		}
	}
	return Kmax;
}

int findKmin(vector<uavData> uavDatas)
{
	int Kmin = uavDatas[0].targetPosition.size();
	for (int i = 1; i < uavDatas.size(); ++i)
	{
		if (uavDatas[i].targetPosition.size()<Kmin)
		{
			Kmin = uavDatas[i].targetPosition.size();
		}
	}
	return Kmin;
}

double getS2(Mat data, Mat tempClasses, vector<int> centerIndex)
{
	vector<Point2f> Kpoint; // 存放聚类中心的坐标
	vector<double> disMeans(centerIndex.size(), 0); // 存放每一类的平均距离
	// 统计每一类的聚类中心坐标
	for (int i = 0; i < centerIndex.size(); ++i)
	{
		Kpoint.push_back(Point2f(data.at<double>(centerIndex[i], 0), data.at<double>(centerIndex[i], 1)));
	}
	

	for (int i = 0; i < tempClasses.rows; ++i) // 遍历每一条数据
	{
		int nClass = tempClasses.at<int>(i, 0) - 1; // 该数据是哪一类
		Point2f nP(data.at<double>(i, 0), data.at<double>(i, 1)); // 该点的坐标
		double dis = getDistance(nP, Kpoint[nClass]); // 计算该点到对应聚类中心的距离
		disMeans[nClass] += dis;
	}
	// 每条数据对应的类转换为vector
	vector<int> Classes(tempClasses.begin<int>(), tempClasses.end<int>());
	for (int i = 0; i < centerIndex.size(); ++i)
	{
		disMeans[i] = disMeans[i] / count(Classes.begin(), Classes.end(), i+1);
	}
	double disMean = 0;
	for (int i = 0; i < disMeans.size(); ++i)
	{
		disMean += disMeans[i];
	}
	disMean = disMean / disMeans.size();
	double variance = 0.0;
	for (double value : disMeans) 
	{
		double diff = value - disMean;
		variance += diff * diff;
	}
	variance /= disMeans.size();
	return variance;
}

/*多属性关联*/
void MAA(vector<uavData> uavDatas,Mat& mImage)
{
	Mat data = Mat::zeros(uavDatas[0].targetPosition.size()*uavDatas.size(), 2, CV_64FC1);
	// 将每个无人机的目标位置数据导入聚类矩阵
	int m = 0;
	for (int i = 0; i < uavDatas.size(); ++i)
	{
		for (int j = 0; j < uavDatas[i].targetPosition.size(); ++j)
		{
			int n = 0;
			data.at<double>(m, n++) = uavDatas[i].targetPosition[j].x;
			data.at<double>(m, n) = uavDatas[i].targetPosition[j].y;
			m++;
		}
	}
	// 多属性关联(聚类)
	vector<int> centerIndex; // 存放聚类中心是第几条数据
	float Theta = 0.2; // 该参数应根据雷达误差和目标最近距离综合设定。越小则簇的数量越多,越大则簇越少
	Mat classes;
	//Mat tempClasses = MaxMinDisFun(data, Theta, 0, centerIndex); // 聚类
	//classes = tempClasses;
	int kmax = findKmax(uavDatas);
	int kmin = findKmin(uavDatas);
	for (int stt = 0; stt < data.rows*data.cols; ++stt)
	{
		Mat tempClasses = MaxMinDisFun(data, Theta, stt, centerIndex); // 聚类
		if (centerIndex.size() < kmin || centerIndex.size() > kmax)
			continue;
		// 方差
		double s2 = 0;
		s2 = getS2(data, tempClasses, centerIndex);
		if (s2 > 0.8)
			classes = tempClasses;
		else
			continue;
	/*	classes = tempClasses;*/
		// 绘图
		map<int, int> cnt;
		classes.convertTo(classes, CV_8UC1);
		// 绘制聚类中心(白色)
		for (int i = 0; i < centerIndex.size(); ++i)
		{
			int cc = centerIndex[i];
			Point2f tg;
			tg.x = data.at<double>(cc, 0);
			tg.y = mImage.rows - data.at<double>(cc, 1);
			circle(mImage, tg, 10, colors[20], -1);
		}
		// 同一类目标绘制同一标签
		int k = 0;
		uchar* pcls = classes.data;
		for (int i = 0; i < uavDatas.size(); ++i)
		{
			for (int j = 0; j < uavDatas[i].targetPosition.size(); ++j)
			{
				uavDatas[i].targetPosition[j].y = mImage.rows - uavDatas[i].targetPosition[j].y;
				putText(mImage, std::to_string(pcls[k++]), uavDatas[i].targetPosition[j] + Point2f(-30, 0),
					cv::FONT_HERSHEY_SIMPLEX, 2, CV_RGB(255, 255, 255), 2);
				//circle(mImage, uavDatas[i].targetPosition[j], 10, colors[i], -1);
				if (!cnt.count(pcls[k-1])) // 第一次出现(即聚类中心)画白色
				{
					cnt.insert({ pcls[k - 1], 0 });
					cnt[pcls[k - 1]]++;
					circle(mImage, uavDatas[i].targetPosition[j], 10, colors[20], -1);
				}
				else
				{
					circle(mImage, uavDatas[i].targetPosition[j], 10, colors[i], -1);
					cnt[pcls[k - 1]]++;
				}	
			}
		}
		centerIndex.clear();
		mImage = Mat::zeros(mImage.rows, mImage.cols, mImage.type());
	}
	
	waitKey();
		
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
	// 根据我方无人机坐标和目标信息，计算我方无人机测到的雷达信息(含雷达误差)
	calculateRadarInfo(uavDatas, targetLocation, targetSpeed);
	// 将无人机获得的目标位置，转化到同一坐标系
	Point2f center = (uavDatas[0].position + uavDatas[1].position) / 2;
	Mat sysLocationImage = Mat::zeros(Size(5000, 5000), CV_8UC3);
	coordinatCvt(sysLocationImage, uavDatas, center);
	// 画出目标在大地坐标系位置
	drawPoint(RealLocationImage, targetLocation, Scalar(0, 255, 255), 0);
	// 画出目标在系统坐标系的位置
	/*for (int i = 0; i < uavDatas.size(); ++i)
	{
		int textIndex = 0;
		drawPoint(sysLocationImage, uavDatas[i].targetPosition, colors[i], 10);
	}*/
	// 初始状态生成完成
	int ntime = 60; // 60秒
	for (int i = 0; i < ntime; ++i)
	{
		// 多属性关联
		MAA(uavDatas, sysLocationImage);
		
		// 状态更新
	}


	
	waitKey();
	return 0;
}

// 笛卡尔坐标系 获取两点间的距离
inline double getDistance(Point2f p1, Point2f p2)
{
	double dx = (p1.x - p2.x)*(p1.x - p2.x);
	double dy = (p1.y - p2.y)*(p1.y - p2.y);
	return sqrt(dx + dy);
}

// 笛卡尔坐标系 获取两点之间的角度
inline double getdisTheta(Point2f uav, Point2f target)
{
	double dx = (uav.x - target.x);
	double dy = (uav.y - target.y);
	double angle = atan2(dy, dx); // -180到+180
	return angle; // 弧度
	//return angle * 180 / CV_PI; // 角度
}

// 笛卡尔坐标系 获得目标径向速度（向我方行驶）
double getRSpeed(Point2f targetSpeed, double theta)
{
	double speed = 0;
	double spdTheta = CV_PI - targetSpeed.y + theta;
	speed = abs(targetSpeed.x * cos(spdTheta));
	if (spdTheta <= CV_PI / 2)
		return speed;// 靠近无人机
	else
		return -speed;// 远离无人机
}
