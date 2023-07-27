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
	Point2f position = Point2f(-1, -1); // �ҷ�������λ��
	Point2f uavSpeed = Point2f(-1, -1); // �ҷ��������ٶ�Point2f(�ٶȴ�С, �ٶȷ���);
	vector<double> distance; // Ŀ���״����
	vector<double> disTheta; // Ŀ���״﷽λ�� 0-360
	vector<Point2f> targetPosition; // Ŀ��ת����ͬһ����ϵ�������
	vector<double> speed; // Ŀ�꾶���ٶ�
						  /* �������� */
	vector<vector<double>> HRRP; // Ŀ��һά������
								 /* ... */
	uavData(Point2f position, Point2f uavSpeed)
	{
		this->position = position;
		this->uavSpeed = uavSpeed;
	}
};

/* 
��������: ����ҪĿ��Ϊ���ģ��ڵѿ�������ϵ�У�����Ŀ��������λ����Ϣ
vector<Point2f> targetLocation: ���Ŀ������λ��
vector<Point2f> targetSpeed: ��������ٶ�
double colsMidIndex: ��ҪĿ��ĺ�����
double rowsIndex: ��ҪĿ���������
*/
void targetInfoInit(vector<Point2f>& targetLocation, vector<Point2f>& targetSpeed, double colsMidIndex, double rowsIndex)
{
	// λ����Ϣ
	targetLocation.push_back(Point2f(colsMidIndex, rowsIndex)); //��ҪĿ��
	int angle1 = 60, dist1 = 500;
	targetLocation.push_back(Point2f(colsMidIndex - dist1, rowsIndex)); // ����ҪĿ��500��(60�ȵȷ��Ų�)  1
	targetLocation.push_back(Point2f(colsMidIndex - dist1 * cos(1.0 * angle1 * CV_PI / 180), rowsIndex - dist1 * sin(1.0 * angle1 * CV_PI / 180))); // 2
	targetLocation.push_back(Point2f(colsMidIndex + dist1 * cos(1.0 * angle1 * CV_PI / 180), rowsIndex - dist1 * sin(1.0 * angle1 * CV_PI / 180))); // 3
	targetLocation.push_back(Point2f(colsMidIndex + dist1, rowsIndex)); // 4
	//int angle2 = 45, dist2 = 1000; 
	//targetLocation.push_back(Point2f(colsMidIndex - dist2, rowsIndex)); // ���1ǧ��(45�ȵȷ��Ų�) 1
	//targetLocation.push_back(Point2f(colsMidIndex - dist2 * cos(1.0 * angle2 * CV_PI / 180), rowsIndex + dist2 * sin(1.0 * angle2 * CV_PI / 180))); // 2
	//targetLocation.push_back(Point2f(colsMidIndex - dist2 * cos(1.0 * angle2 * 2 * CV_PI / 180), rowsIndex + dist2 * sin(1.0 * angle2 * 2 * CV_PI / 180))); // 3
	//targetLocation.push_back(Point2f(colsMidIndex + dist2 * cos(1.0 * angle2 * CV_PI / 180), rowsIndex + dist2 * sin(1.0 * angle2 * CV_PI / 180))); // 4
	//targetLocation.push_back(Point2f(colsMidIndex + dist2, rowsIndex)); // 5
	//int angle3 = 36, dist3 = 2000;
	//targetLocation.push_back(Point2f(colsMidIndex - dist3, rowsIndex)); // ���2ǧ��(36�ȵȷ��Ų�) 1
	//targetLocation.push_back(Point2f(colsMidIndex - dist3 * cos(1.0 * angle3 * CV_PI / 180), rowsIndex + dist3 * sin(1.0 * angle3 * CV_PI / 180))); // 2
	//targetLocation.push_back(Point2f(colsMidIndex - dist3 * cos(1.0 * angle3 * 2 * CV_PI / 180), rowsIndex + dist3 * sin(1.0 * angle3 * 2 * CV_PI / 180))); // 3
	//targetLocation.push_back(Point2f(colsMidIndex + dist3 * cos(1.0 * angle3 * 2 * CV_PI / 180), rowsIndex + dist3 * sin(1.0 * angle3 * 2 * CV_PI / 180))); // 4
	//targetLocation.push_back(Point2f(colsMidIndex + dist3 * cos(1.0 * angle3 * CV_PI / 180), rowsIndex + dist3 * sin(1.0 * angle3 * CV_PI / 180))); // 5
	//targetLocation.push_back(Point2f(colsMidIndex + dist3, rowsIndex)); // 6			

	// �ٶ���Ϣ
	targetSpeed.push_back(Point(15, 270 * CV_PI / 180)); targetSpeed.push_back(Point(15, 270 * CV_PI / 180));
	targetSpeed.push_back(Point(15, 270 * CV_PI / 180)); targetSpeed.push_back(Point(15, 270 * CV_PI / 180));
	targetSpeed.push_back(Point(15, 270 * CV_PI / 180));
}

/*
��������: ����ҪĿ��Ϊ���ģ��ڵѿ�������ϵ�У������ҷ����˻�

double colsMidIndex: ��ҪĿ��ĺ�����
double rowsIndex: ��ҪĿ���������
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
��������: �������ɵ�Ŀ�����ꡢ�ٶȵ���Ϣ������ʵ���״��ȡ�ĺ��в���������Ϣ
vector<uavData>& uavDatas: �ҷ����˻�������
vector<Point2f> targetLocation: ���Ŀ������λ��
vector<Point2f> targetSpeed: ��������ٶ�
*/
void calculateRadarInfo(vector<uavData>& uavDatas, vector<Point2f> targetLocation, vector<Point2f> targetSpeed)
{
	for (int m = 0; m < uavDatas.size(); ++m)
	{
		random_device rd;
		mt19937 gen(rd());
		shuffle(targetLocation.begin(), targetLocation.end(), gen);
		int radarDisError = 30; // ����radarDisError�׵Ĳ�����
		int radarAngError = 0.2; // ����radarAngError�ȵĲ�����
		for (int i = 0; i < targetLocation.size(); ++i)
		{
			// �״����
			double ethe = ((static_cast<double>(rand()) / RAND_MAX) * radarAngError * 2 - radarAngError) * CV_PI / 180;
			double edis = (rand() % radarDisError * 2 - radarDisError) + rand() / (RAND_MAX * 1.f);
			/*double ethe = 0;
			double edis = 0;*/
			uavDatas[m].distance.push_back(getDistance(uavDatas[m].position, targetLocation[i]) + edis);// ���˻���Ŀ�����
			uavDatas[m].disTheta.push_back(getdisTheta(uavDatas[m].position, targetLocation[i]) + ethe);// ���˻���Ŀ��Ƕ�(�ѿ�������ϵ0-180)
		}
	}

}



// ��ʵ�ʵѿ�������ϵ�ĵ㣬��opencv��Mat��ʽͼ������ʾ
int drawPoint(Mat& input, vector<Point2f> vPoints)
{
	



	return 0;
}


int main()
{
	// 5000x5000���أ� 1���� = 1��
	Mat RealLocationImage = Mat::zeros(Size(5000, 5000), CV_8UC3);
	int colsMidIndex = RealLocationImage.cols / 2; // ��Ŀ��x��������
	int rowsIndex = RealLocationImage.rows - 100; // ��Ŀ��y��������
	int textSize = 2; // ��ͼ���ִ�С
	int pointSize = 10; // ��ͼĿ���С
	// ����Ŀ����Ϣ λ�� �ٶ�
	vector<Point2f> targetLocation;
	vector<Point2f> targetSpeed; // �ٶȶ���Ϊ15m/s,��y������
	targetInfoInit(targetLocation, targetSpeed, colsMidIndex, rowsIndex); // ����Ŀ����Ϣ
	// �ҷ����˻������ʼ��
	vector<uavData> uavDatas;
	uavInfoInit(uavDatas, colsMidIndex, rowsIndex);
	

	// ��ʾ��ͼ��
	//int kk = 0;
	//for (int i = 0; i < targetLocation.size(); ++i)
	//{
	//	circle(RealLocationImage, targetLocation[i], 10, Scalar(0, 255, 255), -1);
	//	putText(RealLocationImage, std::to_string(kk++), targetLocation[i] + Point2f(-5, 30),
	//		cv::FONT_HERSHEY_SIMPLEX, textSize, CV_RGB(255, 230, 255), 1.8);
	//}
		
	
	
	// �����״��ȡ����Ϣ
	calculateRadarInfo(uavDatas, targetLocation, targetSpeed);
	//for (int m = 0; m < uavDatas.size(); ++m)
	//{
	//	random_device rd;
	//	mt19937 gen(rd());
	//	shuffle(targetLocation.begin(), targetLocation.end(), gen);
	//	int radarDisError = 30; // ����radarDisError�׵Ĳ�����
	//	int radarAngError = 0.2; // ����radarAngError�ȵĲ�����
	//	for (int i = 0; i < targetLocation.size(); ++i)
	//	{
	//		// �״����
	//		double ethe = ((static_cast<double>(rand()) / RAND_MAX) * radarAngError * 2 - radarAngError) * CV_PI / 180;
	//		double edis = (rand() % radarDisError * 2 - radarDisError) + rand() / (RAND_MAX * 1.f);
	//		/*double ethe = 0;
	//		double edis = 0;*/
	//		uavDatas[m].distance.push_back(getDistance(uavDatas[m].position, targetLocation[i]) + edis);// ���˻���Ŀ�����
	//		uavDatas[m].disTheta.push_back(getdisTheta(uavDatas[m].position, targetLocation[i]) + ethe);// ���˻���Ŀ��Ƕ�(�ѿ�������ϵ0-180)
	//	}
	//}
	
	// �Է�����1��2��ʼλ���е�Ϊԭ�㣬��������ϵת��
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
			// ת��Ϊͬһ����ϵ
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

			// ��ͼ
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

// ��ȡ�����ľ���
double getDistance(Point2f p1, Point2f p2)
{
	double dx = (p1.x - p2.x)*(p1.x - p2.x);
	double dy = (p1.y - p2.y)*(p1.y - p2.y);
	return sqrt(dx + dy);
}

// ��ȡ����֮��ĽǶ�
double getdisTheta(Point2f uav, Point2f target)
{
	double dx = (uav.x - target.x);
	double dy = (uav.y - target.y);
	double angle = atan2(dy, dx); // -180��+180
	return angle; // ����
				  //return angle * 180 / CV_PI; // �Ƕ�
}