//#include <iostream>
//#include <ctime>
//#include <random>
//#include <opencv2/opencv.hpp>
//
//using namespace std;
//using namespace cv;
//
//int main()
//{
//	double dx = 1.732;
//	double dy = 1;
//	double aa = atan2(-dy, -dx) * 180 / CV_PI;
//	cout << aa << endl;
//	waitKey();
//
//}
//
//#define  _CRT_SECURE_NO_WARNINGS
//#include <stdio.h>
//#include <stdlib.h>
//#include <opencv2/opencv.hpp>
//#include <iostream>
////#include "CoordinateConvertAPI.h"
//#define NODENUM 8
//#define TGTNUM 10
//using namespace cv;
//using namespace std;
//
//double uniform_rand(double lowerBndr, double upperBndr)
//{
//	return lowerBndr + ((double)rand() / (RAND_MAX + 1.0)) * (upperBndr - lowerBndr);
//}
////para:平均值，方差， return：误差
//double gauss_rand(double mean, double sigma)
//{
//	double x, y, r2;
//	do {
//		x = -1.0 + 2.0 * uniform_rand(0.0, 1.0);
//		y = -1.0 + 2.0 * uniform_rand(0.0, 1.0);
//		r2 = x * x + y * y;
//	} while (r2 > 1.0 || r2 == 0.0);
//	return mean + sigma * y * sqrt(-2.0 * log(r2) / r2);
//}
//
//
//void LLH2XYZ_test()
//{
//	//1 地理坐标转SX(OSG)小场景坐标系:(lon,lat,height)->(x,y,z)
//	//测试用例说明
//	//红军飞机坐标(地理坐标)
//	double srctgtpos_red[3];
//	srctgtpos_red[0] = 120.7239478;  //lon
//	srctgtpos_red[1] = 23.05812005;  //lat
//	srctgtpos_red[2] = 8000.0;       //height
//									 //蓝军飞机坐标(地理坐标)
//	double srctgtpos_blue[3];
//	srctgtpos_blue[0] = 121.3004868;  //lon
//	srctgtpos_blue[1] = 22.75562431;  //lat
//	srctgtpos_blue[2] = 8000.0;      //height
//
//									 //物方地心坐标下的速度
//	double srctgtspd[3] = { 0.0,0.0,0.0 };
//	//红军sx小场景坐标
//	double despos_red[3] = { 0,0,0 };
//	//蓝军sx小场景坐标
//	double despos_blue[3] = { 0,0,0 };
//	//物方OSG小场景速度
//	double desspd[3] = { 0,0,0 };
//	//场景地图基准地理坐标(台海全岛)中心
//	double mappos[3] = { 120.926735,23.864941,8000.0 };
//	//场景地图基准坐标对应的SX(OSG)小场景坐标点
//	double maposgpos[3] = { -19550.37,24292.38,8000.0 };
//	//计算红军对应的SX(OSG)小场景坐标点
//	//cc_geoToSXOSG(srctgtpos_red, 1, srctgtspd, mappos, maposgpos, despos_red, desspd);
//	//计算蓝军对应的SX(OSG)小场景坐标点
//	//cc_geoToSXOSG(srctgtpos_blue, 1, srctgtspd, mappos, maposgpos, despos_blue, desspd);
//	//计算红蓝军的在SX(OSG)小场景中的距离
//	double disRedToBlue = sqrt(abs(despos_red[0] - despos_blue[0]) * abs(despos_red[0] - despos_blue[0])
//		+ abs(despos_red[1] - despos_blue[1]) * abs(despos_red[1] - despos_blue[1]));
//	cout << despos_red[0] << "," << despos_red[1] << "," << despos_red[2] << endl;
//	cout << despos_blue[0] << "," << despos_blue[1] << "," << despos_blue[2] << endl;
//	cout << "dist:" << disRedToBlue << endl;
//
//
//	//2 SX(OSG)场景坐标转地理坐标:(x,y,z)->(lon,lat,height)
//	//测试用例说明
//	//红军飞机坐标（SX(OSG)场景坐标）
//	double srctgtpos_red2[3];
//	srctgtpos_red2[0] = despos_red[0];
//	srctgtpos_red2[1] = despos_red[1];
//	srctgtpos_red2[2] = despos_red[2];
//	//蓝军飞机坐标（SX(OSG)场景坐标）
//	double srctgtpos_blue2[3];
//	srctgtpos_blue2[0] = despos_blue[0];
//	srctgtpos_blue2[1] = despos_blue[1];
//	srctgtpos_blue2[2] = despos_blue[2];
//	//物方SX（OSG）场景下的速度
//	double srcTgtSpd_OSG[3] = { 0.0,0.0,0.0 };
//	//物方地理坐标系下速度
//	double desTgtSpd_GEO[3] = { 0.0,0.0,0.0 };
//	//红军地理坐标系下坐标
//	double despos_red_GEO[3] = { 0,0,0 };
//	//蓝军地理坐标系下坐标
//	double despos_blue_GEO[3] = { 0,0,0 };
//	//计算红军对应的地理坐标系下坐标
//	/*cc_SXOSGToGeo(true, srctgtpos_red2, srcTgtSpd_OSG, mappos, maposgpos, despos_red_GEO, desTgtSpd_GEO);*/
//	//计算蓝军对应的地理坐标系下坐标
//	/*cc_SXOSGToGeo(true, srctgtpos_blue2, srcTgtSpd_OSG, mappos, maposgpos, despos_blue_GEO, desTgtSpd_GEO);*/
//	std::cout << despos_red_GEO[0] << " ";
//	std::cout << despos_red_GEO[1] << " ";
//	std::cout << despos_red_GEO[2] << endl;
//	std::cout << despos_blue_GEO[0] << " ";
//	std::cout << despos_blue_GEO[1] << " ";
//	std::cout << despos_blue_GEO[2] << endl;
//
//	cin.get();
//
//}
//
//double Dist(Point2d a, Point2d b) {
//	return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
//}
//
//double posAddErr(Point2d& pos, double sigma)
//{
//	double disErr = gauss_rand(0, sigma);
//	double absDisErr = fabs(disErr);
//	Point2d errDir;
//	errDir.x = -1 + 2.0 * rand() / (RAND_MAX);
//	errDir.y = -1 + 2.0 * rand() / (RAND_MAX);
//	errDir /= sqrt(errDir.dot(errDir));
//
//	pos += absDisErr * errDir;
//	return disErr;
//}
//
//struct cmp
//{
//	bool operator()(const Point2d& pt1, const Point2d& pt2)
//	{
//		return pt1.y < pt2.y;
//	}
//};
//
//bool cmp2(const Point2d& pt1, const Point2d& pt2)
//{
//	return pt1.x < pt2.x;
//}
//
//void mysort(Point2d a[], int nCnt) {
//	Point2d baseVec = Point2d(1, 0);
//	vector<double> sins;
//	for (int i = 1; i < nCnt; ++i) {
//		Point2d vec = a[i] - a[0];
//		sins.push_back(vec.cross(baseVec) / (sqrt(vec.dot(vec))));
//	}
//	for (int i = 0; i < sins.size(); ++i) {
//		for (int j = i + 1; j < sins.size(); ++j) {
//			if (sins[i] > sins[j]) {
//				swap(sins[i], sins[j]);
//				swap(a[j + 1], a[i + 1]);
//			}
//		}
//	}
//}
////单个目标（飞行器位置，目标位置，目标方向（弧度），目标距离，飞行器数量，是否叠加误差）
//void tgtDir_Dis(Point2d* airPos, Point2d tgtpos, double *dir, double *dis, int nCnt, bool adderr = true) {
//	vector<Point2d> tgtDir;
//	vector<double>temp_bta;
//	vector<double>bta;
//	for (int i = 0; i < nCnt; ++i) {
//		tgtDir.push_back(tgtpos - airPos[i]);
//		temp_bta.push_back(fabs(tgtDir[i].cross(Point2d(1, 0))) / sqrt(tgtDir[i].dot(tgtDir[i])));//sin值
//		if (tgtDir[i].y > 0 && tgtDir[i].x > 0)//第一象限
//			bta.push_back(asin(temp_bta[i]));
//		else if (tgtDir[i].y > 0 && tgtDir[i].x < 0)//第二象限
//			bta.push_back(CV_PI - asin(temp_bta[i]));
//		else if (tgtDir[i].y < 0 && tgtDir[i].x < 0)//第三象限
//			bta.push_back(CV_PI + asin(temp_bta[i]));
//		else
//			bta.push_back(-asin(temp_bta[i]));//第四象限
//		dir[i] = bta[i];
//		if (adderr) {
//			//协同探测信息叠加误差
//			dir[i] += (-0.25 + 0.5 * rand() / RAND_MAX) * CV_PI / 180;
//			dis[i] += -15 + 30.0 * rand() / RAND_MAX;
//		}
//	}
//}
//
///***
//* @brief       : 多飞行器多目标_角度遍历（飞行器位置，目标位置，各目标相对各飞行器的绝对方向，各目标相对各飞行器的距离，飞行器数量，目标数量，是否叠加误差）
//* @in param    : Point2d* airPos, Point2d *tgtpos，double dis[][TGTNUM]，int nCnt, int tgtNum, bool adderr = true
//* @out param   : double dir[][TGTNUM], double dis[][TGTNUM]
//* @return      : none
//*/
//
//void tgtsDirDis(Point2d* airPos, Point2d *tgtpos, double dir[][TGTNUM], double dis[][TGTNUM], int nCnt, int tgtNum, bool adderr = true) {
//	vector<vector<Point2d>> tgtDirs;
//	vector<vector<double>>temp_btas;
//	vector<vector<double>>btas;
//	for (int i = 0; i < nCnt; ++i) {
//		vector<Point2d> tgtDir;
//		vector<double>temp_bta;
//		vector<double>bta;
//
//		for (int j = 0; j < tgtNum; ++j) {
//			tgtDir.push_back(tgtpos[j] - airPos[i]);
//			temp_bta.push_back(fabs(tgtDir[j].cross(Point2d(1, 0))) / sqrt(tgtDir[j].dot(tgtDir[j])));//sin值
//			if (tgtDir[j].y > 0 && tgtDir[j].x > 0)//第一象限
//				bta.push_back(asin(temp_bta[j]));
//			else if (tgtDir[j].y > 0 && tgtDir[j].x < 0)//第二象限
//				bta.push_back(CV_PI - asin(temp_bta[j]));
//			else if (tgtDir[j].y < 0 && tgtDir[j].x < 0)//第三象限
//				bta.push_back(CV_PI + asin(temp_bta[j]));
//			else
//				bta.push_back(-asin(temp_bta[j]));//第四象限
//			dir[i][j] = bta[j];
//			if (adderr) {
//				//协同探测信息叠加误差
//				dir[i][j] += (-0.25 + 0.5 * rand() / RAND_MAX) * CV_PI / 180;
//				dis[i][j] += -15 + 30.0 * rand() / RAND_MAX;
//			}
//
//		}
//	}
//}
//
//void relpos();
//void relpos2();
//void rotTest();
//
//int main()
//{
//	//relpos2();
//	rotTest();
//	waitKey(0);
//	return 0;
//}
//
//void relpos() {
//	srand((unsigned)time(0));
//	int iter_num = 1000;
//	const int nodeNum = 8;//节点数量
//	ofstream out("rloutput" + to_string(nodeNum) + "_err.csv");
//	out << "id" << "," << "真值x(m)" << "," << "真值y(m)" << "," << "惯导x(m)" << "," << "惯导y(m)" << ","
//		<< "结果x(m)" << "," << "结果y(m)" << "," << "定位精度提升(倍)" << endl;
//	int sideNum = 0;
//	double sum_imp = 0;
//	for (int i = 0; i < nodeNum; ++i) {
//		sideNum += i;
//	}
//
//	//蒙特卡洛模拟1000次
//	for (int iter = 0; iter < iter_num; ++iter) {
//
//		Point2d realpos[nodeNum];
//		vector<double> d;
//		int cr_num = 0;
//		//随机生成节点真值
//		while (1) {
//			if (cr_num > nodeNum - 1)
//				break;
//			realpos[cr_num] = Point2d(40000.0 * rand() / RAND_MAX - 20000, 40000.0 * rand() / RAND_MAX - 20000);
//			bool ismeet = true;
//			for (int i = 0; i < cr_num; ++i) {
//				double dis = Dist(realpos[cr_num], realpos[i]);
//				if (dis < 10000) {
//					ismeet = false;
//					break;
//				}
//			}
//			if (!ismeet)
//				continue;
//			else
//				cr_num++;
//		}
//		//对集群节点排序
//		sort(realpos, realpos + nodeNum, cmp2);
//		mysort(realpos, nodeNum);
//		Point2d errpos[nodeNum];
//		for (int i = 0; i < nodeNum; ++i) {
//			errpos[i] = realpos[i];
//			posAddErr(errpos[i], 1500);//叠加误差
//		}
//		for (int i = 0; i < nodeNum; ++i) {
//			for (int j = i + 1; j < nodeNum; ++j) {
//				d.push_back(Dist(realpos[i], realpos[j]));//根据真值计算测距距离
//			}
//		}
//		Point2d tgtpos;
//		double Dis1 = 0;
//		double Dis2 = 0;
//		while (1) {
//			tgtpos = Point2d(200000.0 * rand() / RAND_MAX - 100000, 200000.0 * rand() / RAND_MAX - 100000);
//			Dis1 = Dist(tgtpos, realpos[0]);
//			Dis2 = Dist(tgtpos, realpos[1]);
//			if (Dis1 > 80000 && Dis2 > 80000) {
//				break;
//			}
//		}
//		//目标向量
//		Point2d vec_tgtToplt1 = tgtpos - realpos[0];
//		Point2d vec_tgtToplt2 = tgtpos - realpos[1];
//		double bta1 = 0;
//		double bta2 = 0;
//		double test_cross = vec_tgtToplt1.cross(Point2d(1, 0));
//		double temp_bta1 = fabs(vec_tgtToplt1.cross(Point2d(1, 0))) / sqrt(vec_tgtToplt1.dot(vec_tgtToplt1));//sin值
//																											 //double temp_bta11 = vec_tgtToplt1.dot(Point2d(1, 0)) / sqrt(vec_tgtToplt1.dot(vec_tgtToplt1));//cos值
//
//		double temp_bta2 = fabs(vec_tgtToplt2.cross(Point2d(1, 0))) / sqrt(vec_tgtToplt2.dot(vec_tgtToplt2));//sin值
//																											 //double temp_bta22 = vec_tgtToplt2.dot(Point2d(1, 0)) / sqrt(vec_tgtToplt2.dot(vec_tgtToplt2));//cos值
//		if (vec_tgtToplt1.y > 0 && vec_tgtToplt1.x > 0)//第一象限
//			bta1 = asin(temp_bta1);
//		else if (vec_tgtToplt1.y > 0 && vec_tgtToplt1.x < 0)//第二象限
//			bta1 = CV_PI - asin(temp_bta1);
//		else if (vec_tgtToplt1.y < 0 && vec_tgtToplt1.x < 0)//第三象限
//			bta1 = CV_PI + asin(temp_bta1);
//		else
//			bta1 = -asin(temp_bta1);//第四象限
//
//		if (vec_tgtToplt2.y > 0 && vec_tgtToplt2.x > 0)//第一象限
//			bta2 = asin(temp_bta2);
//		else if (vec_tgtToplt2.y > 0 && vec_tgtToplt2.x < 0)//第二象限
//			bta2 = CV_PI - asin(temp_bta2);
//		else if (vec_tgtToplt2.y < 0 && vec_tgtToplt2.x < 0)//第三象限
//			bta2 = CV_PI + asin(temp_bta2);
//		else
//			bta2 = -asin(temp_bta2);//第四象限
//
//									//协同探测信息叠加误差
//		bta1 += (-0.25 + 0.5 * rand() / RAND_MAX) * CV_PI / 180;
//		bta2 += (-0.25 + 0.5 * rand() / RAND_MAX) * CV_PI / 180;
//		Dis1 += -15 + 30.0 * rand() / RAND_MAX;
//		Dis2 += -15 + 30.0 * rand() / RAND_MAX;
//
//		double agl1 = bta1 / CV_PI * 180;
//		double agl2 = bta2 / CV_PI * 180;
//
//		double SINS = sin(bta1 - bta2) / d[0] * Dis2;
//		//添加协同探测误差后目标三角形过约束(ensure)
//		if (SINS > 1)
//			SINS = 1;
//		else if (SINS < -1)
//			SINS = -1;
//		double theta = 0;
//
//		double bta[nodeNum - 2] = { 0 };
//		double cosA[nodeNum - 2] = { 0 };
//
//		Point2d RESpoint[nodeNum], TgtPoint[nodeNum];
//		double Tgt_theta = 0;
//		//Mat board = Mat::zeros(1000, 1000, CV_8UC3);
//		Mat board = Mat(1000, 1000, CV_8UC3, Scalar(128, 128, 128));
//		//目标函数值
//		double MinTarget = 1000000000000;
//		int SCALE = 200;
//
//		for (int yjx = 0; yjx < 2; yjx++) {
//			if (yjx == 0)
//				theta = asin(SINS) + bta1;
//			else
//				theta = 2 * bta1 - theta - CV_PI;
//
//			double agl_theta = theta / CV_PI * 180;
//
//			for (int jx = 0; jx < 2; jx++) {
//				if (jx == 0)
//					for (int i = 0; i < nodeNum - 2; ++i) {
//						cosA[i] = (d[0] * d[0] + d[i + 1] * d[i + 1] - d[nodeNum - 1 + i] * d[nodeNum - 1 + i]) / (2 * d[0] * d[i + 1]);
//						bta[i] = acos(cosA[i]);
//					}
//				else
//					for (int i = 0; i < nodeNum - 2; ++i) {
//						cosA[i] = (d[0] * d[0] + d[i + 1] * d[i + 1] - d[nodeNum - 1 + i] * d[nodeNum - 1 + i]) / (2 * d[0] * d[i + 1]);
//						bta[i] = -acos(cosA[i]);
//					}
//				//对角度的遍历计算会比较慢，可采用迭代法
//				Mat display = Mat::zeros(1000, 1000, CV_8UC3);
//
//				RESpoint[0].x = -d[0] * cos(theta) + errpos[0].x + errpos[1].x;
//				RESpoint[0].y = -d[0] * sin(theta) + errpos[0].y + errpos[1].y;
//
//				for (int i = 0; i < nodeNum - 2; ++i) {
//					RESpoint[0].x += errpos[i + 2].x - d[i + 1] * cos(bta[i] + theta);
//					RESpoint[0].y += errpos[i + 2].y - d[i + 1] * sin(bta[i] + theta);
//				}
//				RESpoint[0].x /= nodeNum;
//				RESpoint[0].y /= nodeNum;
//
//				//allpoint[0] = Point2d(0, 0);
//				RESpoint[1] = Point2d(RESpoint[0].x + d[0] * cos(theta), RESpoint[0].y + d[0] * sin(theta));
//				double L = sqrt(pow(RESpoint[0].x - errpos[0].x, 2) + pow(RESpoint[0].y - errpos[0].y, 2)) +
//					sqrt(pow(RESpoint[1].x - errpos[1].x, 2) + pow(RESpoint[1].y - errpos[1].y, 2));
//				for (int i = 2; i < nodeNum; ++i) {
//					RESpoint[i] = Point2d(RESpoint[0].x + d[i - 1] * cos(bta[i - 2] + theta), RESpoint[0].y + d[i - 1] * sin(bta[i - 2] + theta));
//					//if (i >= 3) {
//					//	double dis1 = Dist(RESpoint[i], RESpoint[2]);
//					//	double errdis = fabs(dis1 - d[nodeNum - 1 + i - 2]);
//					//	if(errdis>1)
//					//		RESpoint[i] = Point2d(RESpoint[0].x + d[i - 1] * cos(bta[i - 2] - theta), RESpoint[0].y + d[i-1] * sin(theta-bta[i - 2] ));
//					//}
//					L += sqrt(pow(RESpoint[i].x - errpos[i].x, 2) + pow(RESpoint[i].y - errpos[i].y, 2));
//				}
//				//RESpoint[3] = Point2d(RESpoint[0].x + d[2] * cos(theta), RESpoint[0].y + d[2] * sin(theta));
//
//				if (L < MinTarget) {
//					MinTarget = L;
//					for (int i = 0; i < nodeNum; ++i) {
//						TgtPoint[i] = RESpoint[i];
//					}
//					Tgt_theta = theta * 180 / CV_PI;
//				}
//
//				//绘图
//				//for (int i = 0; i < nodeNum; i++) {
//				//	if ((i + 1) == nodeNum)
//				//		line(display, Point(RESpoint[i].x / SCALE + board.rows / 2, RESpoint[i].y / SCALE + board.cols / 2), Point(RESpoint[0].x / SCALE + board.rows / 2, RESpoint[0].y / SCALE + board.cols / 2),
//				//			Scalar(0, 255, 0), 2);
//				//	else
//				//		line(display, Point(RESpoint[i].x / SCALE + board.rows / 2, RESpoint[i].y / SCALE + board.cols / 2), Point(RESpoint[i + 1].x / SCALE + board.rows / 2, RESpoint[i + 1].y / SCALE + board.cols / 2),
//				//			Scalar(0, 255, 0), 2);
//				//	circle(display, Point(errpos[i].x / SCALE + board.rows / 2, errpos[i].y / SCALE + board.cols / 2), 5, Scalar(0, 0, 255), 5);
//				//	circle(display, Point(realpos[i].x / SCALE + board.rows / 2, realpos[i].y / SCALE + board.cols / 2), 5, Scalar(255, 255, 255), 5);
//				//	putText(display, to_string(i + 1), Point(errpos[i].x / SCALE + board.rows / 2 + 10, errpos[i].y / SCALE + board.cols / 2 + 10), FONT_HERSHEY_COMPLEX, 2, cv::Scalar(255, 0, 0), 3, 8, 0);
//				//	putText(display, to_string(i + 1), Point(RESpoint[i].x / SCALE + board.rows / 2, RESpoint[i].y / SCALE + board.cols / 2), FONT_HERSHEY_COMPLEX, 2, cv::Scalar(0, 255, 255), 3, 8, 0);
//
//				//}
//				//line(display, Point(tgtpos.x / SCALE + board.rows / 2, tgtpos.y / SCALE + board.cols / 2), Point(RESpoint[0].x / SCALE + board.rows / 2, RESpoint[0].y / SCALE + board.cols / 2),
//				//	Scalar(250, 255, 0), 2);
//				//line(display, Point(tgtpos.x / SCALE + board.rows / 2, tgtpos.y / SCALE + board.cols / 2), Point(RESpoint[1].x / SCALE + board.rows / 2, RESpoint[1].y / SCALE + board.cols / 2),
//				//	Scalar(250, 255, 0), 2);
//
//			}
//		}
//		////结果转为经纬高
//		//double TgtLLH[3][3] = { 0 };
//		double dis1[nodeNum] = { 0 };
//		double dis2[nodeNum] = { 0 };
//		double sum1 = 0;
//		double sum2 = 0;
//		for (int i = 0; i < nodeNum; i++) {
//			if ((i + 1) == nodeNum)
//				line(board, Point(TgtPoint[i].x / SCALE + board.rows / 2, TgtPoint[i].y / SCALE + board.cols / 2), Point(TgtPoint[0].x / SCALE + board.rows / 2, TgtPoint[0].y / SCALE + board.cols / 2),
//					Scalar(0, 255, 0), 2);
//			else
//				line(board, Point(TgtPoint[i].x / SCALE + board.rows / 2, TgtPoint[i].y / SCALE + board.cols / 2), Point(TgtPoint[i + 1].x / SCALE + board.rows / 2, TgtPoint[i + 1].y / SCALE + board.cols / 2),
//					Scalar(0, 255, 0), 2);
//			circle(board, Point(errpos[i].x / SCALE + board.rows / 2, errpos[i].y / SCALE + board.cols / 2), 5, Scalar(0, 0, 255), 5);
//			circle(board, Point(realpos[i].x / SCALE + board.rows / 2, realpos[i].y / SCALE + board.cols / 2), 5, Scalar(255, 255, 255), 5);
//
//			putText(board, to_string(i + 1), Point(errpos[i].x / SCALE + board.rows / 2 + 10, errpos[i].y / SCALE + board.cols / 2 + 10), FONT_HERSHEY_COMPLEX, 2, cv::Scalar(255, 0, 0), 3, 8, 0);
//			putText(board, to_string(i + 1), Point(TgtPoint[i].x / SCALE + board.rows / 2, TgtPoint[i].y / SCALE + board.cols / 2), FONT_HERSHEY_COMPLEX, 2, cv::Scalar(0, 255, 255), 3, 8, 0);
//			line(board, Point(tgtpos.x / SCALE + board.rows / 2, tgtpos.y / SCALE + board.cols / 2), Point(TgtPoint[0].x / SCALE + board.rows / 2, TgtPoint[0].y / SCALE + board.cols / 2),
//				Scalar(250, 255, 0), 2);
//			line(board, Point(tgtpos.x / SCALE + board.rows / 2, tgtpos.y / SCALE + board.cols / 2), Point(TgtPoint[1].x / SCALE + board.rows / 2, TgtPoint[1].y / SCALE + board.cols / 2),
//				Scalar(250, 255, 0), 2);
//
//			//double tgtPos[3] = { TgtPoint[i].x ,TgtPoint[i].y,10 };
//
//			//cc_SXOSGToGeo(true, tgtPos, srctgtspd, mapMid, mapgpMid, TgtLLH[i], desspd);
//			//cout << TgtLLH[i][0] << "," << TgtLLH[i][1] << endl;
//			dis1[i] = Dist(errpos[i], realpos[i]);
//			dis2[i] = Dist(TgtPoint[i], realpos[i]);
//
//			out << i + 1 << "," << realpos[i].x << "," << realpos[i].y << "," << errpos[i].x << "," << errpos[i].y
//				<< "," << TgtPoint[i].x << "," << TgtPoint[i].y << ","
//				<< dis1[i] / dis2[i] << endl;
//			sum1 += dis1[i] * dis1[i];
//			sum2 += dis2[i] * dis2[i];
//
//		}
//		double imp = sqrt(sum1 / nodeNum) / sqrt(sum2 / nodeNum);
//		if (imp < 1)
//			continue;
//		sum_imp += imp;
//		out << "精度总提升：" << "," << imp << "," << "飞行器1，2与x轴夹角：" << "," << Tgt_theta << endl;
//		//cout << imp<<endl;
//	}
//	out << "平均提升：" << "," << sum_imp / iter_num;
//	cout << "平均提升：" << sum_imp / iter_num;
//	out.close();
//
//
//}
////协同探测信息叠加误差后遍历
//void relpos2() {
//	srand((unsigned)time(0));
//	int iter_num = 500;
//	const int nodeNum = 8;//节点数量
//	ofstream out("rloutput" + to_string(nodeNum) + "_err.csv");
//	out << "id" << "," << "真值x(m)" << "," << "真值y(m)" << "," << "惯导x(m)" << "," << "惯导y(m)" << ","
//		<< "结果x(m)" << "," << "结果y(m)" << "," << "定位精度提升(倍)" << endl;
//	int sideNum = 0;
//	double sum_imp = 0;
//	for (int i = 0; i < nodeNum; ++i) {
//		sideNum += i;
//	}
//
//	//蒙特卡洛模拟
//	for (int iter = 0; iter < iter_num; ++iter) {
//
//		Point2d realpos[nodeNum];
//		vector<double> d;
//		int cr_num = 0;
//		//随机生成节点真值
//		while (1) {
//			if (cr_num > nodeNum - 1)
//				break;
//			realpos[cr_num] = Point2d(40000.0 * rand() / RAND_MAX - 20000, 40000.0 * rand() / RAND_MAX - 20000);
//			bool ismeet = true;
//			for (int i = 0; i < cr_num; ++i) {
//				double dis = Dist(realpos[cr_num], realpos[i]);
//				if (dis < 10000) {
//					ismeet = false;
//					break;
//				}
//			}
//			if (!ismeet)
//				continue;
//			else
//				cr_num++;
//		}
//		//对集群节点排序
//		sort(realpos, realpos + nodeNum, cmp2);
//		mysort(realpos, nodeNum);
//		Point2d errpos[nodeNum];
//		for (int i = 0; i < nodeNum; ++i) {
//			errpos[i] = realpos[i];
//			posAddErr(errpos[i], 1500);//叠加误差
//		}
//		for (int i = 0; i < nodeNum; ++i) {
//			for (int j = i + 1; j < nodeNum; ++j) {
//				d.push_back(Dist(realpos[i], realpos[j]));//根据真值计算测距距离
//			}
//		}
//		Point2d tgtpos;
//		double Dis1 = 0;
//		double Dis2 = 0;
//		while (1) {
//			tgtpos = Point2d(200000.0 * rand() / RAND_MAX - 100000, 200000.0 * rand() / RAND_MAX - 100000);
//			Dis1 = Dist(tgtpos, realpos[0]);
//			Dis2 = Dist(tgtpos, realpos[1]);
//			if (Dis1 > 80000 && Dis2 > 80000) {
//				break;
//			}
//		}
//		//目标向量
//		Point2d vec_tgtToplt1 = tgtpos - realpos[0];
//		Point2d vec_tgtToplt2 = tgtpos - realpos[1];
//		double bta1 = 0;
//		double bta2 = 0;
//		double test_cross = vec_tgtToplt1.cross(Point2d(1, 0));
//		double temp_bta1 = fabs(vec_tgtToplt1.cross(Point2d(1, 0))) / sqrt(vec_tgtToplt1.dot(vec_tgtToplt1));//sin值
//																											 //double temp_bta11 = vec_tgtToplt1.dot(Point2d(1, 0)) / sqrt(vec_tgtToplt1.dot(vec_tgtToplt1));//cos值
//
//		double temp_bta2 = fabs(vec_tgtToplt2.cross(Point2d(1, 0))) / sqrt(vec_tgtToplt2.dot(vec_tgtToplt2));//sin值
//																											 //double temp_bta22 = vec_tgtToplt2.dot(Point2d(1, 0)) / sqrt(vec_tgtToplt2.dot(vec_tgtToplt2));//cos值
//		if (vec_tgtToplt1.y > 0 && vec_tgtToplt1.x > 0)//第一象限
//			bta1 = asin(temp_bta1);
//		else if (vec_tgtToplt1.y > 0 && vec_tgtToplt1.x < 0)//第二象限
//			bta1 = CV_PI - asin(temp_bta1);
//		else if (vec_tgtToplt1.y < 0 && vec_tgtToplt1.x < 0)//第三象限
//			bta1 = CV_PI + asin(temp_bta1);
//		else
//			bta1 = -asin(temp_bta1);//第四象限
//
//		if (vec_tgtToplt2.y > 0 && vec_tgtToplt2.x > 0)//第一象限
//			bta2 = asin(temp_bta2);
//		else if (vec_tgtToplt2.y > 0 && vec_tgtToplt2.x < 0)//第二象限
//			bta2 = CV_PI - asin(temp_bta2);
//		else if (vec_tgtToplt2.y < 0 && vec_tgtToplt2.x < 0)//第三象限
//			bta2 = CV_PI + asin(temp_bta2);
//		else
//			bta2 = -asin(temp_bta2);//第四象限
//
//									//协同探测信息叠加误差
//		bta1 += (-0.25 + 0.5 * rand() / RAND_MAX) * CV_PI / 180;
//		bta2 += (-0.25 + 0.5 * rand() / RAND_MAX) * CV_PI / 180;
//		Dis1 += -15 + 30.0 * rand() / RAND_MAX;
//		Dis2 += -15 + 30.0 * rand() / RAND_MAX;
//
//		double agl1 = bta1 / CV_PI * 180;
//		double agl2 = bta2 / CV_PI * 180;
//
//		double SINS = sin(bta1 - bta2) / d[0] * Dis2;
//		//添加协同探测误差后目标三角形过约束(ensure)
//		if (SINS > 1)
//			SINS = 1;
//		else if (SINS < -1)
//			SINS = -1;
//		double temp_theta = 0;
//
//		double bta[nodeNum - 2] = { 0 };
//		double cosA[nodeNum - 2] = { 0 };
//
//		Point2d RESpoint[nodeNum], TgtPoint[nodeNum];
//		double Tgt_theta = 0;
//		//Mat board = Mat::zeros(1000, 1000, CV_8UC3);
//		Mat board = Mat(1000, 1000, CV_8UC3, Scalar(128, 128, 128));
//		//目标函数值
//		double MinTarget = 1000000000000;
//		int SCALE = 200;
//
//		for (int yjx = 0; yjx < 2; yjx++) {
//			if (yjx == 0)
//				temp_theta = asin(SINS) + bta1;
//			else
//				temp_theta = 2 * bta1 - temp_theta - CV_PI;
//
//
//			double agl_theta = temp_theta / CV_PI * 180;
//			double tep_theta = agl_theta - 0.5;
//			double theta = 0;
//			for (int cnSt = 0; cnSt < 50; ++cnSt) {
//				theta = (tep_theta + cnSt * 0.02) * CV_PI / 180;
//				for (int jx = 0; jx < 2; jx++) {
//					if (jx == 0)
//						for (int i = 0; i < nodeNum - 2; ++i) {
//							cosA[i] = (d[0] * d[0] + d[i + 1] * d[i + 1] - d[nodeNum - 1 + i] * d[nodeNum - 1 + i]) / (2 * d[0] * d[i + 1]);
//							bta[i] = acos(cosA[i]);
//						}
//					else
//						for (int i = 0; i < nodeNum - 2; ++i) {
//							cosA[i] = (d[0] * d[0] + d[i + 1] * d[i + 1] - d[nodeNum - 1 + i] * d[nodeNum - 1 + i]) / (2 * d[0] * d[i + 1]);
//							bta[i] = -acos(cosA[i]);
//						}
//					Mat display = Mat::zeros(1000, 1000, CV_8UC3);
//
//					RESpoint[0].x = -d[0] * cos(theta) + errpos[0].x + errpos[1].x;
//					RESpoint[0].y = -d[0] * sin(theta) + errpos[0].y + errpos[1].y;
//
//					for (int i = 0; i < nodeNum - 2; ++i) {
//						RESpoint[0].x += errpos[i + 2].x - d[i + 1] * cos(bta[i] + theta);
//						RESpoint[0].y += errpos[i + 2].y - d[i + 1] * sin(bta[i] + theta);
//					}
//					RESpoint[0].x /= nodeNum;
//					RESpoint[0].y /= nodeNum;
//
//					//allpoint[0] = Point2d(0, 0);
//					RESpoint[1] = Point2d(RESpoint[0].x + d[0] * cos(theta), RESpoint[0].y + d[0] * sin(theta));
//					double L = sqrt(pow(RESpoint[0].x - errpos[0].x, 2) + pow(RESpoint[0].y - errpos[0].y, 2)) +
//						sqrt(pow(RESpoint[1].x - errpos[1].x, 2) + pow(RESpoint[1].y - errpos[1].y, 2));
//					for (int i = 2; i < nodeNum; ++i) {
//						RESpoint[i] = Point2d(RESpoint[0].x + d[i - 1] * cos(bta[i - 2] + theta), RESpoint[0].y + d[i - 1] * sin(bta[i - 2] + theta));
//						//if (i >= 3) {
//						//	double dis1 = Dist(RESpoint[i], RESpoint[2]);
//						//	double errdis = fabs(dis1 - d[nodeNum - 1 + i - 2]);
//						//	if(errdis>1)
//						//		RESpoint[i] = Point2d(RESpoint[0].x + d[i - 1] * cos(bta[i - 2] - theta), RESpoint[0].y + d[i-1] * sin(theta-bta[i - 2] ));
//						//}
//						L += sqrt(pow(RESpoint[i].x - errpos[i].x, 2) + pow(RESpoint[i].y - errpos[i].y, 2));
//					}
//					//RESpoint[3] = Point2d(RESpoint[0].x + d[2] * cos(theta), RESpoint[0].y + d[2] * sin(theta));
//
//					if (L < MinTarget) {
//						MinTarget = L;
//						for (int i = 0; i < nodeNum; ++i) {
//							TgtPoint[i] = RESpoint[i];
//						}
//						Tgt_theta = theta * 180 / CV_PI;
//					}
//
//					//绘图
//					//for (int i = 0; i < nodeNum; i++) {
//					//	if ((i + 1) == nodeNum)
//					//		line(display, Point(RESpoint[i].x / SCALE + board.rows / 2, RESpoint[i].y / SCALE + board.cols / 2), Point(RESpoint[0].x / SCALE + board.rows / 2, RESpoint[0].y / SCALE + board.cols / 2),
//					//			Scalar(0, 255, 0), 2);
//					//	else
//					//		line(display, Point(RESpoint[i].x / SCALE + board.rows / 2, RESpoint[i].y / SCALE + board.cols / 2), Point(RESpoint[i + 1].x / SCALE + board.rows / 2, RESpoint[i + 1].y / SCALE + board.cols / 2),
//					//			Scalar(0, 255, 0), 2);
//					//	circle(display, Point(errpos[i].x / SCALE + board.rows / 2, errpos[i].y / SCALE + board.cols / 2), 5, Scalar(0, 0, 255), 5);
//					//	circle(display, Point(realpos[i].x / SCALE + board.rows / 2, realpos[i].y / SCALE + board.cols / 2), 5, Scalar(255, 255, 255), 5);
//					//	putText(display, to_string(i + 1), Point(errpos[i].x / SCALE + board.rows / 2 + 10, errpos[i].y / SCALE + board.cols / 2 + 10), FONT_HERSHEY_COMPLEX, 2, cv::Scalar(255, 0, 0), 3, 8, 0);
//					//	putText(display, to_string(i + 1), Point(RESpoint[i].x / SCALE + board.rows / 2, RESpoint[i].y / SCALE + board.cols / 2), FONT_HERSHEY_COMPLEX, 2, cv::Scalar(0, 255, 255), 3, 8, 0);
//
//					//}
//					//line(display, Point(tgtpos.x / SCALE + board.rows / 2, tgtpos.y / SCALE + board.cols / 2), Point(RESpoint[0].x / SCALE + board.rows / 2, RESpoint[0].y / SCALE + board.cols / 2),
//					//	Scalar(250, 255, 0), 2);
//					//line(display, Point(tgtpos.x / SCALE + board.rows / 2, tgtpos.y / SCALE + board.cols / 2), Point(RESpoint[1].x / SCALE + board.rows / 2, RESpoint[1].y / SCALE + board.cols / 2),
//					//	Scalar(250, 255, 0), 2);
//
//				}
//			}
//		}
//
//		////结果转为经纬高
//		//double TgtLLH[3][3] = { 0 };
//		double dis1[nodeNum] = { 0 };
//		double dis2[nodeNum] = { 0 };
//		double sum1 = 0;
//		double sum2 = 0;
//		for (int i = 0; i < nodeNum; i++) {
//			if ((i + 1) == nodeNum)
//				line(board, Point(TgtPoint[i].x / SCALE + board.rows / 2, TgtPoint[i].y / SCALE + board.cols / 2), Point(TgtPoint[0].x / SCALE + board.rows / 2, TgtPoint[0].y / SCALE + board.cols / 2),
//					Scalar(0, 255, 0), 2);
//			else
//				line(board, Point(TgtPoint[i].x / SCALE + board.rows / 2, TgtPoint[i].y / SCALE + board.cols / 2), Point(TgtPoint[i + 1].x / SCALE + board.rows / 2, TgtPoint[i + 1].y / SCALE + board.cols / 2),
//					Scalar(0, 255, 0), 2);
//			circle(board, Point(errpos[i].x / SCALE + board.rows / 2, errpos[i].y / SCALE + board.cols / 2), 5, Scalar(0, 0, 255), 5);
//			circle(board, Point(realpos[i].x / SCALE + board.rows / 2, realpos[i].y / SCALE + board.cols / 2), 5, Scalar(255, 255, 255), 5);
//
//			putText(board, to_string(i + 1), Point(errpos[i].x / SCALE + board.rows / 2 + 10, errpos[i].y / SCALE + board.cols / 2 + 10), FONT_HERSHEY_COMPLEX, 2, cv::Scalar(255, 0, 0), 3, 8, 0);
//			putText(board, to_string(i + 1), Point(TgtPoint[i].x / SCALE + board.rows / 2, TgtPoint[i].y / SCALE + board.cols / 2), FONT_HERSHEY_COMPLEX, 2, cv::Scalar(0, 255, 255), 3, 8, 0);
//			line(board, Point(tgtpos.x / SCALE + board.rows / 2, tgtpos.y / SCALE + board.cols / 2), Point(TgtPoint[0].x / SCALE + board.rows / 2, TgtPoint[0].y / SCALE + board.cols / 2),
//				Scalar(250, 255, 0), 2);
//			line(board, Point(tgtpos.x / SCALE + board.rows / 2, tgtpos.y / SCALE + board.cols / 2), Point(TgtPoint[1].x / SCALE + board.rows / 2, TgtPoint[1].y / SCALE + board.cols / 2),
//				Scalar(250, 255, 0), 2);
//
//			//double tgtPos[3] = { TgtPoint[i].x ,TgtPoint[i].y,10 };
//
//			//cc_SXOSGToGeo(true, tgtPos, srctgtspd, mapMid, mapgpMid, TgtLLH[i], desspd);
//			//cout << TgtLLH[i][0] << "," << TgtLLH[i][1] << endl;
//			dis1[i] = Dist(errpos[i], realpos[i]);
//			dis2[i] = Dist(TgtPoint[i], realpos[i]);
//
//			out << i + 1 << "," << realpos[i].x << "," << realpos[i].y << "," << errpos[i].x << "," << errpos[i].y
//				<< "," << TgtPoint[i].x << "," << TgtPoint[i].y << ","
//				<< dis1[i] / dis2[i] << endl;
//			sum1 += dis1[i] * dis1[i];
//			sum2 += dis2[i] * dis2[i];
//
//		}
//		double imp = sqrt(sum1 / nodeNum) / sqrt(sum2 / nodeNum);
//		//if (imp < 1)
//		//	continue;
//		sum_imp += imp;
//		out << "精度总提升：" << "," << imp << "," << "飞行器1，2与x轴夹角：" << "," << Tgt_theta << endl;
//		cout << iter << endl;
//	}
//	out << "平均提升：" << "," << sum_imp / iter_num;
//	cout << "平均提升：" << sum_imp / iter_num;
//	out.close();
//
//}
//
//void rotTest() {
//	srand((unsigned)time(0));
//	int iter_num = 500;
//	const int nodeNum = NODENUM;//节点数量
//	const int tgtNum = TGTNUM;
//	Point2d realpos[nodeNum];
//	vector<double> d;
//	int cr_num = 0;
//	//随机生成节点坐标真值
//	while (1) {
//		if (cr_num > nodeNum - 1)
//			break;
//		realpos[cr_num] = Point2d(40000.0 * rand() / RAND_MAX - 20000, 40000.0 * rand() / RAND_MAX - 20000);
//		bool ismeet = true;
//		for (int i = 0; i < cr_num; ++i) {
//			double dis = Dist(realpos[cr_num], realpos[i]);
//			if (dis < 10000) {
//				ismeet = false;
//				break;
//			}
//		}
//		if (!ismeet)
//			continue;
//		else
//			cr_num++;
//	}
//	//对集群节点排序
//	sort(realpos, realpos + nodeNum, cmp2);
//	mysort(realpos, nodeNum);
//	//Point2d errpos[nodeNum];
//	//for (int i = 0; i < nodeNum; ++i) {
//	//	errpos[i] = realpos[i];
//	//	posAddErr(errpos[i], 1500);//叠加误差
//	//}
//	for (int i = 0; i < nodeNum; ++i) {
//		for (int j = i + 1; j < nodeNum; ++j) {
//			d.push_back(Dist(realpos[i], realpos[j]));//根据真值计算测距距离
//		}
//	}
//	Point2d tgtpos[tgtNum];
//	double Dis[nodeNum][tgtNum] = { 0 };
//	double Dir[nodeNum][tgtNum] = { 0 };
//	int genTgtNum = 0;
//	//随机生成目标坐标真值
//	while (1) {
//		if (genTgtNum > tgtNum - 1)
//			break;
//		tgtpos[genTgtNum] = Point2d(200000.0 * rand() / RAND_MAX - 100000, 50000.0 * rand() / RAND_MAX - 100000);
//		bool flag = true;
//		for (int i = 0; i < nodeNum; i++) {
//			Dis[i][genTgtNum] = Dist(tgtpos[genTgtNum], realpos[i]);
//			if (Dis[i][genTgtNum] < 80000) {
//				flag = false;
//				break;
//			}
//
//		}
//		if (!flag)
//			continue;
//		else
//			genTgtNum++;
//	}
//	//根据真值计算目标方向与距离并叠加误差（可看做相对信息）
//	tgtsDirDis(realpos, tgtpos, Dir, Dis, nodeNum, tgtNum);
//	double bta[nodeNum - 2] = { 0 };
//	double cosA[nodeNum - 2] = { 0 };
//	Point2d RESpoint[nodeNum];
//	double theta = 0;
//	int SCALE = 250;//缩放尺度
//	for (int jx = 0; jx < 2; jx++) {//组网存在镜像
//		for (int st = 0; st < 360; st++) {
//			theta = st * CV_PI / 180;
//			if (jx == 0)
//				for (int i = 0; i < nodeNum - 2; ++i) {
//					cosA[i] = (d[0] * d[0] + d[i + 1] * d[i + 1] - d[nodeNum - 1 + i] * d[nodeNum - 1 + i]) / (2 * d[0] * d[i + 1]);
//					bta[i] = acos(cosA[i]);
//				}
//			else
//				for (int i = 0; i < nodeNum - 2; ++i) {
//					cosA[i] = (d[0] * d[0] + d[i + 1] * d[i + 1] - d[nodeNum - 1 + i] * d[nodeNum - 1 + i]) / (2 * d[0] * d[i + 1]);
//					bta[i] = -acos(cosA[i]);
//				}
//			Mat display = Mat::zeros(1000, 1000, CV_8UC3);
//			//节点1在相对坐标系的位置
//			RESpoint[0].x = 0;
//			RESpoint[0].y = 0;
//			//根据组网测距信息与组网姿态角度计算其他节点相对坐标
//			RESpoint[1] = Point2d(RESpoint[0].x + d[0] * cos(theta), RESpoint[0].y + d[0] * sin(theta));
//			for (int i = 2; i < nodeNum; ++i) {
//				RESpoint[i] = Point2d(RESpoint[0].x + d[i - 1] * cos(bta[i - 2] + theta), RESpoint[0].y + d[i - 1] * sin(bta[i - 2] + theta));
//			}
//			//绘图
//			for (int i = 0; i < nodeNum; i++) {
//				if ((i + 1) == nodeNum)
//					line(display, Point(RESpoint[i].x / SCALE + display.rows / 2, RESpoint[i].y / SCALE + display.cols / 2), Point(RESpoint[0].x / SCALE + display.rows / 2, RESpoint[0].y / SCALE + display.cols / 2),
//						Scalar(0, 255, 0), 2);
//				else
//					line(display, Point(RESpoint[i].x / SCALE + display.rows / 2, RESpoint[i].y / SCALE + display.cols / 2), Point(RESpoint[i + 1].x / SCALE + display.rows / 2, RESpoint[i + 1].y / SCALE + display.cols / 2),
//						Scalar(0, 255, 0), 2);
//				putText(display, to_string(i + 1), Point(RESpoint[i].x / SCALE + display.rows / 2, RESpoint[i].y / SCALE + display.cols / 2), FONT_HERSHEY_COMPLEX, 2, cv::Scalar(0, 255, 255), 3, 8, 0);
//				Point2d tempTgtPos;
//				for (int j = 0; j < tgtNum; j++) {
//					tempTgtPos.x = RESpoint[i].x + Dis[i][j] * cos(Dir[i][j]);
//					tempTgtPos.y = RESpoint[i].y + Dis[i][j] * sin(Dir[i][j]);
//					line(display, Point(RESpoint[i].x / SCALE + display.rows / 2, RESpoint[i].y / SCALE + display.cols / 2), Point(tempTgtPos.x / SCALE + display.rows / 2, tempTgtPos.y / SCALE + display.cols / 2),
//						Scalar(255, 255, 0), 2);
//					circle(display, Point(tempTgtPos.x / SCALE + display.rows / 2, tempTgtPos.y / SCALE + display.cols / 2), 3, Scalar(255, 0, 0), 3);
//
//				}
//			}
//			imshow("rol", display);
//			waitKey(100);
//		}
//	}
//
//}
//
