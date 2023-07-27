//#define  _CRT_SECURE_NO_WARNINGS
//#include <iostream>
//#include <vector>
//#include <cmath>
//#include <ctime>
//#include <cstdlib>
//#include <opencv2/opencv.hpp>
//
//using namespace std;
//using namespace cv;
//
//struct Mode
//{
//	int di;                      // 数据点的维度（例如：位置坐标和速度坐标总共有4个维度）
//	vector<double> datas;        // 属性数据 存储向量如位置、速度
//
//	Mode() : di(4) {} // 构造函数初始化
//};
//
//typedef vector<vector<Mode>> ModeVec;
//
//const int N = 1000;
//const double eps = 1e-2;
//const double eps_dis = 1e-6;
//
//double getDistance(Mode &m1, Mode &m2)//计算欧式距离没有开方
//{
//	int di = m1.di;
//	double ans = 0;
//	for (int i = 0; i < di; i++)
//		ans += (m1.datas[i] - m2.datas[i]) * (m1.datas[i] - m2.datas[i]);
//	return ans;
//}
//
//void init_c(vector<Mode> &p, int n, int clusternum, vector<Mode> &c)
//{
//	int di = p[0].di;        // 获取数据点的维度
//
//							 // 初始化聚类中心
//	for (int i = 0; i < clusternum; i++)
//	{
//		c.emplace_back();
//		c[i].di = di;
//		c[i].datas.resize(di, 0);
//	}
//
//	// 计算聚类中心的初始值
//	for (int i = 0; i < clusternum; i++)
//	{
//		// 将合适的数据点的属性赋到对应的聚类中心数据向量上
//		for (int j = 0; j < di; j++)
//		{
//			// 添加随机偏移，范围为[-0.01, 0.01]
//			double randomOffset = (double(rand() % 200) - 100) / 10000.0;
//			c[i].datas[j] = p[i].datas[j] * (1.0 + randomOffset);
//		}
//	}
//}
//
//void comp_dis(vector<Mode> &p, vector<Mode> &c, int n, int clusternum, vector<vector<double>> &dis)
//{
//	for (int i = 0; i < n; i++)
//		for (int j = 0; j < clusternum; j++)
//			dis[i][j] = getDistance(p[i], c[j]);
//}
//
//void comp_u(vector<vector<double>> &dis, int n, int clusternum, vector<vector<double>> &u)
//{
//	for (int i = 0; i < n; i++)
//	{
//		double tp = 0;
//		for (int j = 0; j < clusternum; j++)
//		{
//			if (dis[i][j] < eps_dis)
//			{
//				for (int k = 0; k < clusternum; k++)
//					u[i][k] = 0;
//				u[i][j] = 1;
//				return;
//			}
//			tp += 1 / dis[i][j];
//		}
//		tp = 1 / tp;
//		for (int j = 0; j < clusternum; j++)
//			u[i][j] = tp * (1 / dis[i][j]);
//	}
//}
//
//void update_c(vector<Mode> &p, vector<vector<double>> &u, int n, int clusternum, vector<Mode> &c)
//{
//	int di = p[0].di;
//	for (int j = 0; j < clusternum; j++)
//	{
//		c[j].di = di;
//		c[j].datas.clear();
//		c[j].datas.resize(di, 0);
//		double tp = 0;
//		for (int i = 0; i < n; i++)
//		{
//			for (int k = 0; k < di; k++)
//				c[j].datas[k] += u[i][j] * u[i][j] * p[i].datas[k];
//			tp += u[i][j] * u[i][j];
//		}
//		for (int k = 0; k < di; k++)
//			c[j].datas[k] /= tp;
//	}
//}
//
//double comp_obj_func(vector<vector<double>> &u, vector<vector<double>> &dis, int n, int clusternum, int di)
//{
//	double sum = 0;
//	for (int i = 0; i < n; i++)
//		for (int j = 0; j < clusternum; j++)
//			sum += u[i][j] * u[i][j] * dis[i][j];
//	return sum;
//}
//
//void FCM(vector<Mode> &p, int n, int di, int clusternum, vector<vector<Mode>> &ans, double membershipThreshold)
//{
//	int index = 0;
//	double sum = 0, psum;
//	vector<Mode> c(clusternum); //聚类中心
//	vector<vector<double>> dis(n, vector<double>(clusternum)); //距离
//	vector<vector<double>> u(n, vector<double>(clusternum));   //隶属度矩阵
//	init_c(p, n, clusternum, c);
//	comp_dis(p, c, n, clusternum, dis);
//	while (true)
//	{
//		index++;
//		printf("第%d次迭代----------------------------------------\n隶属度矩阵：\n", index);
//		comp_u(dis, n, clusternum, u);
//		for (int i = 0; i < n; i++)
//		{
//			for (int j = 0; j < clusternum; j++)
//			{
//				printf("%f ", u[i][j]);
//			}
//			printf("\n");
//		}
//		update_c(p, u, n, clusternum, c);
//		comp_dis(p, c, n, clusternum, dis);
//		psum = sum;
//		sum = comp_obj_func(u, dis, n, clusternum, di);
//		printf("J函数值=%f\n", sum);
//		if (fabs(psum - sum) < eps)
//			break;
//	}
//	ans.clear();
//	for (int i = 0; i <= clusternum; i++)//生成c+1类 第c+1类用于存储低于设定关联阈值的点
//	{
//		vector<Mode> m;
//		ans.push_back(m);
//	}
//	for (int i = 0; i < n; i++)
//	{
//		double tp = -1;
//		int index = -1; // 初始化为-1，表示数据点为新的独立聚类
//		for (int j = 0; j < clusternum; j++)
//		{
//			if (u[i][j] > tp && u[i][j] >= membershipThreshold)
//			{
//				tp = u[i][j];
//				index = j;
//			}
//		}
//
//		if (index != -1)
//		{
//			// 高于一定阈值的数据点，归入现有聚类
//			ans[index].push_back(p[i]);
//		}
//		else
//		{
//			// 低于一定阈值的数据点，作为新的独立聚类
//			ans[clusternum].push_back(p[i]);
//		}
//	}
//
//}
//
//int main()
//{
//	int n = 30;          // 设置航迹点的个数
//	int dimension = 2;   // 设置航迹维度
//	int clusternum = 6;  // 设置聚类数目
//	const double membershipThreshold = 0.8; // 定义隶属度阈值
//
//	vector<Mode> p;
//	srand(time(NULL));
//
//	for (int i = 0; i < n; i++)
//	{
//		p.emplace_back();
//		p[i].di = dimension;
//		if (i >= clusternum)//使仿真的数据为clusternum个一组
//		{
//			for (int j = 0; j < dimension; j++)
//			{
//				double temp = p[i - clusternum].datas[j] + (double(rand() % 200) - 100) / 3.0;
//				p[i].datas.push_back(temp);
//			}
//		}
//		else
//		{
//			for (int j = 0; j < dimension; j++)
//			{
//				double temp = (rand() % 500); // 随机生成航迹
//				p[i].datas.push_back(temp);
//			}
//		}
//
//	}
//
//	vector<vector<Mode>> ans;
//	vector<Mode> outliers;
//	FCM(p, n, dimension, clusternum, ans, membershipThreshold);
//
//	for (int i = 0; i < clusternum; i++)
//	{
//		printf("目标%d：\n", i + 1);
//		for (int j = 0; j < ans[i].size(); j++)
//		{
//			printf("(");
//			for (int k = 0; k < dimension; k++)
//			{
//				if (k == 0)
//					printf("%.2f", ans[i][j].datas[k]);
//				else
//					printf(",%.2f", ans[i][j].datas[k]);
//			}
//			printf(") ");
//		}
//		printf("\n");
//	}
//
//	// 创建一个图像来绘制航迹点
//	int imgWidth = 500;
//	int imgHeight = 500;
//	Mat image(imgHeight, imgWidth, CV_8UC3, Scalar(0, 0, 0));
//
//	// 计算航迹点的缩放因子
//	double maxX = 0, maxY = 0;
//	for (int i = 0; i < n; i++)
//	{
//		maxX = max(maxX, p[i].datas[0]);
//		maxY = max(maxY, p[i].datas[1]);
//	}
//	double scaleX = imgWidth / (maxX + 1);
//	double scaleY = imgHeight / (maxY + 1);
//
//	//在图像上绘制航迹点
//	//颜色
//	Scalar colors[] = { Scalar(0, 0, 255),    // 红色
//		Scalar(0, 255, 0),    // 绿色
//		Scalar(0, 128, 255),  // 橙色
//		Scalar(0, 255, 255),  // 黄色
//		Scalar(255, 255, 0),  // 青色
//		Scalar(255, 0, 255),  // 洋红色
//		Scalar(0, 0, 128),    // 深红色
//		Scalar(0, 128, 0),    // 深绿色
//		Scalar(255, 0, 128),  // 紫色
//		Scalar(128, 128, 0),  // 深黄色
//		Scalar(128, 0, 128),  // 深紫色
//		Scalar(128, 128, 0),  // 深青色
//		Scalar(128, 255, 0),  // 青绿色
//		Scalar(255, 0, 0),    // 深蓝色
//		Scalar(128, 128, 255), // 浅红色
//		Scalar(255, 0, 0),    // 蓝色
//		Scalar(128, 255, 128), // 浅绿色
//		Scalar(255, 128, 128), // 浅蓝色
//		Scalar(128, 255, 255), // 浅黄色
//		Scalar(255, 128, 255)  // 浅洋红色
//
//	};
//
//	for (int i = 0; i <= clusternum; i++)
//	{
//		for (int j = 0; j < ans[i].size(); j++)
//		{
//			Point pt(ans[i][j].datas[0] * scaleX, ans[i][j].datas[1] * scaleY);
//
//			if (i == clusternum)
//			{
//				// 对于独立的一类，用黑色点画出,隶属度小于阈值会被判定为独立
//				circle(image, pt, 5, Scalar(255, 255, 255), -1);
//			}
//			else
//			{
//				// 对于非独立的聚类，用对应的颜色画出
//				circle(image, pt, 5, colors[i], -1);
//			}
//
//		}
//	}
//
//	imshow("当前时刻航迹点关联", image);
//	waitKey(0);
//
//	return 0;
//}
