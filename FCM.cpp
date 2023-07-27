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
//	int di;                      // ���ݵ��ά�ȣ����磺λ��������ٶ������ܹ���4��ά�ȣ�
//	vector<double> datas;        // �������� �洢������λ�á��ٶ�
//
//	Mode() : di(4) {} // ���캯����ʼ��
//};
//
//typedef vector<vector<Mode>> ModeVec;
//
//const int N = 1000;
//const double eps = 1e-2;
//const double eps_dis = 1e-6;
//
//double getDistance(Mode &m1, Mode &m2)//����ŷʽ����û�п���
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
//	int di = p[0].di;        // ��ȡ���ݵ��ά��
//
//							 // ��ʼ����������
//	for (int i = 0; i < clusternum; i++)
//	{
//		c.emplace_back();
//		c[i].di = di;
//		c[i].datas.resize(di, 0);
//	}
//
//	// ����������ĵĳ�ʼֵ
//	for (int i = 0; i < clusternum; i++)
//	{
//		// �����ʵ����ݵ�����Ը�����Ӧ�ľ�����������������
//		for (int j = 0; j < di; j++)
//		{
//			// ������ƫ�ƣ���ΧΪ[-0.01, 0.01]
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
//	vector<Mode> c(clusternum); //��������
//	vector<vector<double>> dis(n, vector<double>(clusternum)); //����
//	vector<vector<double>> u(n, vector<double>(clusternum));   //�����Ⱦ���
//	init_c(p, n, clusternum, c);
//	comp_dis(p, c, n, clusternum, dis);
//	while (true)
//	{
//		index++;
//		printf("��%d�ε���----------------------------------------\n�����Ⱦ���\n", index);
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
//		printf("J����ֵ=%f\n", sum);
//		if (fabs(psum - sum) < eps)
//			break;
//	}
//	ans.clear();
//	for (int i = 0; i <= clusternum; i++)//����c+1�� ��c+1�����ڴ洢�����趨������ֵ�ĵ�
//	{
//		vector<Mode> m;
//		ans.push_back(m);
//	}
//	for (int i = 0; i < n; i++)
//	{
//		double tp = -1;
//		int index = -1; // ��ʼ��Ϊ-1����ʾ���ݵ�Ϊ�µĶ�������
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
//			// ����һ����ֵ�����ݵ㣬�������о���
//			ans[index].push_back(p[i]);
//		}
//		else
//		{
//			// ����һ����ֵ�����ݵ㣬��Ϊ�µĶ�������
//			ans[clusternum].push_back(p[i]);
//		}
//	}
//
//}
//
//int main()
//{
//	int n = 30;          // ���ú�����ĸ���
//	int dimension = 2;   // ���ú���ά��
//	int clusternum = 6;  // ���þ�����Ŀ
//	const double membershipThreshold = 0.8; // ������������ֵ
//
//	vector<Mode> p;
//	srand(time(NULL));
//
//	for (int i = 0; i < n; i++)
//	{
//		p.emplace_back();
//		p[i].di = dimension;
//		if (i >= clusternum)//ʹ���������Ϊclusternum��һ��
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
//				double temp = (rand() % 500); // ������ɺ���
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
//		printf("Ŀ��%d��\n", i + 1);
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
//	// ����һ��ͼ�������ƺ�����
//	int imgWidth = 500;
//	int imgHeight = 500;
//	Mat image(imgHeight, imgWidth, CV_8UC3, Scalar(0, 0, 0));
//
//	// ���㺽�������������
//	double maxX = 0, maxY = 0;
//	for (int i = 0; i < n; i++)
//	{
//		maxX = max(maxX, p[i].datas[0]);
//		maxY = max(maxY, p[i].datas[1]);
//	}
//	double scaleX = imgWidth / (maxX + 1);
//	double scaleY = imgHeight / (maxY + 1);
//
//	//��ͼ���ϻ��ƺ�����
//	//��ɫ
//	Scalar colors[] = { Scalar(0, 0, 255),    // ��ɫ
//		Scalar(0, 255, 0),    // ��ɫ
//		Scalar(0, 128, 255),  // ��ɫ
//		Scalar(0, 255, 255),  // ��ɫ
//		Scalar(255, 255, 0),  // ��ɫ
//		Scalar(255, 0, 255),  // ���ɫ
//		Scalar(0, 0, 128),    // ���ɫ
//		Scalar(0, 128, 0),    // ����ɫ
//		Scalar(255, 0, 128),  // ��ɫ
//		Scalar(128, 128, 0),  // ���ɫ
//		Scalar(128, 0, 128),  // ����ɫ
//		Scalar(128, 128, 0),  // ����ɫ
//		Scalar(128, 255, 0),  // ����ɫ
//		Scalar(255, 0, 0),    // ����ɫ
//		Scalar(128, 128, 255), // ǳ��ɫ
//		Scalar(255, 0, 0),    // ��ɫ
//		Scalar(128, 255, 128), // ǳ��ɫ
//		Scalar(255, 128, 128), // ǳ��ɫ
//		Scalar(128, 255, 255), // ǳ��ɫ
//		Scalar(255, 128, 255)  // ǳ���ɫ
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
//				// ���ڶ�����һ�࣬�ú�ɫ�㻭��,������С����ֵ�ᱻ�ж�Ϊ����
//				circle(image, pt, 5, Scalar(255, 255, 255), -1);
//			}
//			else
//			{
//				// ���ڷǶ����ľ��࣬�ö�Ӧ����ɫ����
//				circle(image, pt, 5, colors[i], -1);
//			}
//
//		}
//	}
//
//	imshow("��ǰʱ�̺��������", image);
//	waitKey(0);
//
//	return 0;
//}
