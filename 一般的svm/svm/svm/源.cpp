/*********************************************************
**CopyRight by zhanghuanhuan, BNU in china, BeiJing ******
**********************************************************/

#include<iostream>
#include<vector>
#include<math.h>
using namespace std;
double zer0 = 0.000001;
double zerob = 0.000000001;
struct arry
{
	double x1;
	double x2;
	int  y;
};
struct W{
	double w1;
	double w2;
};

//求内积
double Inner_Product(arry a, arry b)
{
	return a.x1*b.x1 + a.x2*b.x2;
}

//求Max
double Svm_max(double a[][2], vector<arry> x)
{
	int i, j;
	double sumai = 0;
	double sum = 0;
	for (i = 0; i < x.size(); i++)
	{
		sumai = sumai + a[i][1];
	}
	for (i = 0; i < x.size(); i++)
	{
		for (j = 0; j < x.size(); j++)
			sum = sum + a[i][1] * a[j][1] * x[i].y*x[j].y*Inner_Product(x[i],x[j]);
	}
	return sumai - 0.5*sum;

}

//求C
double Svm_C(vector<arry> x, double a[][2], int m)
{
	double c = 0;
	int i;
	for (i = 0; i < m; i++)
	{
		if (a[i][0] == 0)
			c = c + a[i][1] * x[i].y;
	}
	return -c;
}

//计算系数项
double Svm_k_bi_b0(vector<arry> x, double a[][2], int m, arry xi)
{
	double k = 0;
	int i;
	for (i = 0; i < m; i++)
	{
		if (a[i][0] == 0)
			k = k + 0.5*xi.y*a[i][1] * x[i].y*Inner_Product(xi,x[i]);
	}
	return 1-k;
}

double Svm_k_bj_b1(vector<arry> x, double a[][2], int m, arry xj)
{
	double k = 0;
	int i;
	for (i = 0; i < m; i++)
	{
		if (a[i][0] == 0)
			k = k + 0.5*xj.y*a[i][1] * x[i].y*Inner_Product(xj, x[i]);
	}
	return 1 - k;
}

//求yiyjxixj,g
double Svm_g(arry xi, arry xj)
{
	return xi.y*xj.y*Inner_Product(xi, xj);
}


//求n,f
double Svm_ni_fj(arry a)
{
	return Inner_Product(a, a);
 }


//一次项系数

double Svm_B(double b0, double b1, double c, double n, double g, arry xi, arry xj)
{
	double b;
	b = -b0*xj.y / xi.y + b1 + n*xi.y*c / (xi.y*xi.y) - g*c / xi.y;
	return b;
}

//二次项系数
double Svm_A(double g, double n, double f, arry xi, arry xj)
{
	return g*xj.y / xi.y - n*xj.y*xj.y / (2 * xi.y*xi.y) - f / 2;
}

//跟新的aj
double Svm_Update_Aj(vector<arry> x, double a[][2], arry xi, arry xj, int m, double c)
{
	double b0, b1, n, g, f;
	double A, B;
	//c = Svm_C(x, a, m);
	b0 = Svm_k_bi_b0(x, a, m, xi);
	b1 = Svm_k_bj_b1(x, a, m, xj);
	n = Svm_ni_fj(xi);
	f = Svm_ni_fj(xj);
	g = Svm_g(xi,xj);
	A = Svm_A(g,  n,  f, xi,  xj);
	B = Svm_B(b0, b1, c, n, g, xi, xj);
	return -B / (2 * A);
}

//跟新的ai
double Svm_Update_Ai(double aj, double c, arry xi, arry xj)
{
	return (c - aj*xj.y) / xi.y;
}


//跟新的w
W Svm_Update_w(vector<arry> x, double a[][2])
{
	W w;
	w.w1 = w.w2 = 0;
	int i;
	for (i = 0; i < x.size(); i++)
	{
		w.w1 = w.w1 + a[i][1] * x[i].y*x[i].x1;
		w.w2 = w.w2 + a[i][1] * x[i].y*x[i].x2;
	}
	return w;
}

double Svm_Update_b(vector<arry> x, double a[][2], W w, double b)
{
	double s = 0;
	int i;
	int n=0;
	for (i = 0; i < x.size(); i++)
	{
		if (a[i][1] >= 0 && (x[i].y*(w.w1*x[i].x1 + w.w2*x[i].x2 + b) - 1) >= zerob
			&& abs(a[i][1] * (x[i].y*(w.w1*x[i].x1 + w.w2*x[i].x2 + b) - 1))< zerob)
		{
			n++;
			s = s + x[i].y - w.w1*x[i].x1 + w.w2*x[i].x2;
		}
	}
	if (n == 0) return 0;
	else return s / n;
}
double Svm_K12_k11_k22(arry x1, arry  x2)
{
	return Inner_Product(x1, x1) + Inner_Product(x2, x2) - 2 * Inner_Product(x1,x2);
}


//ej-ei
double  Svm_ej_ei(arry xj, arry xi, W w)
{
	return xi.y - xj.y + w.w1*(xj.x1 - xi.x1) + w.w2*(xj.x2 - xi.x2);
}

double Svm_upa_aj(arry xi, arry xj, W w, double aj)
{

	return aj + xj.y*Svm_ej_ei(xi, xj,w)/Svm_K12_k11_k22(xi,xj);
}
//ej-ei is max ei，第二个ａi的固定,a[][10]是为了保证前后的两个是否保证目标函数发生改变。
int Svm_Ej_Mai(vector<arry> x, int j, W w, double aij[][10])
{
	double max;
	int book=-1;
	int i;
	max = 0;
	for (i = 0; i < x.size(); i++)
	{
		if (abs(Svm_ej_ei(x[j], x[i], w)) >max &&aij[j][i]!=1)
		{
			max = abs(Svm_ej_ei(x[j], x[i], w));
			book = i;
		}
		//if (a[i][0] == 1) a[i][0] = 0;
	}
	return book;
}

//第一个不满足kkt的训练样本ｊ和极端情况下的aj通过oldkj来判断
int Svm_First_aj(vector<arry> x, double a[][2], W w, double b,double *oldkj)
{
	int book=-2;
	int i;
	double  ab;
	double temp;
	double max = 0;
	double gx = 0;
	for (i = 0; i < x.size(); i++)
	{
		
		gx = x[i].y*(w.w1*x[i].x1 + w.w2*x[i].x2 + b) - 1;
		if (a[i][1] < zer0 && (abs(gx)<zer0||gx>0)
			|| a[i][1] > zer0&&abs(gx) < zer0)
			continue;
		else 
		{
			if (oldkj[i] != 1)
			{
				book = i;
				break;
			}
			continue;
		}
		
	}
	if (book<10&&book!=-2)
	return book;
	else return -2;
}
double Svm_ai(vector<arry> x, double a[][2],double ajold, int  ki, int kj)
{
	return a[ki][1] + x[ki].y*x[kj].y*(ajold - a[kj][1]);
}

double Svm_binew(vector<arry> x, double a[][2],int ai)
{
	int i = 0;
	double sum = 0;
	for (i = 0; i < x.size(); i++)
		sum = sum + a[i][1] * x[i].y*Inner_Product(x[i], x[ai]);
	return x[ai].y - sum;

}
int main()
{
	vector<arry> x;
	int count = 0;
	x.reserve(sizeof(arry)*10);
	int m = 10;
	int inter = 0;
	int i;
	double bb = 0;
	double a[10][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 },
	{ 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 } };
	arry b[10]; 
	b[0] = { 0, 0, 1 };
	b[1] = { -2, -2, 1 };
	b[2] = { 0, 4, -1 };
	b[3] = { 4, 0, -1};
	b[4] = { 2.5, 2, -1 };
	b[5] = { -3, -3, 1 };
	b[6] = { 3, 3, -1 };
	b[7] = { -5, -6, 1 };
	b[8] = { 100, 80, -1 };
	b[9] = { -100, -100, 1 };
	W w;
	double oldbb = 0;
	for (i = 0; i < 10; i++)
	{   
		x.push_back(b[i]);
	}
	w = Svm_Update_w(x, a);
	int ki = 0;
	int kj = 1;
	int tempki=0;
	int tempkj=0;
	W oldw;
	oldw = {2,2};
	double aj, ai;
	double c = 0;
	double h = 0;
	double ajold = 0;
	double akid = 0;
	double akjd = 0;
	double akjaki[10][10] = {0};
	double oldkj[10] = {0};
	while (inter <120)
	{  
		//考虑极端情况，aj中下一个ai都不满足条件，更新aj,再次寻找aj,一直到找到为止。
		//找到合适的aj后，然后跟新把aj的指示数全部归0
		do
		{
			kj = Svm_First_aj(x, a, w, bb, oldkj);
			if (kj == -2) break;
			ki = Svm_Ej_Mai(x, kj, w, akjaki);
			if (ki == -1)  oldkj[kj] = 1;
			if (ki != -1)
			{
				for (i = 0; i < x.size(); i++)
					oldkj[kj] = 0;
			}
		} while (ki == -1);
		if (kj == -2) break;
			a[ki][0] = 1;
			a[kj][0] = 1;
			c = c = Svm_C(x,a,10);
		    a[ki][0] = akid;
		    a[kj][0] = akjd;
			ajold = a[kj][1];
			//h = Svm_Update_Aj(x,a,x[ki], x[kj],10, c);
			cout << Svm_Update_Aj(x, a, x[ki], x[kj], 10, c) << endl;
		  h = Svm_upa_aj(x[ki], x[kj], w, a[kj][1]);
			if (x[kj].y != x[ki].y)
			{
				double L;
				L = (a[kj][1] - a[ki][1]) ? (a[kj][1] - a[ki][1]) : 0;
				if (h >= L)
					aj = h;
				else aj = L;
			}
			if (x[kj].y == x[ki].y)
			{
				if (0>=h)
					a[kj][1] = 0;
				if (h> 0&&h <(a[ki][1] + a[kj][1]))
					aj = h;
				if (h>(a[ki][1] + a[kj][1]))
					aj = (a[ki][1] + a[kj][1]);

			}
			a[kj][1] = aj;
			//ai = Svm_ai(x,a,ajold,ki,kj);
			ai = Svm_Update_Ai(aj, c, x[ki], x[kj]);
			a[ki][1] = ai;
			//a[kj][1] = aj;	
			oldw =w;
			oldbb = bb;
			w = Svm_Update_w(x, a);	
			if (aj >= 0)
				bb = Svm_binew(x, a, kj);
			else if(ai>=0)
				bb = Svm_binew(x, a, ki);
			else bb = (Svm_binew(x, a, kj) + Svm_binew(x, a, ki))*0.5;
			
			//if ((abs(w.w1 - oldw.w1) + abs(w.w2 - oldw.w2)) < zerob)
			if (((w.w1*w.w1 - oldw.w1*oldw.w1) + (w.w2*w.w2 - oldw.w2*oldw.w2)) < zerob)
			{
				//a[ki][0] = 1;
				akjaki[kj][ki] = 1;
			}
			else
			{
				for (i = 0; i < x.size();i++)
					akjaki[kj][i] = 0;
			}
			cout << Svm_max(a, x) << endl;
			cout << w.w1 << " " << w.w2 << endl;
			cout << bb << endl;
			inter++;
	}
	cout << Svm_max(a, x) << endl;
	cout << w.w1 << " " << w.w2 << endl;
	cout << bb << endl;
	cout << inter << endl;
	for (i = 0; i < 10; i++)
		cout << a[i][1]<< endl;
	system("pause");
	return 0;
}

