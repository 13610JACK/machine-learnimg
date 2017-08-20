#include<iostream>
#include<vector>
#include<math.h>
using namespace std;
float arfa = 0.001;
struct arry
{
	float x1;
	float x2;
};
struct W
{
	float w1;
	float w2;
};
float error(vector<arry> b, W w)
{
	float errorsum = 0;
	int i;
	for (i = 0; i < b.size(); i++)
		errorsum = errorsum + (b[i].x2 - w.w1*b[i].x1 - w.w2)*(b[i].x2 - w.w1*b[i].x1 - w.w2);
	return errorsum;
}
W updatew(vector<arry> a, W w)
{
	int i;
	W ww;
	float upwsumw1 = 0;
	float upwsumw2 = 0;
	for (i = 0; i < a.size(); i++)
	{
		upwsumw1 = upwsumw1 + a[i].x1*(a[i].x2 - w.w1*a[i].x1 - w.w2);
		upwsumw2 = upwsumw2 + (a[i].x2 - w.w1*a[i].x1 - w.w2);
	}
	ww.w1 = w.w1 + arfa*upwsumw1;
	ww.w2 = w.w2 + arfa*upwsumw2;
	return ww;
}
int main()
{
	vector<arry> a;
	int count = 0;
	a.reserve(sizeof(arry) * 7);
	cout << a.capacity() << endl;
	float error0 = 0, error1 = 0;
	W w;
	w.w1 = 0;
	w.w2 = 0;
	int i;
	arry b;
	//vector<float> b;
	/*b[0] = { 1, 3 };
	b[1] = { 2, 5 };
	b[2] = { 4, 9 };
	b[3] = { 0, 1 };
	b[4] = { 8, 17 };
	b[5] = { 9, 19 };
	b[6] = { 1.5, 4 };*/
	for (i = 0; i < 7; i++)
	{
		cin >> b.x1 >> b.x2;
		a.push_back(b);
	}
	error0 = error(a, w);
	do
	{
		w = updatew(a, w);
		count++;
		//error1 = error0;
		error0 = error(a, w);
	} while (error0>0.00000001);//abs(error1 - error0) > 0.000001);
	cout << w.w1 << " " << w.w2 << "   " << count;
	system("pause");
	return 0;
}