#include <iostream>
#include <vector>
#include "Degree.h"
using namespace std;


int upperBounds(Degree &degree, vector<int> h, int n, int d, int r) {
	int res = degree.upperBoundUnivariate(r, n, d, h);
	cout << "results ----- degree at round" << r << ": " << res << " - solving time:" << degree.max_time << endl;
	degree.min_time = 100000;
	degree.max_time = 0;
	return res;
}

void findOptimalAffineLayers(Degree& degree, int r, int n, int d, int w) {
	degree.searchParameters(r, n, d, w);
	degree.outputTime();
	degree.min_time = 100000;
	degree.max_time = 0;
}

void result1InPaper(Degree &degree) {
	cout << "results for Chaghri:" << endl;
	findOptimalAffineLayers(degree, 6, 63, 32, 3);
}

void result2InPaper(Degree& degree) {
	cout << "results for MiMC:" << endl;
	findOptimalAffineLayers(degree, 7, 129, 1, 4);
}

void result3InPaper(Degree& degree) {
	cout << "tests upper bounds for different (h1,h2,h3) where (n,d)=(63,32):" << endl;
	int arr[5][2] = { {0,3},{0,6},{0,36} };
	vector<int> h(2);
	for (int i = 0; i < 3; i++) {
		cout << "h1, h2:";
		for (int j = 0; j < 2; j++) {
			h[j] = arr[i][j];
		}
		cout << h[0]<<" "<<h[1]<< endl;
		int n = 63, d = 32;
		for (int r = 1; r < 30; r++) {
			int res = degree.upperBoundUnivariate(r, n, d, h);
			cout << "round " << r << ":" << res << " - solving time:" << degree.max_time << endl;
			degree.min_time = 100000;
			degree.max_time = 0;
			if (res == n)
				break;
		}
		cout << endl;
	}
}

void result4InPaper(Degree& degree) {
	cout << "tests upper bounds for different (h1,h2,h3,h4) where (n,d)=(129,1):" << endl;
	int arr[5][2] = { {0,6},{0,9}, {0,36},{0,40}, {0,63} };
	vector<int> h(2);
	for (int i = 0; i < 5; i++) {
		cout << "h1, h2:";
		for (int j = 0; j < 2; j++) {
			h[j] = arr[i][j];
		}
		cout << h[0] << " " << h[1] << endl;
		int n = 129, d = 1;
		for (int r = 1; r < 30; r++) {
			int res = degree.upperBoundUnivariate(r, n, d, h);
			cout << "round " << r << ":" << res << " - solving time:" << degree.max_time << endl;
			degree.min_time = 100000;
			degree.max_time = 0;
			if (res == n)
				break;
		}
		cout << endl;
	}
}

int main() {
	Degree degree;

	int cmd = 1;
	cout << "please input your command:" << endl;
	cout << "0 --------> exit" << endl;
	cout << "1 --------> find optimal (h1, h2, h3) for (n, d, r)=(63,32,6)" << endl;
	cout << "2 --------> find optimal (h1, h2, h3, h4) for (n, d, r)=(129,1,7)" << endl;
	cout << "3 --------> test upper bounds for different (h1,h2) where (n,d)=(63,32)" << endl;
	cout << "4 --------> test upper bounds for different (h1,h2) where (n,d)=(129,1)" << endl;
	cout << "5 --------> find optimal (h1, ..., h_w) for your choice of (n,d,r,w)" << endl;
	cout << "6 --------> test upper bounds for your choice of (h1, ..., h_w) and (n,d,r)" << endl;
	while (cin >> cmd) {
		if (cmd == 1)
			result1InPaper(degree);
		else if (cmd == 2)
			result2InPaper(degree);
		else if (cmd == 3)
			result3InPaper(degree);
		else if (cmd == 4)
			result4InPaper(degree);
		else if (cmd == 5) {
			int n = 0, d = 0, r = 0, w = 0;
			cout << "input n,d,r,w (seperated with space):";
			cin >> n >> d >> r >> w;
			findOptimalAffineLayers(degree, r, n, d, w);
		}
		else if (cmd == 6) {
			degree.env.set(GRB_IntParam_OutputFlag, 1);//you set set it to 0 if you do not want to see the details of the solving procedure
			int n = 0, d = 0, r = 0, w = 0;
			cout << "input n,d,r,w (seperated with space):";
			cin >> n >> d >> r >> w;
			vector<int> h(w);
			cout << "please input h1, ..., h" << w << " (seperated with space):";
			for (int i = 0; i < w; i++) {
				cin >> h[i];
			}
			upperBounds(degree, h, n, d, r);
			cout << endl;
			degree.env.set(GRB_IntParam_OutputFlag, 0);
		}
		else {
			break;
		}

		cout << "please input your command:" << endl;
		cout << "0 --------> exit" << endl;
		cout << "1 --------> find optimal (h1, h2, h3) for (n, d, r)=(63,32,6)" << endl;
		cout << "2 --------> find optimal (h1, h2, h3, h4) for (n, d, r)=(129,1,7)" << endl;
		cout << "3 --------> test upper bounds for different (h1,h2) where (n,d)=(63,32)" << endl;
		cout << "4 --------> test upper bounds for different (h1,h2) where (n,d)=(129,1)" << endl;
		cout << "5 --------> find optimal (h1, ..., h_w) for your choice of (n,d,r,w)" << endl;
		cout << "6 --------> test upper bounds for your choice of (h1, ..., h_w) and (n,d,r)" << endl;
	}
}