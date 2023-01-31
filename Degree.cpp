#include "Degree.h"
#include <vector>
#include <set>
#include <algorithm>
#include <ctime>
using namespace std;

Degree::Degree() {
	env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 4);

	//initialize bi
	int row = 31, col = 32;
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			bi[i][j] = 0;
		}
	}
	bi[0][0] = 1;
	for (int i = 1; i < row; i++) {
		bi[i][0] = 1;
		for (int j = 1; j < i + 1; j++) {
			bi[i][j] = bi[i - 1][j - 1] + bi[i - 1][j];
		}
	}

	candidates.clear();
}

int Degree::constructV(vector<int> &k, vector<bool>& flag, int r, int n) {
	//V_{r,s}^{R} = {(\sum b_ik_i)%n |\sum b_i=r-1}
	//flag is used to mark the possible positions for v
	flag.clear();
	for (int i = 0; i < n; i++)
		flag.push_back(0);

	//compute V_{r,s}^{R}
	bool isFull = false;
	vector<int> b;
	b.clear();
	int s = k.size();
	computeV(flag, k,r,n,isFull,0, r - 1, s, b);

	//output V_{r,s}^{R}
	int nonzero = 0;
	//cout << "the set V_{r,s}^{R}:";
	for (int i = 0; i < n; i++) {
		if (flag[i] == 1) {
			//cout << i << " ";
			nonzero++;
		}
	}
	//cout << "size:"<<nonzero<<endl;

	b.clear();
	return nonzero;
}

void Degree::computeV(vector<bool>& flag, vector<int>& k, int r, int n, 
	bool &isFull, int curLayer, int left, int& totalLayer, vector<int>& b) {
	if (curLayer == totalLayer-1) {
		b.push_back(left);

		int sum = 0;
		int pos = 0;
		for (int i = 0; i < totalLayer; i++) {
			//cout << b[i] << " ";
			sum += b[i];
			pos = (pos + k[i] * b[i]) % n;
		}
		//cout << endl;
		if (sum != r-1) {
			cout << "program errors" << endl;
		}
		flag[pos] = 1;
		b.pop_back();

		isFull = true;
		for (int i = 0; i < n; i++) {
			if (flag[i] == false) {
				isFull = false;
				break;
			}
		}
	}
	else {
		for (int i = 0; i <= left; i++) {
			b.push_back(i);
			computeV(flag, k,r,n, isFull, curLayer + 1, left - i, totalLayer, b);
			b.pop_back();
			if (isFull) {
				break;
			}
		}
	}
}

int Degree::constructVectorA(vector<vector<bool> >& A, vector<int> &Z, vector<int>& k, int r, int n,int d) {
	vector<bool> flag;
	flag.clear();
	constructV(k, flag, r, n);

	//construct A
	for (int i = 0; i <=r; i++) {
		for (int j = 0; j < n; j++) {
			A[i][(j + (r - i) * d) % n] = flag[j];
		}
	}

	//construct Z: all possible positions
	Z.clear();
	for (int j = 0; j < n; j++) {
		bool isFind = false;
		for (int i = 0; i <= r; i++) {
			if (A[i][j]) {
				isFind = true;
				break;
			}
		}
		if (isFind) {
			Z.push_back(j);
		}
	}

	flag.clear();
	return Z.size();
}

int Degree::isExponentialIncrease(vector<int>& k, int r, int n, int d) {
	vector<vector<bool> > A(r + 1);
	for (int i = 0; i<A.size(); i++)
		A[i].resize(n);
	vector<int> Z;
	Z.clear();
	constructVectorA(A, Z, k, r, n, d);

	//initial necessary conditions
	int e = log2(n);
	int max = n;
	if (e >= r)
		max = 1 << r;
	if (Z.size() < max) {
		//cout << "The condition on |Z| does not hold: |Z|="<<Z.size() << endl;
		return -1;
	}
	int nonzero = 0;
	for (int j = 0; j < n; j++) {
		if (A[0][j] == 1)
			nonzero++;
	}

	if (max< n && nonzero < bi[r][r / 2]) {
		//cout << "The condition on |V_{r,s}^R| does not hold: |V_{r,s}^R|="<<nonzero << endl;
		return -1;
	}

	int res = buildModelUnivariate(A, r, n,d);
	//cout << "solution:" << res << endl;
	if (res == max) {
		for (int i = 0; i < k.size(); i++)
			cout << k[i] << " ";
		cout << ": ";
		cout << "The maximal degree " << max << " is reached!   ";
		outputTime();
		
		/*candidates.push_back(k);
		if (candidates.size() == 80) {
			outputCandidate();
			outputTime();
			system("pause");
		}*/
		return 1;
	}

	return 0;//the necessary condition does not hold
}

void Degree::initializeVector(GRBModel& model, vector<GRBVar>& v, int cols) {
	v.clear();
	v.resize(cols);
	for (int i = 0; i < cols; i++)
		v[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
}
void Degree::initializeVector(GRBModel& model, vector<vector<GRBVar> >& v, int rows, int cols) {
	v.clear();
	v.resize(rows);
	for (int i = 0; i < rows; i++){
		v[i].resize(cols);
		for (int j = 0; j < cols; j++)
			v[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
	}
}

int Degree::buildModelUnivariate(vector<vector<bool> >& A, int r, int n,int d) {
	clock_t  t;
	t = clock();
	
	vector<vector<GRBVar> > up, dw;
	vector<GRBVar> X;

	//GRBEnv env;
	//env.set(GRB_IntParam_OutputFlag, 0);
	//env.set(GRB_IntParam_Threads, 4);
	GRBModel model = GRBModel(env);

	initializeVector(model, X, n);
	initializeVector(model, up, r + 1, n);
	initializeVector(model, dw, r + 1, n);

	//up relation between up[i+1] and dw[i]
	for (int i = 1; i <r+1; i++) {
		for (int j = 0; j < n; j++) {
			up[i][j] = dw[i - 1][(j + d) % n];
		}
	}

	//on up[0] and dw[r]
	for (int i = 0; i < n; i++) {
		model.addConstr(up[0][i] == 0);
		model.addConstr(dw[r][i] == 0);
	}

	//on positions that cannot be chosen
	for (int i = 0; i < r + 1; i++) {
		for (int j = 0; j < n; j++) {
			if (A[i][j] == 0) {
				model.addConstr(up[i][j] == 0);
				model.addConstr(dw[i][j] == 0);
			}
		}
	}

	/*//on the mutual choice of up[i][j] and dw[i][j]
	for (int i = 0; i < r + 1; i++) {
		for (int j = 0; j < n; j++) {
			model.addConstr(up[i][j] + dw[i][j] <= 1);
		}
	}*/

	//on the total number of choices
	for (int i = 0; i < r; i++) {
		GRBLinExpr sum = 0;
		for (int j = 0; j < n; j++)
			sum = sum + dw[i][j];
		model.addConstr(sum <= bi[r - 1][i]);//we start the index from 0 (in the paper, we start from 1)
	}

	//on X
	for (int i = 0; i < n; i++) {
		GRBLinExpr sum = 0;
		for (int j = 0; j < r + 1; j++) {
			model.addConstr(X[i] >= up[j][i]);
			model.addConstr(X[i] >= dw[j][i]);
			//model.addConstr(X[i] >= up[j][i] + dw[j][i]);
			sum = sum + up[j][i] + dw[j][i];
		}
		model.addConstr(X[i] <= sum);
	}

	GRBLinExpr obj = 0;
	for (int i = 0; i < n; i++) {
		obj = obj + X[i];
	}

	model.setObjective(obj, GRB_MAXIMIZE);
	model.optimize();

	t = clock() - t;

	if (model.get(GRB_IntAttr_Status) == 3) {
		return -2;
	}

	float ft = ((float)t) / CLOCKS_PER_SEC;
	if (ft > max_time)
		max_time = ft;
	if (ft < min_time)
		min_time = ft;
	//cout << "time:" << ((float)t) / CLOCKS_PER_SEC << "seconds" << endl;

	return model.getObjective().getValue();
}

void Degree::searchParameters(int r, int n, int d, int s) {
	vector<int> k;
	k.clear();
	int validNum = 0;
	k.push_back(0);//k_0 = 0
	iterativeEnu(k, 1, s, r, n, d, validNum);
	cout << "total valid candidates:" << validNum << endl;
}

void Degree::iterativeEnu(vector<int>& k, int curLayer, int totalLayer, int r, int n, int d, int& validNum) {
	if (curLayer == totalLayer) {
		//outputVec(k);
		bool isEqui=checkEquivalence(k, n);
		if (!isEqui) {
			//outputVec(k);
			
			int isValid = isExponentialIncrease(k, r, n, d);
			//cout << endl;
			if (isValid == 1) {
				validNum++;
			}
		}
	}
	else {
		int start = k[k.size()-1] + 1;
		for (int i = start; i < n; i++) {
			k.push_back(i);
			iterativeEnu(k, curLayer + 1, totalLayer, r, n, d, validNum);
			k.pop_back();
		}
	}
}

bool Degree::checkEquivalence(vector<int>& k, int n) {
	vector<int> eq(k.size());
	//cout << endl;
	//outputVec(k);
	for (int i = 1; i < k.size(); i++) {
		for (int j = 0; j < k.size(); j++) {
			eq[j] = (k[j] - k[i] + n) % n;
		}
		sort(eq.begin(), eq.end());
		bool isFind = true;
		for (int i = 1; i < k.size(); i++) {
			if (eq[i] == k[i]) {
				continue;
			}
			else if (eq[i] > k[i]) {
				isFind = false;
				break;
			}
			else {
				break;
			}
		}
		if (isFind)
			return true;
		//outputVec(eq);
	}
	return false;
}







/*
*compute useful upper bounds for attackers
*/

void Degree::modAddition(GRBModel& model, vector<GRBVar>& a, int sa, vector<GRBVar>& b, int sb, vector<GRBVar>& c, vector<GRBVar>& c1, vector<GRBVar>& q, vector<GRBVar>& d, int n) {
	model.addConstr(c[0] == 0);
	for (int i = 0; i < n; i++) {
		model.addConstr(a[(i - sa + n) % n] + b[(i - sb + n) % n] + c[i] == q[i] + 2 * c[i + 1]);
	}
	model.addConstr(c1[0] == c[n]);
	for (int i = 0; i < n; i++) {
		model.addConstr(2 * c1[i + 1] + d[i] == q[i] + c1[i]);
	}
}

void Degree::addition(GRBModel& model, vector<GRBVar>& a, vector<GRBVar>& b, vector<GRBVar>& c, vector<GRBVar>& s, int la,int sumLen) {
	model.addConstr(c[0] == 0);
	/*if (n < sumLen) {
		cout << "errors! Please modify addition()!" << endl;
	}*/
	for (int i = 0; i < la; i++) {
		model.addConstr(a[i] + b[i] + c[i] == s[i] + 2 * c[i + 1]);
	}
	for (int i = la; i < sumLen; i++) {
		model.addConstr(b[i] + c[i] == s[i] + 2 * c[i + 1]);
	}
}

int Degree::upperBoundUnivariate(int r, int n, int d, vector<int>& k) {
	vector<vector<bool> > A(r + 1);
	for (int i = 0; i < A.size(); i++)
		A[i].resize(n);
	vector<int> Z;
	Z.clear();
	constructVectorA(A, Z, k, r, n, d);

	int pi = buildModelUnivariate(A, r, n, d);//the number of possible nonzero positions in N_r
	//cout << "pi:" << pi << endl;
	//cout << "Z.size:" << Z.size() << endl;

	int res = solveGeneralOP(Z, pi, r, n);
	return res;
}

int Degree::solveGeneralOP(vector<int>& Z, int pi, int r, int n) {
	clock_t  t;
	t = clock();
	
	GRBModel model = GRBModel(env);
	int k = Z.size();

	int sumLen = ceil(log2f(n));
	sumLen = sumLen + r + 1;//we need to sum up at most n variables smaller than 2^r.

	vector<vector<GRBVar> > a;
	vector<vector<GRBVar> > ms;//modular sum
	vector<vector<GRBVar> > cs;//common sum
	vector<vector<GRBVar> > mcarry;//carry for the modular addition
	vector<vector<GRBVar> > ccarry;//carry for the common addition
	vector<GRBVar> u;//use to denote the choice of pi positions from Z

	initializeVector(model, a, k, n);
	initializeVector(model,ms, k * 2 + 1, n);
	initializeVector(model,cs, k + 1, sumLen);
	initializeVector(model,mcarry, k * 2, n + 1);
	initializeVector(model,ccarry, k, sumLen + 1);
	initializeVector(model, u, k);

	for (int i = n - 1; i >= r + 1; i--) {
		for (int j = 0; j < k; j++) {
			model.addConstr(a[j][i] == 0);
		}
	}

	for (int i = 0; i < n; i++) {
		model.addConstr(ms[0][i] == 0);
	}

	for (int i = 0; i < sumLen; i++) {
		model.addConstr(cs[0][i] == 0);
	}

	for (int i = 0; i < k; i++) {
		modAddition(model, ms[i * 2], 0, a[i], Z[i], mcarry[i * 2], mcarry[i * 2 + 1], ms[i * 2 + 1], ms[i * 2 + 2], n);
		//addition(model, cs[i], a[i], ccarry[i], cs[i + 1], r+1, sumLen);
		addition(model, a[i], cs[i], ccarry[i], cs[i + 1], r + 1, sumLen);
	}

	//extra condition on the common addition sum
	for (int i = sumLen-1; i >= r + 1; i--) {
		model.addConstr(cs[k][i] == 0);
	}
	for (int i = 0; i < r; i++) {
		model.addConstr(cs[k][i] + cs[k][r] <= 1);
	}

	//extra condition on the choice
	GRBLinExpr uSum = 0;
	for (int i = 0; i < k; i++) {
		uSum = uSum + u[i];
		GRBLinExpr sum = 0;
		for (int j = 0; j < n; j++) {
			sum = sum + a[i][j];
			model.addConstr(u[i] >= a[i][j]);
		}
		model.addConstr(u[i] <= sum);
	}
	model.addConstr(uSum <= pi);

	GRBLinExpr obj = 0;
	for (int i = 0; i < n; i++) {
		obj = obj + ms[k*2][i];
	}

	model.setObjective(obj, GRB_MAXIMIZE);
	model.optimize();

	t = clock() - t;

	float ft = ((float)t) / CLOCKS_PER_SEC;
	if (ft > max_time)
		max_time = ft;
	if (ft < min_time)
		min_time = ft;

	int res = model.getObjective().getValue();

	return res;
}











void Degree::outputVec(vector<int>& vec) {
	for (int i = 0; i < vec.size(); i++)
		cout << vec[i] << " ";
	cout << endl;
}

void Degree::outputCandidate() {
	for (int i = 0; i < candidates.size(); i++) {
		for (int j = 0; j < candidates[i].size(); j++) {
			if (j == 0)
				cout << "(";
			cout << candidates[i][j];
			if (j < candidates[i].size() - 1)
				cout << ",";
			else if (j == candidates[i].size() - 1)
				cout << "),";
		}
		if ((i + 1) % 8 == 0) {
			cout << endl;
		}
	}
}

void Degree::outputTime() {
	cout << "min_time:" << min_time << "seconds;       ";
	cout << "max_time:" << max_time << "seconds" << endl;
}

