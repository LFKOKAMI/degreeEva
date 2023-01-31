#ifndef _DEGREE_H_
#define _DEGREE_H_

#include <vector>
#include "gurobi_c++.h"
using namespace std;

class Degree {
	int bi[31][32];//we only precompute a table for binomial{a}{b} where a<=30
	vector<vector<int> > candidates;
public:
	GRBEnv env;
	float max_time = 0;
	float min_time = 100000;
	Degree();
	void outputCandidate();
	void searchParameters(int r, int n, int d, int s);
	void iterativeEnu(vector<int>& k, int curLayer, int totalLayer, int r, int n, int d, int& validNum);
	int constructV(vector<int> &k, vector<bool>& flag,int r,int n);
	void computeV(vector<bool>& flag, vector<int>& k, int r, int n, bool& isFull, int curLayer, int left, int& totalLayer, vector<int>& b);
	int constructVectorA(vector<vector<bool> >&A, vector<int>& Z, vector<int>& k, int r, int n,int d);
	int isExponentialIncrease(vector<int>& k, int r, int n, int d);
	int buildModelUnivariate(vector<vector<bool> >& A, int r, int n,int d);
	void initializeVector(GRBModel& model, vector<GRBVar>& v, int cols);
	void initializeVector(GRBModel& model, vector<vector<GRBVar> >& v, int rows, int cols);

	//upper bounding the degree
	void modAddition(GRBModel& model, vector<GRBVar>& a, int sa, vector<GRBVar>& b, int sb, vector<GRBVar>& c, vector<GRBVar>& c1, vector<GRBVar>& q, vector<GRBVar>& d, int n);
	void addition(GRBModel& model, vector<GRBVar>& a, vector<GRBVar>& b, vector<GRBVar>& c, vector<GRBVar>& s, int n,int sumLen);
	int upperBoundUnivariate(int r, int n, int d, vector<int>& k);
	int solveGeneralOP(vector<int>& Z, int pi, int r, int n);

	void outputVec(vector<int>& vec);
	bool checkEquivalence(vector<int>& k,int n);

	void outputTime();
};

#endif
