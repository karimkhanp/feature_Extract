#pragma once
#include <vector>
using namespace std;


class ButterFilter
{
public:
	ButterFilter();
	ButterFilter(int kind);
	~ButterFilter();
public:
	int len;
	double a[2];
	double b[2];
	double zi;

public:
	void FiltData(const vector<double> sig, vector<double>& sout, const int data_len);
	void BiDirectionalFilter(const vector<double> sig, vector<double>& sout, const int data_len);
	void signalExtend(const vector<double> sig, vector<double>& sout, const int data_len, const int nfact);
	void signalShrink(const vector<double> sig, vector<double>& sout, const int data_len, const int nfact);
	void dataRollingOver(const vector<double> sig, vector<double> &sout, const int data_len);
};





