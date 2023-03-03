

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
// #include <windows.h>
#include <string.h>
#include <vector>
using namespace std;
#include "ECGEngine.h"

vector<double> test_target(const char* input_file_path)
{
	vector<double> bufVec;
	bufVec.reserve(10000);
	FILE* fp = fopen(input_file_path, "r");
	
	if (fp)
	{
		double temp;
		while (1 == fscanf(fp,"%lf,",&temp))
		{
			bufVec.push_back(temp);
		}
		
	}
	else
	{
		printf("Error  *Input file opening error*. There is no?");
	}

	return bufVec;
}
int main(int argc, char** argv)
{
	vector<double> inBuf;
	inBuf = test_target("100_85_N.csv");
	double skew_val = 0, iqr_val = 0, pw_skew = 0, pw_skew_mean = 0;
	int pw_skew_abs = 0;
	qrs_skew_q(inBuf, skew_val, iqr_val);
	pwave_attr(inBuf, 300, pw_skew, pw_skew_abs, pw_skew_mean);
	printf("pw_skew:%f, pw_skew_abs:%d, pw_skew_mean:%f, sk: %f, iqr:%f\n", pw_skew, pw_skew_abs, pw_skew_mean, skew_val, iqr_val);
	return 0;
}