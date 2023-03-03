#include "butterFilter.h"
#include <iostream>



void ButterFilter::FiltData(const vector<double> sig, vector<double>& sout, const int data_len)
{
    double scale = sig[0];
	for (int i = 0; i < data_len; i++)
	{
		double tmp = 0;
		tmp += b[0] * sig[i];
		for (int j = 1; j < len; j++)
		{
            if(i - j < 0){ 
                if (zi){
                    tmp += zi * scale;
                }
                break;
            }
			tmp += b[j] * sig[i - j];
			tmp -= a[j] * sout[i - j];
		}
		sout[i] = tmp;
	}
}

// CHANGED
void ButterFilter::BiDirectionalFilter(const vector<double> sig, vector<double> &sout, const int data_len)
{
	int nfact = 3 * len;
	if (nfact>data_len)
		nfact = data_len - 1;
    int data_len_ext = data_len + 2*nfact;
	vector<double> s_ext(data_len_ext);
	vector<double> tmp1(data_len_ext);
	vector<double> tmp2(data_len_ext);
	
    signalExtend(sig, s_ext, data_len, nfact);
	FiltData(s_ext, tmp1, data_len_ext);
	dataRollingOver(tmp1, tmp2, data_len_ext);
	FiltData(tmp2, tmp1, data_len_ext);
	dataRollingOver(tmp1, tmp2, data_len_ext);
    signalShrink(tmp2, sout, data_len, nfact);


}

// CHANGED
void ButterFilter::signalExtend(const vector<double> sig, vector<double>& sout, const int data_len, const int nfact)
{
	for (int i = 0; i < nfact; i++)
		sout[i] = 2 * sig[0] - sig[nfact - i];

	for (int i = nfact; i > 0; i--)
		sout[data_len+nfact+i-1] = 2 * sig[data_len-1] - sig[data_len-1-i];

    for (int i = nfact; i < nfact + data_len; i++)
        sout[i] = sig[i-nfact];
}

void ButterFilter::signalShrink(const vector<double> sig, vector<double>& sout, const int data_len, const int nfact)
{
	for (int i = 0; i < data_len; i++)
        sout[i] = sig[i + nfact];
}

void ButterFilter::dataRollingOver(const vector<double> sig, vector<double>& sout, const int data_len)
{
	for (int i = 0; i < data_len; i++)
		sout[data_len - i - 1] = sig[i];
}
ButterFilter::ButterFilter()
{
}

ButterFilter::ButterFilter(int kind)
{
	len = 2;
	if (kind == 1)
	{//low pass filter
		a[0] = 1, a[1] = -0.72654253;
		b[0] = 0.13672874, b[1] = 0.13672874;
		zi = 0.86327126;
	}
	else
	{//high pass filter
		a[0] = 1, a[1] = -0.90040404;
		b[0] = 0.95020202, b[1] = -0.95020202;
		zi = -0.95020202;
	}
}

ButterFilter::~ButterFilter()
{
}
