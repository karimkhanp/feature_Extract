#pragma once
#include <vector>
#include <math.h>
#include <iomanip> 
#include <algorithm> 
#include "Skew.h"
#include "butterFilter.h"
using namespace std;

vector<double> lgth_transform(vector<double> ecgData, int ws)
{
	int lgth = ecgData.size();
	vector<double> diff(lgth);
	diff.assign(lgth, 0);
	//pad 'edge' ws
	vector<double> temp1(ws);
	vector<double> temp2(ws);
	temp1.assign(ws, ecgData.front());
	temp2.assign(ws, ecgData.back());	
	ecgData.insert(ecgData.begin(), temp1.begin(), temp1.end());
	ecgData.insert(ecgData.end(), temp2.begin(), temp2.end());

	for (int i = 0; i < lgth; i++)
	{
		double leftVal = ecgData[i + ws] - ecgData[i];
		double rightVal = ecgData[i + ws] - ecgData[i + ws + ws];
		diff[i] = max(0.0, min(leftVal, rightVal));
		diff[i] = diff[i] * diff[i];
	}
	return diff;
}


vector<double> integrate(vector<double> ecgData, int ws)
{
	int lgth = ecgData.size();
	vector<double> integrate_ecg(lgth);
	integrate_ecg.assign(lgth, 0);
	// symmetric padding with pad-width 'ws/2'
	int num = (int)ceil(ws / 2);
	vector<double> temp1(num);
	vector<double> temp2(num);
	copy(ecgData.begin(), ecgData.begin() + num, temp1.begin());
	reverse(temp1.begin(), temp1.end());
	copy(ecgData.end() - num, ecgData.end(), temp2.begin());
	reverse(temp2.begin(), temp2.end());
	ecgData.insert(ecgData.begin(), temp1.begin(), temp1.end());
	ecgData.insert(ecgData.end(), temp2.begin(), temp2.end());
	for (int i = 0; i < lgth; i++)
	{
		for (int j = 0; j < ws; j++)
		{
			integrate_ecg[i] += ecgData[i + j];
		}
		integrate_ecg[i] /= ws;
	}
	return integrate_ecg;
}


vector<double> find_peak(vector<double> data, int ws)
{
	int lgth = data.size();
	vector<double> true_peaks;
	for (int i = 0; i < lgth - ws + 1; i++)
	{
		vector<double> temp(ws);
		copy(data.begin() + i, data.begin() + i + ws, temp.begin());
		int index = (int)((ws - 1) / 2);
		bool peak = true;
		for (int j = 0; j < index; j++)
		{
			if ((temp[index - j] <= temp[index - j - 1]) || (temp[index + j] <= temp[index + j + 1]))
			{
				peak = false;break;
			}
		}
		if (peak) {
			true_peaks.push_back((int)(i + 1 + (ws - 1) / 2));
		}
	}
	return true_peaks;
}


vector<int> find_R_peaks(vector<double> ecgData, vector<double> peaks, int ws)
{
	int num_peak = peaks.size();
	int num_ecg = ecgData.size();
	vector<int> R_peaks;
	for (int i = 0; i < num_peak; i++)
	{ 
		double temp = peaks[i];
		double offset = temp - 2 * ws;
		if((offset > 0) && (temp < num_ecg))
		{
			vector<double> temp_ecg(2 * ws);
			copy(ecgData.begin() + offset, ecgData.begin() + temp, temp_ecg.begin());
			int dist = distance(temp_ecg.begin(), max_element(temp_ecg.begin(), temp_ecg.end())) + offset;
			R_peaks.push_back(dist);
		}
	}
	return R_peaks;
}


vector<int> find_S_point(vector<double> ecgData, vector<int> R_peaks)
{
	int num_peak = R_peaks.size();
	int num_ecg = ecgData.size();
	vector<int> S_point;
	for (int i = 0; i < num_peak; i++)
	{
		try
		{
			int cnt = R_peaks[i];
			if (cnt + 1 >= num_ecg)
				break;
			while (ecgData[cnt] > ecgData[cnt + 1])
			{
				cnt += 1;
				if(cnt >= num_ecg) break;
			}
			S_point.push_back(cnt);
		}
		
		catch (const std::range_error&)
		{
			printf("RangeError");
		}
	}
	return S_point;
}


vector<int> fine_Q_point(vector<double> ecgData, vector<int> R_peaks)
{
	int num_peak = R_peaks.size();
	int num_ecg = ecgData.size();
	vector<int> Q_point;
	for (int i = 0; i < num_peak; i++)
	{
		int cnt = R_peaks[i];
		if (cnt - 1 < 0)
			break;
		while (ecgData[cnt] > ecgData[cnt - 1])
		{
			cnt -= 1;
			if (cnt < 0) break;
		}
		Q_point.push_back(cnt);
	}
	return Q_point;
}


void find_qrs(vector<double> ecgData, int fs, bool QS, vector<int>& R_peaks, vector<int>& Q_point, vector<int>& S_point)
{
	double meanVal = 0, sum = 0.0f;
	for (int i = 0; i < ecgData.size(); i++)
	{
		sum += ecgData.at(i);
	}
	meanVal = sum / ecgData.size();
	for (int i = 0; i < ecgData.size(); i++)
	{
		ecgData[i] -= meanVal;
	}
	vector<double> ecg_lgth_transform = lgth_transform(ecgData, (int)(fs / 20));
	vector<double> ecg_integrate;
	int ws = (int)(fs / 8);
	ecg_integrate = integrate(ecg_lgth_transform, ws);
	for (int i = 0; i < ecg_integrate.size(); i++)
		ecg_integrate[i] /= ws;
	ws = (int)(fs / 6);
	ecg_integrate = integrate(ecg_integrate, ws);
	ws = (int)(fs / 36);
	ecg_integrate = integrate(ecg_integrate, ws);
	ws = (int)(fs / 72);
	ecg_integrate = integrate(ecg_integrate, ws);
	vector<double> peaks = find_peak(ecg_integrate, (int)(fs / 10));
	R_peaks = find_R_peaks(ecgData, peaks, (int)(fs / 40));
	
	if (QS)
	{
		S_point = find_S_point(ecgData, R_peaks);
		Q_point = fine_Q_point(ecgData, R_peaks);
	}
	return;
}


void swap(int a, int b, vector<double> &nums) {
	double temp = nums[a];
	nums[a] = nums[b];
	nums[b] = temp;
	return;
}

void sortArray(int n, vector<double> &nums) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n - 1; j++) {
			if (nums[i] < nums[j])
				swap(i, j, nums);
		}
	}
	return;
}

double getMedian(int n, vector<double> nums) {
	if (n % 2 == 0)
		return (double)(nums[n / 2] + nums[(n / 2) - 1]) / 2;
	else
		return nums[n / 2];
}

double iqr(vector<double> nums){
	
	int	n = nums.size();
	sortArray(n, nums);

	vector<double> lowerNums, higherNums;
	for (int i = 0; i < n; i++) {
		if (n % 2 == 1) {
			if (i < n / 2)
				lowerNums.push_back(nums[i]);
			else if (i > n / 2)
				higherNums.push_back(nums[i]);
		}
		else {
			if (i <= (n / 2) - 1)
				lowerNums.push_back(nums[i]);
			else
				higherNums.push_back(nums[i]);
		}
	}

	double firstQuartile = getMedian(lowerNums.size(), lowerNums);
	double thirdQuartile = getMedian(higherNums.size(), higherNums);

	return thirdQuartile - firstQuartile;
}
void qrs_skew_q(vector<double> ecgData, double& iqrValue, double& skewValue)
{
	vector<int> R_peaks(100), S_point(100), Q_point(100);
	find_qrs(ecgData, 320, true, R_peaks, Q_point, S_point);
	vector<double> diff_list;
	if (R_peaks.size() >= S_point.size())
	{
		for (int i = 0; i < R_peaks.size(); i++)
		{
			if (S_point[i] <= R_peaks[i]) continue;
			vector<double> temp(S_point[i] - R_peaks[i]);
			copy(ecgData.begin() + R_peaks[i], ecgData.begin() + S_point[i], temp.begin());
			double maxValue = *max_element(temp.begin(), temp.end());
			double minValue = *min_element(temp.begin(), temp.end());
			diff_list.push_back(maxValue - minValue);
		}
	}
	double sum, mean, var, dev,skew, kurt;
	computeStats(diff_list.begin(), diff_list.end(), sum, mean, var, dev, skew, kurt);
	iqrValue = iqr(diff_list);
	skewValue = skew;
	return;
}


vector<double> correlate1d(vector<double> a, vector<double> b)
{
	vector<double> corre_vec(a.size());
	if (a.size() < b.size())
	{
		vector<double> temp = a;
		a = b;
		b = temp;
	}

	int len = b.size();
	int lenOrigin = a.size();
	int num1 = len / 2, num2 = len / 2;
	if (len % 2 == 0) num2 -= 1;
	vector<double> temp1(num1);
	vector<double> temp2(num2);
	copy(a.begin(), a.begin() + num1, temp1.begin());
	reverse(temp1.begin(), temp1.end());
	copy(a.end() - num2, a.end(), temp2.begin());
	reverse(temp2.begin(), temp2.end());
	a.insert(a.begin(), temp1.begin(), temp1.end());
	a.insert(a.end(), temp2.begin(), temp2.end());

	for (int i = 0; i < lenOrigin; i++)
	{
		for (int j = 0; j < b.size(); j++)
		{
			corre_vec[i] += a[i + j] * b[j];
		}
	}
	return corre_vec;
}
vector<double> gaussian_filter(vector<double> sig, double sigma, double truncate = 4.0)
{
	vector<double> so(sig.size());
	int lw = (int)(truncate * sigma + 0.5);
	vector<double> weights(lw * 2 + 1);
	double sum = 0;
	for (int i = -lw; i <= lw; i++)
	{
		weights[i + lw] = exp(-0.5 / (sigma * sigma) * i * i);
		sum += weights[i + lw];
	}
	for (int i = 0; i < weights.size(); i++)
	{
		weights[i] /= sum;
	}
	reverse(weights.begin(), weights.end());
	so = correlate1d(sig, weights);
	return so;
}
double average(const vector<double>& data, const int size)
{
	int i;
	double sum = 0, mean;
	for (i = 0; i < size; i++)
	{
		sum = sum + data[i];
	}
	mean = sum / size;
	return mean;
}

double stddev(const vector<double>& data, const int size)
{
	double sum = 0, mean, std0;
	int i;
	mean = average(data, size);
	for (i = 0; i < size; i++)
		sum += pow(data[i] - mean, 2);
	std0 = pow(sum / (size - 1), 0.5);
	return std0;
}
vector<double> shannon_energy(vector<double> a, vector<double> b)
{
	if (a.size() < b.size())
	{
		vector<double> temp = a;
		a = b;
		b = temp;
	}
	vector<double> shannon(a.size());

	int len = b.size();
	int lenOrigin = a.size();

	a.insert(a.begin(), len >> 1, 0);

	if (len % 2 == 0) a.insert(a.end(), (len >> 1) - 1, 0);
	else a.insert(a.end(), len >> 1, 0);
	
	for (int i = 0; i < lenOrigin; i++)
	{
		for (int j = 0; j < b.size(); j++)
		{
			shannon[i] += a[i + j] * b[len - j - 1];
		}
	}
	return shannon;
}

vector<int> peak_detection(vector<double> sig, double rate, double ransac_window_size = 5.0, double lowfreq = 5.0, double highfreq = 15.0)
{
	vector <int> peaks;

	peaks.clear();
	ransac_window_size = (int)(ransac_window_size * rate);

	ButterFilter butter1(1), butter2(2);
	vector<double> tmp1(sig.size()), tmp(sig.size());

	butter1.BiDirectionalFilter(sig, tmp1, sig.size());
	butter2.BiDirectionalFilter(tmp1, tmp, sig.size());

	vector<double> decg_power(sig.size() - 1);

	for (int i = 0; i < tmp.size() - 1; i++)
	{
		decg_power[i] = tmp[i + 1] - tmp[i];
		decg_power[i] *= decg_power[i];
	}

	vector<double> thresholds;
	vector<double> mx_powers;

	for (int i = 0; i < (int)(decg_power.size() / ransac_window_size); i++)
	{
		vector<double> sample(ransac_window_size);
		copy(decg_power.begin() + i * ransac_window_size, decg_power.begin() + (i + 1) * ransac_window_size, sample.begin());
		thresholds.push_back(0.5 * stddev(sample, sample.size()));
		mx_powers.push_back(*max_element(sample.begin(), sample.end()));
	}

	double threshold, max_power;
	threshold = getMedian(thresholds.size(), thresholds);
	max_power = getMedian(mx_powers.size(), mx_powers);

	vector<double> squ_dec_power(sig.size() - 1), shan_energy(sig.size() - 1);

	for (int i = 0; i < decg_power.size(); i++)
	{
		if (decg_power[i] < threshold) decg_power[i] = 0;
		else decg_power[i] /= max_power;

		if (decg_power[i] > 1.0) decg_power[i] = 1.0;
		squ_dec_power[i] = decg_power[i] * decg_power[i];

		if(squ_dec_power[i] > 0) shan_energy[i] = -squ_dec_power[i] * log(squ_dec_power[i]);
		else shan_energy[i] = 0;
	}
	
	int mean_window_len = (int)(rate * 0.125 + 1);

	vector<double> mean_vec(mean_window_len);
	mean_vec.assign(mean_window_len, 1.0 / mean_window_len);

	vector<double> lp_energy;
	lp_energy = shannon_energy(shan_energy, mean_vec);
	lp_energy = gaussian_filter(lp_energy, rate / 8.0);

	vector<double> lp_energy_diff(lp_energy.size() - 1);

	for (int i = 0; i < lp_energy.size() - 1; i++)
		lp_energy_diff[i] = lp_energy[i + 1] - lp_energy[i];

	for (int i = 0; i < lp_energy_diff.size() - 1; i++)
	{
		if (lp_energy_diff[i] > 0 && lp_energy_diff[i + 1] < 0)
			peaks.push_back(i - 1);
	}

	return peaks;
}

void pwave_attr(vector<double> ecg_arr, double fs, double& skew, int& leng, double& skew_mean)
{
	vector<int> peaks = peak_detection(ecg_arr, fs);
	vector<int> RR(peaks.size() - 1),lpwp(peaks.size() - 1);

	for (int i = 0; i < peaks.size() - 1; i++)
	{
		RR[i] = peaks[i + 1] - peaks[i];
		lpwp[i] = peaks[i] + 0.71 * RR[i];
	}

	vector<int> rpwp(peaks.size() - 1);
	vector<double> pwave_arr;

	for (int i = 1; i < peaks.size() ; i++)
	{
		rpwp[i-1] = peaks[i] - 0.07 * RR[i - 1] -11;
	}

	vector<double> maz;

	for (int i = 0;i < lpwp.size(); i++)
	{
		vector<double> temp(rpwp[i] - lpwp[i]);

		
		for (int k = lpwp[i]; k < rpwp[i]; k++)
		{
			temp[k - lpwp[i]] = fabs(ecg_arr[k]);
		}
		
		maz.push_back(*max_element(temp.begin(), temp.end()));
		
	}

	leng = maz.size();

	double sum, mean, var, dev, kurt;

	computeStats(maz.begin(), maz.end(), sum, mean, var, dev, skew, kurt);

	double mean_maz = 0;

	for (int i = 0; i < maz.size(); i++)
	{
		mean_maz += maz[i];
	}

	mean_maz /= maz.size();

	double max_maz = *max_element(maz.begin(), maz.end());

	skew_mean = (mean_maz / max_maz) * 100;
	return;
}
