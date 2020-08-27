/* NucEnergy.cpp
 * Author: Rohan Shah
 */

#include <getopt.h>
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <regex>
#include <map>
#include <math.h>
#include <vector>
#include <numeric>

using namespace std;

#define PI 3.1419265
typedef vector<vector<vector<double>>> weightVec;

map<char, int> baseToIndex;
weightVec weights;


string cleanSeq(string seq)
{
	transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
	replace(seq.begin(), seq.end(), 'U', 'T');
	regex e ("[^GATC]");
	regex_replace(seq, e, "");
	return seq;
}

weightVec weightCalc(int w, double b, double p)
{
	vector<int> pos(w);
	int start = 0;
	iota(pos.begin(), pos.end(), start);

	vector<double> bpl(w), bpl3(w), bmi(w), bmi3(w), bnul(w);
	transform(pos.begin(), pos.end(), bpl.begin(), [=](int i){return (0.25 + b * sin(2*PI*i/p));});
	transform(pos.begin(), pos.end(), bpl3.begin(), [=](int i){return (0.25 + b/3 * sin(2*PI*i/p));});
	transform(pos.begin(), pos.end(), bmi.begin(), [=](int i){return (0.25 - b * sin(2*PI*i/p));});
	transform(pos.begin(), pos.end(), bmi3.begin(), [=](int i){return (0.25 - b/3 * sin(2*PI*i/p));});
	fill(bnul.begin(), bnul.end(), 0.25);

	vector<double> aa = bpl, ac = bmi3, ag = bmi3, at = bmi3;
	vector<double> ca = bnul, cc = bnul, cg = bnul, ct = bnul;
	vector<double> ga = bpl3, gc = bmi, gg = bpl3, gt = bpl3;
	vector<double> ta = bpl, tc = bmi, tg = bmi, tt = bpl;

	vector< vector<double> > a = {aa, ac, ag, at};
	vector< vector<double> > c = {ca, cc, cg, ct};
	vector< vector<double> > g = {ga, gc, gg, gt};
	vector< vector<double> > t = {ta, tc, tg, tt};

	vector< vector< vector<double> > > acgt = {a, c, g, t};

	return acgt;
}

double indexE(double pf, double pr)
{
    return -(pf * log(pf) + pr * log(pr)) / (pr + pf);
}

vector<double> calcEnergy(string seq, int w, double b, double p)
{
	vector<double> pf, pr;

	for (unsigned int i = 0; i < seq.length() - w; i++)
	{
		double psf = 1, psr = 1;

		for (int s = 0; s < w; s++)
		{
			int baseSense = baseToIndex[seq.at(i+s)];
			int nextBaseSense = baseToIndex[seq.at(i+s+1)];
			int baseAntiSense = 3-baseToIndex[seq.at(i+w-s)];
			int nextBaseAntiSense = 3-baseToIndex[seq.at(i+w-s-1)];

			psf *= weights[baseSense][nextBaseSense][s]*4;
			psr *= weights[baseAntiSense][nextBaseAntiSense][s]*4;
		}

		pf.push_back(psf);
		pr.push_back(psr);
	}

	pr.push_back(pr[0]);
	pr.erase(pr.begin());

	vector<double> energyVec(pf.size());
	transform(pf.begin(), pf.end(), pr.begin(), energyVec.begin(), indexE);

	return energyVec;
}

void fragEnergy(string input_name, string output_name, int w, double b, double p)
{
	ifstream input_file(input_name);
	ofstream output_file(output_name);

	string chr, seq;
	int start, stop;

	while (input_file >> chr >> start >> stop >> seq)
	{
		string seqclean = cleanSeq(seq);
		vector<double> energyVec = calcEnergy(seqclean, w, b, p);
		auto minElementIter = min_element(energyVec.begin(), energyVec.end());
		double minEnergy = *minElementIter;
		int minElementIndex = minElementIter - energyVec.begin();

		output_file << chr << "\t" << start << "\t" << stop << "\t" << minEnergy << "\t" << minElementIndex << "\n";
	}

	input_file.close();
	output_file.close();
}

int main(int argc, char **argv)
{
	string input_name, output_name;
	int w = 74, opt;
	double b = 0.2, p = 10;

	while ((opt = getopt(argc, argv, "i:o:w:b:p:")) != -1)
	{
		switch (opt)
		{
		case 'i':
			input_name = optarg;
			break;
		case 'o':
			output_name = optarg;
			break;
		case 'w':
			w = stoi(optarg);
			break;
		case 'b':
			b = stod(optarg);
			break;
		case 'p':
			p = stod(optarg);
			break;
		}
	}

	baseToIndex['A'] = 0; baseToIndex['C'] = 1; baseToIndex['G'] = 2; baseToIndex['T'] = 3;
	weights = weightCalc(w, b, p);
	fragEnergy(input_name, output_name, w, b, p);

	return 0;
}
