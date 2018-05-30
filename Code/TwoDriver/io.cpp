#include "shared.h"
#include <iostream>
#include <sstream>
#include <string>

void GetRecombinationMap(int chr, vector<int>& pos, vector<double>& r) {
	ifstream rec_file;
	stringstream ss;
	ss << chr;
	string str = ss.str();
	string name;
	name="../../Data/Rec"+str+".dat";
	rec_file.open(name.c_str());
	double x;
	int n;
	for (int i=0;i<10000;i++) {
		if (!(rec_file >> n)) break;
		pos.push_back(n);
		if (!(rec_file >> x)) break;
		r.push_back(x);
	}
}

void GetLocusList(int chr, vector<int>& loci) {
	ifstream loc_file;
	stringstream ss;
	ss << chr;
	string str = ss.str();
	string name;
	name="../../Data/Loci"+str+".in";
        cout << name << "\n";
	loc_file.open(name.c_str());
	int n;
	for (int i=0;i<100000;i++) {
		if (!(loc_file >> n)) break;
		loci.push_back(n);
	}
}

void GetData (int chr, data& d1, data& d2, data& d3, data& d4) {
	ifstream dat_file;
	stringstream ss;
	ss << chr;
	string str = ss.str();
	string name;
	name="../../Data/04_BC4_CAVR_pre-IVM-P2.sync_poly50_95_filter3.snp"+str+".dat";
	dat_file.open(name.c_str());
	int n;
	for (int i=0;i<100000;i++) {
		if (!(dat_file >> n)) break;
		if (!(dat_file >> n)) break;
		d1.loc.push_back(n);
		if (!(dat_file >> n)) break;
		d1.num.push_back(n);
		if (!(dat_file >> n)) break;
		d1.tot.push_back(n);
	}
	dat_file.close();
	name="../../Data/06_BC4_CAVR_post-IVM.sync_poly50_95_filter3.snp"+str+".dat";
        dat_file.open(name.c_str());
	for (int i=0;i<100000;i++) {
		if (!(dat_file >> n)) break;
		if (!(dat_file >> n)) break;
		d2.loc.push_back(n);
		if (!(dat_file >> n)) break;
		d2.num.push_back(n);
		if (!(dat_file >> n)) break;
		d2.tot.push_back(n);
	}
	dat_file.close();
	name="../../Data/08_BC4-passage_CAVR_post-IVM-pass3.sync_poly50_95_filter3.snp"+str+".dat";
        dat_file.open(name.c_str());
	for (int i=0;i<100000;i++) {
		if (!(dat_file >> n)) break;
		if (!(dat_file >> n)) break;
		d3.loc.push_back(n);
		if (!(dat_file >> n)) break;
		d3.num.push_back(n);
		if (!(dat_file >> n)) break;
		d3.tot.push_back(n);
	}
	dat_file.close();
	name="../../Data/09_BC4-passage_CAVR_post-IVM-pass4.sync_poly50_95_filter3.snp"+str+".dat";
        dat_file.open(name.c_str());
	for (int i=0;i<100000;i++) {
		if (!(dat_file >> n)) break;
		if (!(dat_file >> n)) break;
		d4.loc.push_back(n);
		if (!(dat_file >> n)) break;
		d4.num.push_back(n);
		if (!(dat_file >> n)) break;
		d4.tot.push_back(n);
	}
	dat_file.close();

}

