#include "shared.h"
#include <iostream>
#include <string>

void MakePermutation (int n, vector<int>& p, gsl_rng *rgen) {
	p.clear();
	vector<int> a;
	for(int i=0;i<n;i++) {
		a.push_back(i);
		p.push_back(-1);
	}
	for (int i=0;i<n;i++) {
		int r=floor(gsl_rng_uniform(rgen)*a.size());
		p[i]=a[r];
		a.erase(a.begin()+r);
	}
}

void MakePartPermutation (int n, int m, vector<int>& p, gsl_rng *rgen) {
	//cout << "MPP " << n << " " << m << "\n";
	p.clear();
	int verb=0;
	int add=0;
	if (m>0&&n<m+1) {
		int diff=m-n+1;
		for (int i=0;i<diff;i++) {
			int r=floor(gsl_rng_uniform(rgen)*n);
			p.push_back(r);
			add++;
		}
		m=n-1;
	//	cout << "Now " << n << " " << m << "\n";
		verb=1;

	}
	vector<int> a;
	for (int i=0;i<n;i++) {
		a.push_back(i);
	}
	for (int i=0;i<m;i++) {
		p.push_back(-1);
	}
	for (int i=0;i<m;i++) {
		int r=floor(gsl_rng_uniform(rgen)*a.size());
	//	if (verb==1) {
	//		cout << "i= " << i << " " << i+add << " r " << r << " " << a.size() << " " << p.size() << "\n";
	//	}
		p[i+add]=a[r];
		a.erase(a.begin()+r);
	}
}

void RunBackcross (int n_offspring, vector<int> pos, vector<double> r, vector<haplo>& population, vector<haplo>& population_back, gsl_rng *rgen) {
	while (population_back.size()<population.size()) {
		population_back.push_back(population_back[0]);
	}
//	cout << "RunBackcross " << population.size() << " " << population_back.size() << " " << n_offspring << "\n";
	vector<haplo> new_population;
	vector<int> rec_points;
	for (int p=0;p<8;p++) {
		vector<int> b;
		//Random pairs of worms
		MakePermutation(population.size(),b,rgen);
		for (int i=0;i<population.size();i++) {
			for (int j=0;j<n_offspring;j++) {
				//cout << "j " << j << " " << n_offspring << "\n";
				FindRecombinationPoints(r,pos,rec_points,rgen);
				EvaluateRecombination(0,population,population_back,i,b[i],rec_points,new_population,rgen);
			}
		}
	}
	//cout << "Size " << new_population.size() << "\n";
	population=new_population;
//	cout << "End sub\n";
}

void FindRecombinationPoints(vector<double>& r, vector<int>& pos, vector<int>& rec_points, gsl_rng *rgen) {
//	cout << "FindRecombinationPoints\n";
	rec_points.clear();
	for (int i=0;i<r.size();i++) {
		int nr=gsl_ran_poisson(rgen,r[i]);
		for (int j=0;j<nr;j++) {
			int rp=gsl_rng_uniform(rgen)*pos[i];
			rec_points.push_back(rp);
			//cout << "Recombination at " << rp << "\n";
		}
	}
//	cout << "End sub\n";
}

void EvaluateRecombination(int verb, vector<haplo>& population1, vector<haplo>& population2, int in1, int in2, vector<int>& rec_points, vector<haplo>& new_population, gsl_rng *rgen) {
	int r1=floor(gsl_rng_uniform(rgen)*2);
	int r2=floor(gsl_rng_uniform(rgen)*2);
	
	vector<int> seq1;
	vector<int> seq2;
	
	if (r1==1) {
		seq1=population1[in1].seq1;
	} else {
		seq1=population1[in1].seq2;
	}

	if (r2==1) {
		seq2=population2[in2].seq1;
	} else {
		seq2=population2[in2].seq2;
	}

	/*cout << "Before\n";
	for (int i=0;i<seq1.size();i++) {
		cout << seq1[i] << " ";
	}
	cout << "\n";
	for (int i=0;i<seq2.size();i++) {
		cout << seq2[i] << " ";
	}
	cout << "\n";*/
	
	if (rec_points.size()>0&&verb==1) {
		cout << "Before\n";
		 for (int i=0;i<seq1.size();i++) {
		 cout << seq1[i] << " ";
		 }
		 cout << "\n";
		 for (int i=0;i<seq2.size();i++) {
		 cout << seq2[i] << " ";
		 }
		 cout << "\n";

		cout << "recombination points\n";
		for (int i=0;i<rec_points.size();i++) {
			cout << rec_points[i] << " ";
		}
		cout << "\n";
	}
	
	for (int i=0;i<rec_points.size();i++) {
		if (verb==1) {
			cout << "Recombination at " << rec_points[i] << "\n";
		}
		vector<int> temp1;
		vector<int> temp2;
		//Find first position
		int pos1=-1;
		for (int j=0;j<seq1.size();j++) {
			if (seq1[j]<rec_points[i]) {
				pos1=j;
			}
		}
		int pos2=-1;
		for (int j=0;j<seq2.size();j++) {
			if (seq2[j]<rec_points[i]) {
				pos2=j;
			}
		}
		if (verb==1) {
			cout << "Pos " << pos1 << " " << pos2 << "\n";
		}
		//If pos1 and pos2 are negative, all sites are after the new recombination point.

		//Do up to point of recombination
		int gen1=1;
		int gen2=1;
		for (int j=0;j<=pos2;j++) {
			temp1.push_back(seq2[j]);
			gen1=1-gen1;
		}
		
		for (int j=0;j<=pos1;j++) {
			temp2.push_back(seq1[j]);
			gen2=1-gen2;
		}
		
		if (verb==1) {
			cout << "Gen " << gen1 << " " << gen2 << "\n";
		}
		
		//Add recombination point if this constitutes a flip in the genotype.  Genotype at this point given by gen1 and gen2.  Neglect recombination at same point twice.
		if (gen1!=gen2) {
			temp1.push_back(rec_points[i]);
			temp2.push_back(rec_points[i]);
		}

		//After point of recombination
		for (int j=pos1+1;j<seq1.size();j++) {
			temp1.push_back(seq1[j]);
		}
		for (int j=pos2+1;j<seq2.size();j++) {
			temp2.push_back(seq2[j]);
		}
		seq1=temp2;
		seq2=temp1;
		
		if (verb==1) {
		cout << "After\n";
		for (int i=0;i<seq1.size();i++) {
		 cout << seq1[i] << " ";
		}
		cout << "\n";
		for (int i=0;i<seq2.size();i++) {
		 cout << seq2[i] << " ";
		}
		cout << "\n";
		}

	}
	
	if (rec_points.size()>0&&verb==1) {
		
		int loc=12064674;
		int var1=1;
		int var2=1;
		for (int j=0;j<seq1.size();j++) {
			if (seq1[j]<loc) {
				var1=1-var1;
			}
			if (seq1[j]>loc) {
				break;
			}
		}
		for (int j=0;j<seq2.size();j++) {
			if (seq2[j]<loc) {
				var2=1-var2;
			}
			if (seq2[j]>loc) {
				break;
			}
		}
		if (var1==1) {
			if (var2==1) {
				cout << "Type 1\n";
			} else {
				cout << "Type s\n";
			}
		} else {
			if (var2==1) {
				cout << "Type s\n";
			} else {
				cout << "Type 0\n";
			}

		}
		cout << "\n";

	}

	
	haplo h;
	h.seq1=seq1;
	h.seq2=seq2;
	h.fitness=0;
	new_population.push_back(h);
}

void SimpleBottleneck (int n, vector<haplo>& population, gsl_rng *rgen) {
//	cout << "SimpleBottleneck\n";
	//Note: Requires n<population.size();
	//cout << "Bottleneck to " << n << "\n";
	vector<int> b;
	vector<haplo> new_pop;
	MakePartPermutation(population.size(),n,b,rgen);
	for (int i=0;i<b.size();i++) {
		new_pop.push_back(population[b[i]]);
	}
	population=new_pop;
//	cout << "End sub\n";
}

void GrowL3toL4Worms (int n, int sel, int loc, double s1, double sh, double s0, vector<haplo>& population, gsl_rng *rgen) {
//	cout << "L3 to L4 worms\n";
	vector<haplo> new_population;
	vector<haplo> pop0;
	vector<haplo> poph;
	vector<haplo> pop1;
	//Divide population by genotype at given locus
	for (int i=0;i<population.size();i++) {
		int var1=1;
		int var2=1;
		for (int j=0;j<population[i].seq1.size();j++) {
			if (population[i].seq1[j]<loc) {
				var1=1-var1;
			}
			if (population[i].seq1[j]>loc) {
				break;
			}
		}
		for (int j=0;j<population[i].seq2.size();j++) {
			if (population[i].seq2[j]<loc) {
				var2=1-var2;
			}
			if (population[i].seq2[j]>loc) {
				break;
			}
		}
		if (var1==0) {
			if (var2==0) {
				population[i].fitness=s0;
				pop0.push_back(population[i]);
			} else {
				population[i].fitness=sh;
				poph.push_back(population[i]);
			}
		} else {
			if (var2==0) {
				population[i].fitness=sh;
				poph.push_back(population[i]);
			} else {
				population[i].fitness=s1;
				pop1.push_back(population[i]);
			}
		}
	}

	//cout << "N1 " << pop1.size() << " Nh " << poph.size() << " N0 " << pop0.size() << "\n";
	double w1=1;
	double wh=1;
	double w0=1;
	if (sel==1) {
		//cout << "Do selection " << s0 << " " << sh << " " << s1 << "\n";
		//cout << "Sizes " << pop1.size() << " " << poph.size() << " " << pop0.size() << "\n";
		w1=pop1.size()*s1;
		wh=poph.size()*sh;
		w0=pop0.size()*s0;
	} else {
		w1=pop1.size();
		wh=poph.size();
		w0=pop0.size();
	}
	double w_tot=w1+wh+w0;
	//cout << "W1 " << w1 << " Wh " << wh<< " W0 " << w0 << " " << w_tot << "\n";;

	//Multinomial draw: select number of individuals from population given fitness of each class
	double prob_select_mirror[3];
	prob_select_mirror[0]=w1/w_tot;
	prob_select_mirror[1]=wh/w_tot;
	prob_select_mirror[2]=w0/w_tot;
	
	//cout << "Mirror 1 " << prob_select_mirror[0] << " h " << prob_select_mirror[1] << " 0 " << prob_select_mirror[2] << "\n";
	
	unsigned int next_gen[3];
	gsl_ran_multinomial(rgen,3,n,prob_select_mirror,next_gen);
	
	//cout << "N1 " << next_gen[0] << " Nh " << next_gen[1] << " N0 " << next_gen[2] << "\n";
	
	//cout << pop1.size() << " " << poph.size() << " " << pop0.size() << "\n";
	
	//Construct new population with individuals according to the genotypic composition of the new population
	vector<int> b;
	MakePartPermutation(pop1.size(),next_gen[0],b,rgen);
	for (int i=0;i<next_gen[0];i++) {
		new_population.push_back(pop1[b[i]]);
	}
	b.clear();
	MakePartPermutation(poph.size(),next_gen[1],b,rgen);
	for (int i=0;i<next_gen[1];i++) {
		new_population.push_back(poph[b[i]]);
	}
	b.clear();
	MakePartPermutation(pop0.size(),next_gen[2],b,rgen);
	for (int i=0;i<next_gen[2];i++) {
		new_population.push_back(pop0[b[i]]);
	}

	population=new_population;
//	cout << "End sub\n";
}

void GrowL3toL4WormsTD (int n, int sel, int loc1, int loc2, double s1, double s2, vector<haplo>& population, gsl_rng *rgen) {
	//	cout << "L3 to L4 worms\n";
	vector<haplo> new_population;
	vector<haplo> pop0;
	vector<haplo> pop1;
	vector<haplo> pop2;
	vector<haplo> pop11;
	vector<haplo> pop12;
	vector<haplo> pop22;
	vector<haplo> pop112;
	vector<haplo> pop122;
	vector<haplo> pop1122;
	//Divide population by genotype at given locus
	for (int i=0;i<population.size();i++) {
		int var11=1;
		int var12=1;
		int var21=1;
		int var22=1;
		for (int j=0;j<population[i].seq1.size();j++) {
			if (population[i].seq1[j]<loc1) {
				var11=1-var11;
			}
			if (population[i].seq1[j]>loc1) {
				break;
			}
		}
		for (int j=0;j<population[i].seq2.size();j++) {
			if (population[i].seq2[j]<loc1) {
				var21=1-var21;
			}
			if (population[i].seq2[j]>loc1) {
				break;
			}
		}
		for (int j=0;j<population[i].seq1.size();j++) {
			if (population[i].seq1[j]<loc2) {
				var12=1-var12;
			}
			if (population[i].seq1[j]>loc2) {
				break;
			}
		}
		for (int j=0;j<population[i].seq2.size();j++) {
			if (population[i].seq2[j]<loc2) {
				var22=1-var22;
			}
			if (population[i].seq2[j]>loc2) {
				break;
			}
		}

		
		population[i].fitness=1;
		if (var11+var12==2) {
			if (var21+var22==2) {
				population[i].fitness=1+s1+s2;
				pop1122.push_back(population[i]);
			} else if (var21+var22==1) {
				population[i].fitness=1+s1+(0.5*s2);
				pop112.push_back(population[i]);
			} else {
				population[i].fitness=1+s1;
				pop11.push_back(population[i]);
			}
		} else if (var11+var12==1) {
			if (var21+var22==2) {
				population[i].fitness=1+(0.5*s1)+s2;
				pop122.push_back(population[i]);
			} else if (var21+var22==1) {
				population[i].fitness=1+(0.5*s1)+(0.5*s2);
				pop12.push_back(population[i]);
			} else {
				population[i].fitness=1+(0.5*s1);
				pop1.push_back(population[i]);
			}
		} else {
			if (var21+var22==2) {
				population[i].fitness=1+s2;
				pop22.push_back(population[i]);
			} else if (var21+var22==1) {
				population[i].fitness=1+(0.5*s2);
				pop2.push_back(population[i]);
			} else {
				population[i].fitness=1;
				pop0.push_back(population[i]);
			}
		}
	}
	
	//cout << "N " << pop0.size() << " " << pop1.size() << " " << pop2.size() << " " << pop11.size() << " " << pop12.size() << " " << pop22.size() << " " << pop112.size() << " " << pop122.size() << " " << pop1122.size() << "\n";
	
	double w0=1;
	double w1=1;
	double w2=1;
	double w11=1;
	double w12=1;
	double w22=1;
	double w112=1;
	double w122=1;
	double w1122=1;
	if (sel==1) {
		//cout << "Do selection " << s0 << " " << sh << " " << s1 << "\n";
		//cout << "Sizes " << pop1.size() << " " << poph.size() << " " << pop0.size() << "\n";
		w0=pop0.size();
		w1=pop1.size()*(1+(0.5*s1));
		w2=pop2.size()*(1+(0.5*s2));
		w11=pop11.size()*(1+s1);
		w12=pop12.size()*(1+(0.5*s1)+(0.5*s2));
		w22=pop22.size()*(1+s2);
		w112=pop112.size()*(1+s1+(0.5*s2));
		w122=pop122.size()*(1+(0.5*s1)+s2);
		w1122=pop1122.size()*(1+s1+s2);
	} else {
		w0=pop0.size();
		w1=pop1.size();
		w2=pop2.size();
		w11=pop11.size();
		w12=pop12.size();
		w22=pop22.size();
		w112=pop112.size();
		w122=pop122.size();
		w1122=pop1122.size();
	}
	double w_tot=w0+w1+w2+w11+w12+w22+w112+w122+w1122;
//	cout << w0 << " " << w1 << " " << w2 << " " << w11 << " " << w12 << " " << w22 << " " << w112 << " " << w122 << " " << w1122 << "\n";
	
	//Multinomial draw: select number of individuals from population given fitness of each class
	double prob_select_mirror[9];
	prob_select_mirror[0]=w0/w_tot;
	prob_select_mirror[1]=w1/w_tot;
	prob_select_mirror[2]=w2/w_tot;
	prob_select_mirror[3]=w11/w_tot;
	prob_select_mirror[4]=w12/w_tot;
	prob_select_mirror[5]=w22/w_tot;
	prob_select_mirror[6]=w112/w_tot;
	prob_select_mirror[7]=w122/w_tot;
	prob_select_mirror[8]=w1122/w_tot;
	
	//cout << "Mirror 1 " << prob_select_mirror[0] << " h " << prob_select_mirror[1] << " 0 " << prob_select_mirror[2] << "\n";
	
	unsigned int next_gen[9];
	gsl_ran_multinomial(rgen,9,n,prob_select_mirror,next_gen);
	
	//cout << "N1 " << next_gen[0] << " Nh " << next_gen[1] << " N0 " << next_gen[2] << "\n";
	
	//cout << pop1.size() << " " << poph.size() << " " << pop0.size() << "\n";
	
	//Construct new population with individuals according to the genotypic composition of the new population
	vector<int> b;
	MakePartPermutation(pop0.size(),next_gen[0],b,rgen);
	for (int i=0;i<next_gen[0];i++) {
		new_population.push_back(pop0[b[i]]);
	}
	b.clear();
	MakePartPermutation(pop1.size(),next_gen[1],b,rgen);
	for (int i=0;i<next_gen[1];i++) {
		new_population.push_back(pop1[b[i]]);
	}
	b.clear();
	MakePartPermutation(pop2.size(),next_gen[2],b,rgen);
	for (int i=0;i<next_gen[2];i++) {
		new_population.push_back(pop2[b[i]]);
	}
	b.clear();
	MakePartPermutation(pop11.size(),next_gen[3],b,rgen);
	for (int i=0;i<next_gen[3];i++) {
		new_population.push_back(pop11[b[i]]);
	}
	b.clear();
	MakePartPermutation(pop12.size(),next_gen[4],b,rgen);
	for (int i=0;i<next_gen[4];i++) {
		new_population.push_back(pop12[b[i]]);
	}
	b.clear();
	MakePartPermutation(pop22.size(),next_gen[5],b,rgen);
	for (int i=0;i<next_gen[5];i++) {
		new_population.push_back(pop22[b[i]]);
	}
	b.clear();
	MakePartPermutation(pop112.size(),next_gen[6],b,rgen);
	for (int i=0;i<next_gen[6];i++) {
		new_population.push_back(pop112[b[i]]);
	}
	b.clear();
	MakePartPermutation(pop122.size(),next_gen[7],b,rgen);
	for (int i=0;i<next_gen[7];i++) {
		new_population.push_back(pop122[b[i]]);
	}
	b.clear();
	MakePartPermutation(pop1122.size(),next_gen[8],b,rgen);
	for (int i=0;i<next_gen[8];i++) {
		new_population.push_back(pop1122[b[i]]);
	}
	b.clear();
	population=new_population;
	//	cout << "End sub\n";
}



void SelfCross (int n_offspring, vector<int>& pos, vector<double>& r, vector<haplo>& population, gsl_rng *rgen) {
//	cout << "Self cross\n";
	vector<haplo> population_m;
	vector<haplo> population_f;
	int half=population.size()/2;
	for (int i=0;i<half;i++) {
		population_m.push_back(population[i]);
	}
	for (int i=half;i<population.size();i++) {
		population_f.push_back(population[i]);
	}
	RunBackcross(n_offspring/10,pos,r,population_m,population_f,rgen);
	population=population_m;
//	cout << "End sub\n";
}

void PrintPopulation(string s, vector<haplo> population) {
	ofstream out_file;
	out_file.open(s.c_str());
	for (int i=0;i<population.size();i++) {
	//	out_file << i << "\n";
		for (int j=0;j<population[i].seq1.size();j++) {
			out_file << population[i].seq1[j] << " ";
		}
		out_file << "\n";
		for (int j=0;j<population[i].seq2.size();j++) {
			out_file << population[i].seq2[j] << " ";
		}
		out_file << "\n";
	}
	out_file.close();
}

void CalculateAF (vector<int>& loci, vector<haplo>& population, vector<double>& afs) {
//	cout << "CalculateAF\n";
	for (int i=0;i<loci.size();i++) {
		afs.push_back(0);
	}
	double n=2*population.size();
	for (int i=0;i<population.size();i++) {
		double a=1;
		int index=0;
		for (int j=0;j<population[i].seq1.size();j++) {
			double x=a/n;
			while (loci[index]<population[i].seq1[j]) {
				afs[index]=afs[index]+x;
				index++;
			}
			a=1-a;
			//cout << "Index " << i << " 1 " << index << " " << loci[index] << "\n";
		}
		double x=a/n;
		while (index<afs.size()) {
			afs[index]=afs[index]+x;
			index++;
		}
		
		a=1;
		index=0;
		for (int j=0;j<population[i].seq2.size();j++) {
			double x=a/n;
			while (loci[index]<population[i].seq2[j]) {
				afs[index]=afs[index]+x;
				index++;
			}
			a=1-a;
			//cout << "Index " << i << " 2 " << index-1 << " " << loci[index-1] << "\n";
		}
		x=a/n;
		while (index<afs.size()) {
			afs[index]=afs[index]+x;
			index++;
		}
	}
//	cout << "End sub\n";
}

void FindLogFact(vector<double>& fact_store,int N){
	double logN=0;
	fact_store.push_back(0);
	for (int i=1;i<=N;i++) {
		logN=logN+log(i);
		fact_store.push_back(logN);
	}
}


double BetBinomCalcNew (double c, int n, int a, double x, vector<double>& fact_store) {
	vector<int> obs;
	vector<double> inf;
	double bin=0;
	if (x==0) {
		x=0.00001;
	}
	if (n>0) {
		bin=fact_store[n];
		bin=bin-fact_store[a];
		int na=n-a;
		bin=bin-fact_store[na];
		double alpha=c*x;
		double beta=c*(1-x);
		double ab=alpha+beta;
		
		bin=bin+gsl_sf_lngamma(a+alpha);
		bin=bin+gsl_sf_lngamma(na+beta);
		bin=bin-gsl_sf_lngamma(n+alpha+beta);

		
		bin=bin+gsl_sf_lngamma(ab);
		bin=bin-gsl_sf_lngamma(alpha);
		bin=bin-gsl_sf_lngamma(beta);

	} else {
		bin=0;
	}
	//	cout << "L " << bin << "\n";
	//	cout << "\n";
	
	return(bin);
}

