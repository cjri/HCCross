//Program to calculate predicted sigma values for loci in a model system, where data is generated by a sampling process
//This code uses a deterministic method to run the underlying population distribution, and samples it at low frequency

#include <iostream>
#include <vector>
#include <list>
#include <deque>
#include <sstream>
using namespace std;

#include "shared.h"

/*#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>*/


vector<haplo> population;
vector<haplo> population_back;

//To go in: bottleneck sizes

int main(int argc, char *argv[]){
	
	//Set up random number generator
	int seed=atoi(argv[1]);
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rgen, seed);
	double c=78.0;
	
	//Need to check dataset - should be pattern across all chromosomes / datasets.
	
	//Initial parameters
	int chr=atoi(argv[2]);
	int n_offspring=100; /*Number of offspring produced by a mating event*/
	int loc1=atoi(argv[3]); /*First locus under selection*/
	int loc2=atoi(argv[4]); /*Second locus under selection*/
	//mu=0.00001;
	double s0=1;
	double s1=atof(argv[5]); /*First selection coefficient*/
	double s2=atof(argv[6]); /*Second selection coefficient*/
	double sh1=0.5*s1;
	double sh2=0.5*s2;
	int reps=atoi(argv[7]);
	//Import recombination rates
	vector<int> pos;
	vector<double> r;
	//Import recombination map
	GetRecombinationMap(chr,pos,r);
	
	//Import locus list
	vector<int> loci;
	GetLocusList(chr,loci);
	loc1=loci[loc1];
	loc2=loci[loc2];
	cout << "Selected loci " << loc1 << " " << loc2 << "\n";
	cout << "Selection coefficients " << s0 << " " << sh1 << " " << s1 << " " << sh2 << " " << s2 << "\n";
	
	//Import CAVR data
	data d1;
	data d2;
	data d3;
	data d4;
	GetData (chr,d1,d2,d3,d4);
	
	vector<double> fact_store;
	FindLogFact(fact_store,10000);

	double best_log=-1e10;
	vector<double> log_listh;
	vector<double> log_list1;
	
	vector<double> log_store;
	
	for (int it=0;it<=0;it++) {
		
		log_store.clear();
		
		for (int rep=0;rep<reps;rep++) {
			string st="Sim";
			stringstream ss;
			ss << rep;
			string str = ss.str();
			string st2="mkdir "+st+str;
			const char* cc=st2.c_str();
			system (cc);

			population.clear();
			population_back.clear();
		
			//Set up population: 100 of each L4 genotype
			for (int i=0;i<50;i++) {
				haplo h;
				population.push_back(h);
			}
			for (int i=0;i<50;i++) {
				haplo h;
				h.seq1.push_back(-1);  //Mark only the recombination points
				h.seq2.push_back(-1);  //Mark only the recombination points
				population_back.push_back(h);
			}
	
			//Step 1: cross in sheep and collect eggs.  Eight sets of pairings.  10^6 eggs.
	
			//Do recombination
	
			//Do backcross
			vector<haplo> new_population;

			//Produce eggs: Eight sets of mating pairs via backcross - get eggs
			RunBackcross(n_offspring,pos,r,population,population_back,rgen);
	
			//Grow L3 worms from eggs - select 5000 or correct number
			SimpleBottleneck(5000,population,rgen);
	
			//Grow worms in sheep under selection.
			GrowL3toL4WormsTD (66,1,loc1,loc2,s1,s2,population,rgen);

	
		//PrintPopulation(population);
	

	
		//cout << "Population size " << population.size() << " worms\n";
	
			//cout << "Output into 1\n";
		//	PrintPopulation(population);

	
			//Backcross: start of loop - get eggs
			RunBackcross(n_offspring,pos,r,population,population_back,rgen);

			//Grow L3 worms from eggs - select 5000 or correct number
			SimpleBottleneck(8000,population,rgen);
	
			//Grow worms in sheep under selection.
			GrowL3toL4WormsTD (45,1,loc1,loc2,s1,s2,population,rgen);

			//cout << "Output into 2\n";
		//	PrintPopulation(population);

	
			//Backcross: second time - get eggs
			RunBackcross(n_offspring,pos,r,population,population_back,rgen);
	
			//Grow L3 worms from eggs - select 5000 or correct number
			SimpleBottleneck(10000,population,rgen);
	
			//Grow worms in sheep under selection.
			GrowL3toL4WormsTD (100,1,loc1,loc2,s1,s2,population,rgen);

				//cout << "Output into 3\n";
			//	PrintPopulation(population);

	
			//Backcross: third time - get eggs
			RunBackcross(n_offspring,pos,r,population,population_back,rgen);
	
			//Grow L3 worms from eggs - select 5000 or correct number
			SimpleBottleneck(10000,population,rgen);
	
			//Grow worms in sheep under selection.
			//cout << "First\n";
			GrowL3toL4WormsTD (90,1,loc1,loc2,s1,s2,population,rgen);

			//cout << "\n";
				//cout << "Output into 4\n";
			//	PrintPopulation(population);

	
			//Backcross: fourth time - get eggs
			RunBackcross(n_offspring,pos,r,population,population_back,rgen);
	
			//Grow L3 worms from eggs - select 5000 or correct number
			SimpleBottleneck(10000,population,rgen);
	
	
			//Output of first loop
			vector<haplo> population1=population;
	
		//	cout << "Output from Loop\n";
		//	PrintPopulation(population);

	
		//First pathway

			//Grow eggs in sheep not under selection.
			GrowL3toL4WormsTD (10000,0,loc1,loc2,s1,s2,population,rgen);

			//Cross against self: divide into male and female
			SelfCross(n_offspring,pos,r,population,rgen);
		
		//Will be just the same as the backcross code, but divide the population randomly into two.
		//Might need fewer offspring per pair to get sensible numbers here.
	
		//To go here...
	
	
			//Grow L3 worms from eggs - select 5000 or correct number
			SimpleBottleneck(10000,population,rgen);

			vector<haplo> population_store=population;
			//No drug
			//cout << "Next\n";
			GrowL3toL4WormsTD (50,0,loc1,loc2,s1,s2,population,rgen);

			//cout << "\n";
				//cout << "Output 1\n";
				string name="Data1.out";
				PrintPopulation(name,population);
	
			//Get allele frequencies
			vector<double> afs1;
			CalculateAF(d1.loc,population,afs1);
			
			population=population_store;
			GrowL3toL4WormsTD (50,1,loc1,loc2,s1,s2,population,rgen);

			//cout << "Output 2\n";
			name="Data2.out";
			PrintPopulation(name,population);

			//Get allele frequencies
			vector<double> afs2;
			CalculateAF(d2.loc,population,afs2);

	
		//Second pathway
		
			population=population1;

			//Grow eggs in sheep under selection.
			GrowL3toL4WormsTD (10000,1,loc1,loc2,s1,s2,population,rgen);

			//Cross against self: divide into male and female
			SelfCross(n_offspring,pos,r,population,rgen);

			//Grow L3 worms from eggs - select 5000 or correct number
			SimpleBottleneck(5000,population,rgen);

			//Second self-cross
			
			//Grow eggs in sheep under selection.
			GrowL3toL4WormsTD (10000,1,loc1,loc2,s1,s2,population,rgen);

			//Cross against self: divide into male and female
			SelfCross(n_offspring,pos,r,population,rgen);

			//Grow L3 worms from eggs - select 5000 or correct number
			SimpleBottleneck(5000,population,rgen);

			population1=population;
	
			//Third self-cross
			//Grow eggs in sheep under selection.
			GrowL3toL4WormsTD (50,1,loc1,loc2,s1,s2,population,rgen);

		//		cout << "Output 3\n";
				name="Data3.out";
				PrintPopulation(name,population);
	
			//Get allele frequencies
			vector<double> afs3;
			CalculateAF(d3.loc,population,afs3);

			population=population1;
	
			//Evaluate third cross to eggs
			GrowL3toL4WormsTD (10000,1,loc1,loc2,s1,s2,population,rgen);

			SelfCross(n_offspring,pos,r,population,rgen);

			//Grow L3 worms from eggs - select 5000 or correct number
			SimpleBottleneck(5000,population,rgen);

			//Fourth self-cross
			//Grow eggs in sheep under selection.
			//cout << "Last\n";
			GrowL3toL4WormsTD (50,1,loc1,loc2,s1,s2,population,rgen);

	
		//cout << "Output 4\n";
			name="Data4.out";
			PrintPopulation(name,population);

			//Get allele frequencies
			vector<double> afs4;
			CalculateAF(d4.loc,population,afs4);
		
		//		cout << "Output freqs\n";
		//Calculate likelihood
			double logL=0;
			ofstream d1_file;
			ofstream d2_file;
			ofstream d3_file;
			ofstream d4_file;
			d1_file.open("Freqs1.out");
			d2_file.open("Freqs2.out");
			d3_file.open("Freqs3.out");
			d4_file.open("Freqs4.out");
			
			string str1="mv Freqs1.out "+st+str;
			string str2="mv Freqs2.out "+st+str;
			string str3="mv Freqs3.out "+st+str;
			string str4="mv Freqs4.out "+st+str;
			//cout << str1 << "\n";
			const char* c1=str1.c_str();
			const char* c2=str2.c_str();
			const char* c3=str3.c_str();
			const char* c4=str4.c_str();
			system (c1);
			system (c2);
			system (c3);
			system (c4);


			for (int i=0;i<d1.loc.size();i++){
				d1_file << i << " " << d1.loc[i] << " " << d1.num[i] << " " << d1.tot[i] << " " << afs1[i] << " ";
				double l=BetBinomCalcNew(c,d1.tot[i],d1.num[i],afs1[i],fact_store);
				d1_file << l << "\n";
				logL=logL+l;
			}
			for (int i=0;i<d2.loc.size();i++){
				d2_file << i << " " << d2.loc[i] << " " << d2.num[i] << " " << d2.tot[i] << " " << afs2[i] << " ";
				double l=BetBinomCalcNew(c,d2.tot[i],d2.num[i],afs2[i],fact_store);
				d2_file << l << "\n";
				logL=logL+l;
			}
			for (int i=0;i<d3.loc.size();i++){
				d3_file << i << " " << d3.loc[i] << " " << d3.num[i] << " " << d3.tot[i] << " " << afs3[i] << " ";
				double l=BetBinomCalcNew(c,d3.tot[i],d3.num[i],afs3[i],fact_store);
				d3_file << l << "\n";
				logL=logL+l;
			}
			for (int i=0;i<d4.loc.size();i++){
				d4_file << i << " " << d4.loc[i] << " " << d4.num[i] << " " << d4.tot[i] << " " << afs4[i] << " ";
				double l=BetBinomCalcNew(c,d4.tot[i],d4.num[i],afs4[i],fact_store);
				d4_file << l << "\n";
				logL=logL+l;
			}
			ofstream log_file;
			log_file.open("Log.out");
			log_file << "Log L " << logL << "\n";
			log_store.push_back(logL);
			log_file.close();
			str4="mv Log.out "+st+str;
			c4=str4.c_str();
			system (c4);
		}
		

		double all_log=0;
		for (int i=0;i<log_store.size();i++) {
			all_log=all_log+log_store[i];
		}
		all_log=all_log/reps;
		cout << "Mean log " << all_log << "\n";
	}
	
	return 0;
}
	
	
