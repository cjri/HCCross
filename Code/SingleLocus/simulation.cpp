//Simulation of H. contortus experiment at one locus; Correct setup now accounted for.

#include <iostream>
#include <vector>
#include <map>
using namespace std;

#include "shared.h"

int main(int argc, const char **argv){
	
	//Pre-calcualtions
	vector<double> logfactorial;
	find_logfact(logfactorial);
	
    //Selection terms for het and hom
	double s1=1; //Homozygous resistant
	double sh=1; //Heterozygous resistant
    double s0=1; //Homozygous susceptible
    
	//This code treats the worm stage deterministically due to there being a huge number of eggs.
	
	int dep=200; //Sequencing depth
	int grid=20;  //Internal approximation grid size
	double cut=1e-20; //Minimum value considered in approximate distributions
	
	int type=atoi(argv[1]); //WRS=0, CAVR=1  N.B. Specified on command line.  Determines numbers of L4 and L3 worms in each generation
	
	//Number of L3 and L4 worms that were collected
	vector<int> nl3;
	vector<int> nl4m; //Male susceptible
	vector<int> nl4f; //Female
	if (type==0) {
		SetupNumbersWRS(nl3,nl4m,nl4f);
	} else {
		SetupNumbersCAVR(nl3,nl4m,nl4f);
	}
	
	//Initial number of resistant alleles in the L4 resistant female worms; worms are diploid; total 50 female worms; hence between 0 and 100
	double init_all=atoi(argv[2]);
	
	//Proportion of alleles in resistant with resistance allele
	double prob=init_all/100;
	
	//Mean proportion of eggs produced of each genotype.  Fraction wh are heterozygote, w0 are susceptible.
	double wh=prob;
	double w0=(1-prob);

	//Distribution of eggs in neutral sheep via deterministic formula N.B. Could add drift here
	cout << "Producing eggs in naive sheep\n";
	cout << wh << " " << w0 << "\n";
	dvec eggdist;
	eggdist.one.push_back(wh);  //Values here are percentage heterozygosity and probability of this...
	eggdist.two.push_back(1);
	
	//Number of L3 worms from egg proportions; total 5000; how many are heterozygous...
	pop2struc L3num;
	cout << "Generate L3 worms - distribution of number of heterozygotes\n";
	ConstructPop2FromEggs(1,1,nl3[0],grid,cut,eggdist,L3num,logfactorial);

	//L3 worms fed to sheep, undergo selection.  No mating, but model survival of worms under selection
	pop2struc L4num;
	cout << "Produce L4 worms in IVM sheep\n";
	GetL4fromL3 (1,1,nl3[0],nl4f[1],sh,s0,L3num,L4num,logfactorial);
	
	//First backcross
	
	//Produce eggs in neutral sheep via cross.  Get distribution of mean frequencies from L4
	cout << "First backcross: Mean frequency in egg population\n";
	Get2EggFrequencies (1,0,1,nl4f[1],s0,sh,L4num,eggdist); //L4num to eggdist
	
	cout << "Generate L3 worms\n";
	ConstructPop2FromEggs(1,1,nl3[1],grid,cut,eggdist,L3num,logfactorial); //eggdist to L3num

	cout << "Produce L4 worms in IVM sheep\n";
	GetL4fromL3 (1,1,nl3[1],nl4f[2],sh,s0,L3num,L4num,logfactorial); //L3num to L4num

	//Second backcross

	cout << "Second backcross: Mean frequency in egg population\n";
	Get2EggFrequencies (1,0,1,nl4f[2],s0,sh,L4num,eggdist); //L4num to eggdist
	
	cout << "Generate L3 worms\n";
	ConstructPop2FromEggs(1,1,nl3[2],grid,cut,eggdist,L3num,logfactorial); //eggdist to L3num
	
	cout << "Produce L4 worms in IVM sheep\n";
	GetL4fromL3 (1,1,nl3[2],nl4f[3],sh,s0,L3num,L4num,logfactorial); //L3num to L4num

	//Third backcross
	
	cout << "Third backcross: Mean frequency in egg population\n";
	Get2EggFrequencies (1,0,1,nl4f[3],s0,sh,L4num,eggdist); //L4num to eggdist
	
	cout << "Generate L3 worms\n";
	ConstructPop2FromEggs(1,1,nl3[3],grid,cut,eggdist,L3num,logfactorial); //eggdist to L3num
	
	cout << "Produce L4 worms in IVM sheep\n";
	GetL4fromL3 (1,1,nl3[3],nl4f[4],sh,s0,L3num,L4num,logfactorial); //L3num to L4num

	//Fourth backcross
	
	cout << "Fourth backcross: Mean frequency in egg population\n";
	Get2EggFrequencies (1,0,1,nl4f[4],s0,sh,L4num,eggdist); //L4num to eggdist
	
	cout << "Generate L3 worms\n";
	ConstructPop2FromEggs(1,1,nl3[4],grid,cut,eggdist,L3num,logfactorial); //eggdist to L3num

	//Population is now split between drug and no-drug sheep
	
	int nl3x=10000; //Assumed number of L3 worms for this part of the experiment

	cout << "Scale up the L3 population to allow for het resistant.  At this point there aren't any yet\n";
	pop3struc L3num2;
	for (int i=0;i<L3num.one.size();i++) {
		L3num2.one.push_back(0);
		L3num2.two.push_back(L3num.one[i]);
		L3num2.three.push_back(L3num.two[i]);
	}
	pop3struc BCrossOut=L3num2; //Record of population at this point

	//First set of outputs
	
	cout << "Produce eggs in naive sheep\n";  //L3num to eggdist - skip collection of L4 worms
	tvec eggdist2;
	
	Get3EggFrequencies (1,0,1,nl3[4],s0,sh,s1,BCrossOut,eggdist2);

	cout << "Collect L3 worms\n";
	ConstructPop3FromEggs(1,1,nl3x,grid,cut,eggdist2,L3num2,logfactorial);

	pop3struc L3post=L3num2; //Record population going into final split
	
	pop3struc L4num2;

	//Set 1: Growth of L3 to L4 worms in naive sheep
	
	cout << "Grow L3 to L4 worms in naive sheep\n"; //Multinomial draw without selection
	int nl4x=25;
	GetL4fromL3_2 (1,0,nl3x,nl4x,s1,sh,s0,L3post,L4num2,logfactorial); //L3num to L4num

	//Convert these statistics to an allele frequency distribution...
	pop3struc L4num_fin1=L4num2;
	int nl4x_1=nl4x;
	
	//Set 2: Growth of L3 to L4 worms in IVM sheep
	
	cout << "Grow L3 to L4 worms in IVM sheep\n"; //Multinomial draw without selection
	nl4x=50;
	GetL4fromL3_2 (1,1,nl3x,nl4x,s1,sh,s0,L3post,L4num2,logfactorial); //L3num to L4num
	pop3struc L4num_fin2=L4num2;
	int nl4x_2=nl4x;

	//Second set of outputs
	
	//Repeated crossing phase
	cout << "Produce eggs in treated sheep: round 1\n";
	Get3EggFrequencies (1,1,1,nl3[4],s0,sh,s1,BCrossOut,eggdist2);
	cout << "Collect L3 worms\n";
	ConstructPop3FromEggs(1,1,nl3x,grid,cut,eggdist2,L3num2,logfactorial);

	cout << "Produce eggs in treated sheep: round 2\n";
	Get3EggFrequencies (1,1,1,nl3x,s0,sh,s1,L3num2,eggdist2);
	cout << "Collect L3 worms\n";
	ConstructPop3FromEggs(1,1,nl3x,grid,cut,eggdist2,L3num2,logfactorial);
	L3post=L3num2; //Record population at this point
	
	//Set 3: Growth of L3 to L4 worms in IVM sheep
	cout << "Grow L3 to L4 worms in IVM sheep\n"; //Multinomial draw without selection
	nl4x=50;
	GetL4fromL3_2 (1,1,nl3x,nl4x,s1,sh,s0,L3post,L4num2,logfactorial); //L3num to L4num
	pop3struc L4num_fin3=L4num2;
	int nl4x_3=nl4x;

	//Repeated crossing phase
	cout << "Produce eggs in treated sheep: round 1\n";
	Get3EggFrequencies (1,1,1,nl3[4],s0,sh,s1,BCrossOut,eggdist2);
	cout << "Collect L3 worms\n";
	ConstructPop3FromEggs(1,1,nl3x,grid,cut,eggdist2,L3post,logfactorial);

	//Set 3: Growth of L3 to L4 worms in IVM sheep
	cout << "Grow L3 to L4 worms in IVM sheep\n"; //Multinomial draw without selection
	nl4x=50;
	GetL4fromL3_2 (1,1,nl3x,nl4x,s1,sh,s0,L3post,L4num2,logfactorial); //L3num to L4num
	pop3struc L4num_fin4=L4num2;
	int nl4x_4=nl4x;

	
	//Output datasets
	cout << "Allele frequencies 1\n";
	ofstream dat1_file;
	dat1_file.open("BC4P1.dat");
	pop2struc allfreq1;
	AlleleFrequencyDistribution(1,nl4x_1,L4num_fin1,allfreq1);
	//cout << allfreq1.one.size() << "\n";
	for (int i=0;i<allfreq1.one.size();i++) {
		dat1_file << (allfreq1.one[i]+0.)/(2.0*nl4x_1) << " " << allfreq1.two[i] << "\n";
	}
	
	cout << "Allele frequencies 2\n";
	ofstream dat2_file;
	dat2_file.open("BC4P1_IVM.dat");
	pop2struc allfreq2;
	AlleleFrequencyDistribution(1,nl4x_2,L4num_fin2,allfreq2);
	for (int i=0;i<allfreq2.one.size();i++) {
		dat2_file << (allfreq2.one[i]+0.)/(2.0*nl4x_2) << " " << allfreq2.two[i] << "\n";
	}
	
	cout << "Allele frequencies 3\n";
	ofstream dat3_file;
	dat3_file.open("BC4P1_S3.dat");
	pop2struc allfreq3;
	AlleleFrequencyDistribution(1,nl4x_3,L4num_fin3,allfreq3);
	for (int i=0;i<allfreq3.one.size();i++) {
		dat3_file << (allfreq3.one[i]+0.)/(2.0*nl4x_3) << " " << allfreq3.two[i] << "\n";
	}
	
	cout << "Allele frequencies 4\n";
	ofstream dat4_file;
	dat4_file.open("BC4P1_S4.dat");
	pop2struc allfreq4;
	AlleleFrequencyDistribution(1,nl4x_4,L4num_fin4,allfreq4);
	for (int i=0;i<allfreq4.one.size();i++) {
		dat4_file << (allfreq4.one[i]+0.)/(2.0*nl4x_4) << " " << allfreq4.two[i] << "\n";
	}
	
	return 0;

	/*
	cout << "Female L4 worms going into first backcross\n";
	ConstructPop2FromEggs(1,0,nl4f[1],grid,cut,eggdist,L4num,logfactorial);

	
	cout << "Doing first backcross\n";
	Get2EggFrequencies (1,0,1,nl4f[1],s0,sh,L4num,eggdist);

	cout << "Collect L3 worms\n";
	ConstructPop2FromEggs(1,1,nl3[1],grid,cut,eggdist,L3num,logfactorial);

	cout << "Produce eggs in IVM sheep\n";
	Get2EggFrequencies (1,1,0,nl3[1],s0,sh,L3num,eggdist);

	cout << "Female L4 worms going into second backcross\n";
	ConstructPop2FromEggs(1,0,nl4f[2],grid,cut,eggdist,L4num,logfactorial);


	cout << "Doing second backcross\n";
	Get2EggFrequencies (1,0,1,nl4f[2],s0,sh,L4num,eggdist);
	
	cout << "Collect L3 worms\n";
	ConstructPop2FromEggs(1,1,nl3[2],grid,cut,eggdist,L3num,logfactorial);
	
	cout << "Produce eggs in IVM sheep\n";
	Get2EggFrequencies (1,1,0,nl3[2],s0,sh,L3num,eggdist);
	
	cout << "Female L4 worms going into third backcross\n";
	ConstructPop2FromEggs(1,0,nl4f[3],grid,cut,eggdist,L4num,logfactorial);

	
	cout << "Doing third backcross\n";
	Get2EggFrequencies (1,0,1,nl4f[3],s0,sh,L4num,eggdist);
	
	cout << "Collect L3 worms\n";
	ConstructPop2FromEggs(1,1,nl3[3],grid,cut,eggdist,L3num,logfactorial);
	
	cout << "Produce eggs in IVM sheep\n";
	Get2EggFrequencies (1,1,0,nl3[3],s0,sh,L3num,eggdist);
	
	cout << "Female L4 worms going into fourth backcross\n";
	ConstructPop2FromEggs(1,0,nl4f[4],grid,cut,eggdist,L4num,logfactorial);


	cout << "Doing fourth backcross\n";
	Get2EggFrequencies (1,0,1,nl4f[4],s0,sh,L4num,eggdist);
	
	cout << "Collect L3 worms\n";
	ConstructPop2FromEggs(1,1,nl3[4],grid,cut,eggdist,L3num,logfactorial);

	//Next, want to convert the L3 worms into an egg population via a sheep with or without the drug
	
	int nl3x=10000; //Assumed number of L3 worms for this part of the experiment*/
	
	/*cout << "Scale up the L3 population; there are no het resistant at this stage\n";
	pop3struc L3num2;
	for (int i=0;i<L3num.one.size();i++) {
		L3num2.one.push_back(0);
		L3num2.two.push_back(L3num.one[i]);
		L3num2.three.push_back(L3num.two[i]);
	}
	
	pop3struc BCrossOut=L3num2;
	
	cout << "Produce eggs in naive sheep\n";
	tvec eggdist2;
	
	Get3EggFrequencies (1,0,1,nl3[4],s0,sh,s1,BCrossOut,eggdist2);*/

	cout << "Collect L3 worms\n";
	ConstructPop3FromEggs(1,1,nl3x,grid,cut,eggdist2,L3num2,logfactorial);

	//pop3struc L3post=L3num2;
	
	cout << "Grow in naive sheep and collect worms\n";
	//Equivalent to calculating egg probabilities then bottlenecking to L4 population size.  The L3 chance of survival is equal to the L3 grows to L4 and makes egg rate
	Get3EggFrequencies (1,0,0,nl3x,s0,sh,s1,L3post,eggdist2);
	 nl4x=25;
	ConstructPop3FromEggs(1,0,nl4x,grid,cut,eggdist2,L3num2,logfactorial); //Need to write this function...

//	ofstream dat1_file;
	dat1_file.open("BC4P1.dat");
	for (int i=0;i<L3num2.one.size();i++) {
		dat1_file << L3num2.one[i] << " " << L3num2.two[i] << " " << L3num2.three[i] << "\n";
	}
	
	cout << "Grow in treated sheep and collect worms\n";
	//Equivalent to calculating egg probabilities then bottlenecking to L4 population size.  The L3 chance of survival is equal to the L3 grows to L4 and makes egg rate
	Get3EggFrequencies (1,1,0,nl3x,s0,sh,s1,L3post,eggdist2);
	nl4x=50;
	ConstructPop3FromEggs(1,0,nl4x,grid,cut,eggdist2,L3num2,logfactorial); //Need to write this function...
	
	//ofstream dat2_file;
	dat2_file.open("BC4P1_IVM.dat");
	for (int i=0;i<L3num2.one.size();i++) {
		dat2_file << L3num2.one[i] << " " << L3num2.two[i] << " " << L3num2.three[i] << "\n";
	}


	cout << "Produce eggs in treated sheep: round 1\n";
	Get3EggFrequencies (1,1,1,nl3[4],s0,sh,s1,BCrossOut,eggdist2);
	cout << "Collect L3 worms\n";
	ConstructPop3FromEggs(1,1,nl3x,grid,cut,eggdist2,L3num2,logfactorial);

	cout << "Produce eggs in treated sheep: round 2\n";
	Get3EggFrequencies (1,1,1,nl3x,s0,sh,s1,L3num2,eggdist2);
	cout << "Collect L3 worms\n";
	ConstructPop3FromEggs(1,1,nl3x,grid,cut,eggdist2,L3num2,logfactorial);

	L3post=L3num2;
	
	//Grow to L4 then sequence
	Get3EggFrequencies (1,1,0,nl3x,s0,sh,s1,L3post,eggdist2);
	nl4x=50;
	ConstructPop3FromEggs(1,0,nl4x,grid,cut,eggdist2,L3num2,logfactorial); //Need to write this function...

	//ofstream dat3_file;
	dat3_file.open("BC4P1_S3.dat");
	for (int i=0;i<L3num2.one.size();i++) {
		dat3_file << L3num2.one[i] << " " << L3num2.two[i] << " " << L3num2.three[i] << "\n";
	}

	cout << "Produce eggs in treated sheep: round 3\n";
	Get3EggFrequencies (1,1,1,nl3x,s0,sh,s1,L3post,eggdist2);
	cout << "Collect L3 worms\n";
	ConstructPop3FromEggs(1,1,nl3x,grid,cut,eggdist2,L3num2,logfactorial);

	//Grow to L4 then sequence
	Get3EggFrequencies (1,1,0,nl3x,s0,sh,s1,L3num2,eggdist2);
	nl4x=50;
	ConstructPop3FromEggs(1,0,nl4x,grid,cut,eggdist2,L3num2,logfactorial); //Need to write this function...
	
	//ofstream dat4_file;
	dat4_file.open("BC4P1_S4.dat");
	for (int i=0;i<L3num2.one.size();i++) {
		dat4_file << L3num2.one[i] << " " << L3num2.two[i] << " " << L3num2.three[i] << "\n";
	}
	
	cout << init_all << "\n";
	
	return 0;


}


