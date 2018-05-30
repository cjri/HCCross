#include "shared.h"
#include <iostream>
#include <string>

void find_logfact(vector<double>& logfactorial) {
	logfactorial.push_back(0);
	logfactorial.push_back(0);
	double lf=0;
	for (int i=2;i<=12000;i++) {
		lf=lf+log(i);
		logfactorial.push_back(lf);
	}
}

void SetupNumbersWRS (vector<int>& nl3, vector<int>& nl4m, vector<int>& nl4f) {
	nl4m.push_back(50);  //Initial cross
	nl4f.push_back(50);
	nl3.push_back(5000);
	
	nl4m.push_back(70);  //First backcross
	nl4f.push_back(90);
	nl3.push_back(8000);
	
	nl4m.push_back(45);  //Second backcross
	nl4f.push_back(55);
	nl3.push_back(10500);
	
	nl4m.push_back(100);  //Third backcross
	nl4f.push_back(73);
	nl3.push_back(10000);
	
	nl4m.push_back(100);  //Fourth backcross
	nl4f.push_back(100);
	nl3.push_back(10000);  //Assumed for transfer to either drug or naive sheep
}

void SetupNumbersCAVR (vector<int>& nl3, vector<int>& nl4m, vector<int>& nl4f) {
	nl4m.push_back(50);  //Initial cross
	nl4f.push_back(50);
	nl3.push_back(5000);
	
	nl4m.push_back(50);  //First backcross
	nl4f.push_back(66);
	nl3.push_back(8000);
	
	nl4m.push_back(57);  //Second backcross
	nl4f.push_back(45);
	nl3.push_back(10000);
	
	nl4m.push_back(100);  //Third backcross
	nl4f.push_back(100);
	nl3.push_back(10000);
	
	nl4m.push_back(86);  //Fourth backcross
	nl4f.push_back(90);
	nl3.push_back(10000);  //Assumed for transfer to either drug or naive sheeps
}

//New routines here
void ConstructPop2FromEggs (int verbose, int coarse, int N, int grid, double cut, dvec eggdist, pop2struc& L3num, vector<double> logfactorial) {
	L3num.one.clear();
	L3num.two.clear();
	
	pop2struc t;
	
	for (int i=0;i<=N;i++) {
		t.one.push_back(i);
		t.two.push_back(0);
	}
	
	//t will contain the full distribution of the number of worms that are heterozygote
	if (eggdist.one.size()==0) {
		L3num.two[0]=1;
	} else {
		for (unsigned int i=0;i<eggdist.one.size();i++) {
			int m=floor(eggdist.one[i]*N)+1;
			//cout << "M = " << m << "\n";
			for (int j=m;j<=N;j++) {
				double k=binomp(N,j,eggdist.one[i],logfactorial);
				if (k>cut) {
					t.two[j]=t.two[j]+(eggdist.two[i]*k);
				} else {
					break;
				}
			}
			for (int j=m-1;j>=0;j--) {
				double k=binomp(N,j,eggdist.one[i],logfactorial);
				if (k>cut) {
					t.two[j]=t.two[j]+(eggdist.two[i]*k);
				} else {
					break;
				}
			}
		}
	}
	
	
	double mean=0;
	if (coarse==1) {
		//Coarse grain full distribution onto a grid of specified number of points
		int min=-1;
		int max=0;
		for (int i=0;i<=N;i++) {
			if (t.two[i]>cut) {
				if (min==-1) {min=i;}
				max=i;
				mean=mean+(t.two[i]*i);
			}
		}
		cout << "Min " << min << " Max " << max << "\n";
		int cse=(max-min)/grid;
		if (cse%2==0) {
			cse++;
		}
		for (int i=min;i<=max;i=i+cse) {
			double id=0;
			for (int j=i;j<i+cse;j++) {
				id=id+t.two[j];
			}
			L3num.one.push_back(i+(cse/2));
			L3num.two.push_back(id);
		}
	} else {
		for (int i=0;i<t.one.size();i++) {
			if (t.two[i]>0) {
				L3num.one.push_back(t.one[i]);
				L3num.two.push_back(t.two[i]);
			}
		}
	}
	if (verbose==1) {
		//cout << "Heterozygous worms\n";
		double sum=0;
		mean=0;
		for (int i=0;i<L3num.one.size();i++) {
			if (L3num.two[i]>0) {
				cout << L3num.one[i] << " " << L3num.two[i] << "\n";
				sum=sum+L3num.two[i];
				mean=mean+(L3num.one[i]*L3num.two[i]);
			}
		}
		cout << "Sum of values is " << sum << "\n";
		//cout << "Mean is " << mean << " sum is " << sum << "\n";
	}
}

void GetL4fromL3 (int verbose, int sel, int NL3, int NL4, double sh, double s0, pop2struc L3num, pop2struc& L4num, vector<double> logfactorial) {
	L4num.one.clear();
	L4num.two.clear();
	pop2struc t;
	for (int i=0;i<=NL4;i++) {
		t.one.push_back(i);
		t.two.push_back(0);
	}
	for (int i=0;i<L3num.one.size();i++) {
		double het_freq=L3num.one[i]/(NL3+0.);
		double scaled_freq=het_freq;
		if (sel==1) {
			double w_tot=(het_freq*sh)+((1-het_freq)*s0);
			scaled_freq=(het_freq*sh)/w_tot;
		}
		double p=L3num.two[i];
		for (int j=0;j<t.one.size();j++) {
			t.two[j]=t.two[j]+(p*binomp(NL4,j,scaled_freq,logfactorial));
		}
	}
	if (verbose==1) {
		double tot=0;
		for (int i=0;i<t.one.size();i++) {
			cout << t.one[i] << " " << t.two[i] << "\n";
			tot=tot+t.two[i];
		}
		cout << "Sum of values is " << tot << "\n";
	}
	L4num=t;

}

void Get2EggFrequencies (int verbose, int sel, int backcross, int N, double s0, double sh, pop2struc L4num, dvec& eggdist) {
	//cout << "New egg distribution\n";
	eggdist.one.clear();
	eggdist.two.clear();
	for (int i=0;i<L4num.one.size();i++) {
		double wh=(L4num.one[i]+0.)/(N+0.);  //Note that N is the maximum possible size recorded in L4num, as implemented
		double w0=(N-L4num.one[i]+0.)/(N+0.);
		if (backcross==1) {
			double ph=wh*0.5;
			double p0=(wh*0.5)+w0;
			wh=ph;
			w0=p0;
		}
		if (sel==1) {
			double sumfit=(sh*wh)+(s0*w0);
			wh=(wh*sh)/sumfit;
			w0=(w0*s0)/sumfit;
		}
		eggdist.one.push_back(wh);
		eggdist.two.push_back(L4num.two[i]);
	}
	if (verbose==1) {
		double tot=0;
		for (int i=0;i<eggdist.one.size();i++) {
			cout << eggdist.one[i] << " " << eggdist.two[i] << "\n";
			tot=tot+eggdist.two[i];
		}
		cout << "Total is " << tot << "\n";
	}
}

void Get3EggFrequencies (int verbose, int sel, int cross, int N, double s0, double sh, double s1, pop3struc L3num, tvec& eggdist) {
	//cout << "New egg distribution\n";
	eggdist.one.clear();
	eggdist.two.clear();
	eggdist.three.clear();
	for (int i=0;i<L3num.one.size();i++) {
		double w1=(L3num.one[i]+0.)/(N+0);
		double wh=(L3num.two[i]+0.)/(N+0.);
		double w0=(N-L3num.one[i]-L3num.two[i]+0.)/(N+0.);
		//Crossing
		if (cross==1) {
			double p1=(w1*w1)+(w1*wh)+(0.25*wh*wh);
			double ph=(w1*wh)+(2*w1*w0)+(0.5*wh*wh)+(wh*w0);
			double p0=(0.25*wh*wh)+(wh*w0)+(w0*w0);
			w1=p1;
			wh=ph;
			w0=p0;
		}
		if (sel==1) {
			//			cout << w1 << " " << wh << " " << w0 << " ";
			double sumfit=(s1*w1)+(sh*wh)+(s0*w0);
			w1=(w1*s1)/sumfit;
			wh=(wh*sh)/sumfit;
			w0=(w0*s0)/sumfit;
			//			cout << w1 << " " << wh << " " << w0 << "\n";
		}
		eggdist.one.push_back(w1);
		eggdist.two.push_back(wh);
		eggdist.three.push_back(L3num.three[i]);
	}
	if (verbose==1) {
		double tot=0;
		for (int i=0;i<eggdist.one.size();i++) {
			cout << eggdist.one[i] << " " << eggdist.two[i] << " " << eggdist.three[i] << "\n";
			tot=tot+eggdist.three[i];
		}
		cout << "Total is " << tot << "\n";
	}
}

void ConstructPop3FromEggs (int verbose, int coarse, int N, int grid, double cut, tvec eggdist, pop3struc& L3num, vector<double> logfactorial) {
	L3num.one.clear();
	L3num.two.clear();
	L3num.three.clear();
	double eps=1e-10;
	for (int i=0;i<eggdist.one.size();i++) {  //This does the full calculation...
		if (eggdist.three[i]>cut) {
			//		cout << "i= "<< i << "\n";
			int m1=floor(eggdist.one[i]*N)+1; //Mean hom res
			int m2=floor(eggdist.two[i]*N)+1; //Mean het
			if (m1>N) {m1=N;}
			if (m2>N) {m2=N;}
			
			vector<double> prb;
			vector<int> counts;
			prb.clear();
			double p1=eggdist.one[i];
			double p2=eggdist.two[i];
			double q=1-eggdist.one[i]-eggdist.two[i];
			if (q<1e-10) {
				q=0;
				p2=1-p1;
			}
			prb.push_back(p1);
			prb.push_back(p2);
			prb.push_back(q);
			//cout << "Probs " << prb[0] << " " << prb[1] << " " << prb[2] << " " << prb[0]+prb[1]+prb[2] << "\n";
			
			//If q is zero, then have a binomial case, and need to avoid breaking the loop too soon
			if (q==0) {
				//cout << "Doing binomial...\n";
				for (int j1=0;j1<=N;j1++) {  //Calculate binomial probs
					double k=binomp(N,j1,eggdist.one[i],logfactorial);
					if (k>cut) {
						//cout << "j1 " << j1 << " j2 " << N-j1 << " k " << k << "\n";
						L3num.one.push_back(j1);
						L3num.two.push_back(N-j1);
						L3num.three.push_back(eggdist.three[i]*k);
					}
				}
			} else {
				for (int j1=m1;j1<=N;j1++) {  //Calculate multinomial probs
					//cout << j1 << "\n";
					for (int j2=m2;j2<=N-j1;j2++) {
						counts.clear();
						counts.push_back(j1);
						counts.push_back(j2);
						counts.push_back(N-j1-j2);
						double k=eggdist.three[i]*multinomp(N,counts,prb,logfactorial);
						//if (k>1) {
						//cout << "Error N " << N << " Counts: " << j1 << " " << j2 << " " << N-j1-j2 << " Prb " << prb[0] << " " << prb[1] << " " << prb[2] << " " << k << "\n";
						//cout << "Error1 N " << N << " Counts: " << j1 << " " << j2 << " " << N-j1-j2 << " " << k << "\n";
						//}
						if (k>cut) {
							//cout << "j1 " << j1 << " j2 " << j2 << " k " << k << "\n";
							L3num.one.push_back(j1);
							L3num.two.push_back(j2);
							L3num.three.push_back(k);
						} else {
							break;
						}
					}
					for (int j2=m2-1;j2>=0;j2--) {
						if (j2<=N-j1) {
							counts.clear();
							counts.push_back(j1);
							counts.push_back(j2);
							counts.push_back(N-j1-j2);
							double k=eggdist.three[i]*multinomp(N,counts,prb,logfactorial);
							//if (k>1) {
							//cout << "Error N " << N << " Counts: " << j1 << " " << j2 << " " << N-j1-j2 << " Prb " << prb[0] << " " << prb[1] << " " << prb[2] << " " << k << "\n";
							//	cout << "Error2 N " << N << " Counts: " << j1 << " " << j2 << " " << N-j1-j2 << " " << k << "\n";
							//}
							
							if (k>cut) {
								L3num.one.push_back(j1);
								L3num.two.push_back(j2);
								L3num.three.push_back(k);
							} else {
								break;
							}
						}
					}
				}
				for (int j1=m1-1;j1>=0;j1--) {  //Calculate multinomial probs
					for (int j2=m2;j2<=N-j1;j2++) {
						counts.clear();
						counts.push_back(j1);
						counts.push_back(j2);
						counts.push_back(N-j1-j2);
						double k=eggdist.three[i]*multinomp(N,counts,prb,logfactorial);
						//if (k>1) {
						//cout << "Error N " << N << " Counts: " << j1 << " " << j2 << " " << N-j1-j2 << " Prb " << prb[0] << " " << prb[1] << " " << prb[2] << " " << k << "\n";
						//	cout << "Error3 N " << N << " Counts: " << j1 << " " << j2 << " " << N-j1-j2 << " " << k << "\n";
						//}
						
						if (k>cut) {
							L3num.one.push_back(j1);
							L3num.two.push_back(j2);
							L3num.three.push_back(k);
						} else {
							break;
						}
					}
					for (int j2=m2-1;j2>=0;j2--) {
						if (j2<=N-j1) {
							counts.clear();
							counts.push_back(j1);
							counts.push_back(j2);
							counts.push_back(N-j1-j2);
							double k=eggdist.three[i]*multinomp(N,counts,prb,logfactorial);
							//if (k>1) {
							//cout << "Error N " << N << " Counts: " << j1 << " " << j2 << " " << N-j1-j2 << " Prb " << prb[0] << " " << prb[1] << " " << prb[2] << " " << k << "\n";
							//	cout << "Error4 N " << N << " Counts: " << j1 << " " << j2 << " " << N-j1-j2 << " " << k << "\n";
							//}
							
							if (k>cut) {
								L3num.one.push_back(j1);
								L3num.two.push_back(j2);
								L3num.three.push_back(k);
							} else {
								break;
							}
						}
					}
				}
			}
		}
	}
	
	double tot=0;
	for (int i=0;i<L3num.one.size();i++) {
		//cout << L3num.one[i] << " " << L3num.two[i] << " " << L3num.three[i] << "\n";
		tot=tot+L3num.three[i];
	}
	cout << "Preliminary Total " << tot << "\n";
	
	
	if (coarse==1) {
		CoarseGrainL3(grid,L3num);
	} else {
		CombineL3(L3num);
	}
	if (verbose==1) {
		double tot=0;
		for (int i=0;i<L3num.one.size();i++) {
			cout << L3num.one[i] << " " << L3num.two[i] << " " << L3num.three[i] << "\n";
			tot=tot+L3num.three[i];
		}
		cout << "Total " << tot << "\n";
	}
	
}

void CoarseGrainL3 (int grid, pop3struc& L3num) {
	cout << "Coarse-grain\n";
	pop3struc temp;
	int min1=1000;
	int min2=1000;
	int max1=0;
	int max2=0;
	for (int i=0;i<L3num.one.size();i++) {
		if (L3num.one[i]>max1) {
			max1=L3num.one[i];
		}
		if (L3num.one[i]<min1) {
			min1=L3num.one[i];
		}
		if (L3num.two[i]>max2) {
			max2=L3num.two[i];
		}
		if (L3num.two[i]<min2) {
			min2=L3num.two[i];
		}
	}
	cout << "Minmax values " << min1 << " " << max1 << " | " << min2 << " " << max2 << "\n";
	
	int step1=(max1-min1)/grid;
	int step2=(max2-min2)/grid;
	
	temp.one.clear();
	temp.two.clear();
	temp.three.clear();
	double tot=0;
	double x=0;
	for (int i=min1;i<=max1+step1;i=i+step1) {
		//cout << "Here i " << i << "\n";
		for (int j=min2;j<=max2+step2;j=j+step2) {
			x=0;
			for (int k=0;k<L3num.one.size();k++) {
				if (L3num.one[k]>=i&&L3num.one[k]<(i+step1)&&L3num.two[k]>=j&&L3num.two[k]<(j+step2)) {
					x=x+L3num.three[k];
					tot=tot+L3num.three[k];
				}
			}
			if (x>0) {
				temp.one.push_back(i+(step1/2));
				temp.two.push_back(j+(step2/2));
				temp.three.push_back(x);
			}
		}
	}
	cout << "Total is " << tot << "\n";
	L3num=temp;
}

void CombineL3 (pop3struc& L3num) {
	cout << "Combine\n";
	pop3struc temp;
	int min1=1000;
	int min2=1000;
	int max1=-1;
	int max2=-1;
	for (int i=0;i<L3num.one.size();i++) {
		if (L3num.one[i]>max1) {
			max1=L3num.one[i];
		}
		if (L3num.one[i]<min1) {
			min1=L3num.one[i];
		}
		if (L3num.two[i]>max2) {
			max2=L3num.two[i];
		}
		if (L3num.two[i]<min2) {
			min2=L3num.two[i];
		}
	}
	cout << "Minmax values " << min1 << " " << max1 << " | " << min2 << " " << max2 << "\n";
	
	temp.one.clear();
	temp.two.clear();
	temp.three.clear();
	double tot=0;
	double x=0;
	for (int i=min1;i<=max1;i++) {
		//	cout << "Here i " << i << "\n";
		for (int j=min2;j<=max2;j++) {
			//		cout << "Here j " << j << " " << min2 << " " << max1 << "\n";
			x=0;
			for (int k=0;k<L3num.one.size();k++) {
				if (L3num.one[k]==i&&L3num.two[k]==j) {
					x=x+L3num.three[k];
					tot=tot+L3num.three[k];
				}
			}
			if (x>0) {
				temp.one.push_back(i);
				temp.two.push_back(j);
				temp.three.push_back(x);
			}
		}
	}
	cout << "Total is " << tot << "\n";
	L3num=temp;
}

void GetL4fromL3_2 (int verbose, int sel, int NL3, int NL4, double s1, double sh, double s0, pop3struc L3num, pop3struc& L4num, vector<double> logfactorial) {
	L4num.one.clear();
	L4num.two.clear();
	L4num.three.clear();
	pop3struc t;
	for (int i=0;i<=NL4;i++) {
		for (int j=0;j<=NL4;j++) {
			t.one.push_back(i);
			t.two.push_back(j);
			t.three.push_back(0);
		}
	}
	vector<int> counts;
	vector<double> prb;
	counts.push_back(1);
	prb.push_back(0);
	for (int i=0;i<L3num.one.size();i++) {
		double hom_freq=L3num.one[i]/(NL3+0.);
		double het_freq=L3num.two[i]/(NL3+0.);
		double sus_freq=1-hom_freq-het_freq;
		if (sel==1) {
			double w_tot=(hom_freq*s1)+(het_freq*sh)+(sus_freq*s0);
			hom_freq=(hom_freq*s1)/w_tot;
			het_freq=(het_freq*sh)/w_tot;
			sus_freq=(sus_freq*s0)/w_tot;
		}
		double p=L3num.three[i];
		for (int j=0;j<t.one.size();j++) {
			counts.clear();
			counts.push_back(t.one[j]);
			counts.push_back(t.two[j]);
			counts.push_back(NL4-t.one[j]-t.two[j]);
			prb.clear();
			prb.push_back(hom_freq);
			prb.push_back(het_freq);
			prb.push_back(sus_freq);
			t.three[j]=t.three[j]+(p*multinomp(NL4,counts,prb,logfactorial));
		}
	}
	
	if (verbose==1) {
		double tot=0;
		for (int i=0;i<t.one.size();i++) {
			cout << t.one[i] << " " << t.two[i] << " " << t.three[i] << "\n";
			tot=tot+t.three[i];
		}
		cout << "Sum of values is " << tot << "\n";
	}
	L4num=t;
	
}

void AlleleFrequencyDistribution (int verbose, int N4, pop3struc L4num, pop2struc& allfreq) {
	allfreq.one.clear();
	allfreq.two.clear();
	for (int i=0;i<=N4*2;i++) {
		allfreq.one.push_back(i);
		allfreq.two.push_back(0);
	}
	for (int i=0;i<L4num.one.size();i++) {
		int count = (L4num.one[i]*2)+L4num.two[i];
		allfreq.two[count]=allfreq.two[count]+L4num.three[i];
	}
	if (verbose==1) {
		double tot=0;
		for (int i=0;i<allfreq.one.size();i++) {
			cout << (allfreq.one[i]+0.)/(2.0*N4) << " " << allfreq.two[i] << "\n";
			tot=tot+allfreq.two[i];
		}
		cout << "Sum of values is " << tot << "\n";
	}
}

double binomp(int n, int count, double prb, vector<double> logfactorial) {
	double mp=logfactorial[n];
	if (prb==0) {
		if (count==0) {
			mp=1;
		} else {
			mp=0;
		}
	} else if (prb==1) {
		if (count==n) {
			mp=1;
		} else {
			mp=0;
		}
	} else if ((count>n)||(count<0)) {
		mp=0;
	} else {
		int nc=n-count;
		mp=mp-logfactorial[count]-logfactorial[nc];
		mp=mp+(count*log(prb));
		mp=mp+(nc*log(1-prb));
		mp=exp(mp);
	}
	return(mp);
}

double multinomp(int n, vector<int> counts, vector<double> prb, vector<double> logfactorial) {
	double mp=logfactorial[n];
	int iszero=0;
	for (unsigned int i=0;i<prb.size();i++) {
		if (prb[i]==0&&counts[i]>0) {
			iszero=1;
		}
	}
	if (iszero==1) {
		mp=0;
	} else {
		int reduce=0;
		for (unsigned int i=0;i<prb.size();i++) {
			if (prb[i]==0) {
				reduce=1;
			}
		}
		
		if (reduce==1) {
			vector<int> counts2;
			vector<int> prb2;
			for (unsigned int i=0;i<prb.size();i++) {
				if (prb[i]!=0) {
					counts2.push_back(counts[i]);
					prb2.push_back(prb[i]);
				}
			}
			for (unsigned int i=0;i<counts2.size();i++) {
				mp=mp-logfactorial[counts2[i]];
				mp=mp+(counts2[i]*log(prb2[i]));
			}
			mp=exp(mp);
			
		} else {
			for (unsigned int i=0;i<counts.size();i++) {
				mp=mp-logfactorial[counts[i]];
				mp=mp+(counts[i]*log(prb[i]));
			}
			mp=exp(mp);
		}
	}
	return(mp);
}





//Cut this routine...
void GetEggDist (int verbose, double my, double mf, double m1, double cut, dvec& eggdist, int grid) {
	//cout << "Get Egg dist \n";
	
	//cout << "Parameters " << my << " " << mf << " " << m1 << "\n";
	
    double sy=sqrt(my);
	double s1=sqrt(m1);
    //Calibrate scale of calculation
    double scale=0.25;
    double next=func_p(mf+(grid*scale),m1,s1,my,sy);
    for (int i=1;i<=1000;i++) {
        if (next<cut) {
            scale=scale/1.1;
            next=func_p(mf+(grid*scale),m1,s1,my,sy);
        } else {
            break;
        }
    }
    dvec ptab;
    double x=mf-(grid*scale);
	
	if (verbose==1) {
		cout << x << " " << grid << " " << scale << "\n";
	}
	
	//Check for probability at frequencies less than zero - exclude these...
	if (x<0) {
		x=0;
		scale=(mf+(grid*scale))/(2*(grid+1));
	}
	while (x<=mf+((grid+0.01)*scale)) {
		ptab.one.push_back(func_r(x));
		ptab.two.push_back(scale*func_p(x,m1,s1,my,sy));
		x=x+scale;
	}
	
	double tot=0;
	for (unsigned int i=0;i<ptab.one.size();i++) {
		tot=tot+ptab.two[i];
	}
	if (verbose==1) {
		cout << "Total here is " << tot << "\n";
	}

	eggdist.one=ptab.one;
	eggdist.two=ptab.two;
	
	if (verbose==1) {
		double mean=0;
		for (unsigned int i=0;i<eggdist.one.size();i++) {
			cout << eggdist.one[i] << " " << eggdist.two[i] << "\n";
			mean=mean+(eggdist.one[i]*eggdist.two[i]);
		}
		cout << "Mean is " << mean << "\n";
	}
}

double func_p (double x, double m1, double s1, double m2, double s2) {
    double pi=3.1415926535897932385;
    double a=func_a(x,s1,s2);
    double b=func_b(x,m1,s1,m2,s2);
    double c=func_c(m1,s1,m2,s2);
    double d=(pow(b,2)-(c*pow(a,2)))/(2*pow(a,2));
    d=exp(d);
    double p1=((b*d)/(sqrt(2*pi)*s1*s2*(pow(a,3))));
    double p2=emod(b/a)-emod(-b/a);
    double p3=1/(pi*s1*s2*pow(a,2));
    double p4=exp(-c/2);
    double p=(p1*p2)+(p3*p4);
    return p;
}

double func_a (double x, double s1, double s2) {
    double a=sqrt((pow(x,2)/pow(s1,2))+1/pow(s2,2));
    return a;
}

double func_b (double x, double m1, double s1, double m2, double s2) {
    double b=(m1*x)/pow(s1,2)+(m2/pow(s2,2));
    return b;
}

double func_c (double m1, double s1, double m2, double s2) {
    double c=(pow(m1,2)/pow(s1,2))+(pow(m2,2)/pow(s2,2));
    return c;
}

double func_r (double x) {
	double rx=x/(x+1);
	return rx;
}






















void CombineL3Struc (L3struc& L3s) {
	L3struc temp;
	int max1=0;
	int max2=0;
	for (int i=0;i<L3s.one.size();i++) {
		if (L3s.one[i]>max1) {
			max1=L3s.one[i];
		}
		if (L3s.two[i]>max2) {
			max2=L3s.two[i];
		}
	}
	cout << "Max values " << max1 << " " << max2 << "\n";
	
	//	int nm1=max1/grid;
	//	int nm2=max2/grid;
	
	temp.one.clear();
	temp.two.clear();
	temp.three.clear();
	double tot=0;
	double x=0;
	for (int i=0;i<=max1;i++) {
		//	cout << "Here i " << i << "\n";
		for (int j=0;j<=max2;j++) {
			x=0;
			for (int k=0;k<L3s.one.size();k++) {
				if (L3s.one[k]==i&&L3s.two[k]==j) {
					x=x+L3s.three[k];
					tot=tot+L3s.three[k];
				}
			}
			if (x>0) {
				temp.one.push_back(i);
				temp.two.push_back(j);
				temp.three.push_back(x);
			}
		}
	}
	cout << "Total is " << tot << "\n";
	L3s=temp;
}






void L3fromEggs (int verbose, int nl3, int grid, dvec eggdist, map<int,double>& L3num, double cut, vector<double> logfactorial) {
	map<int,double> t;
	for (int i=0;i<=nl3;i++) {
		t[i]=0;
	}
	
	if (eggdist.one.size()==0) {
		L3num[0]=1;
	} else {
	
		if (verbose==1) {
			for (int i=0;i<eggdist.one.size();i++) {
				cout << "Here " << eggdist.one[i] << " " << eggdist.two[i] << "\n";
			}
		}

		for (unsigned int i=0;i<eggdist.one.size();i++) {
			int m=floor(eggdist.one[i]*nl3)+1;
			for (int j=m;j<=nl3;j++) {
				double k=binomp(nl3,j,eggdist.one[i],logfactorial);
				//	if (verbose==1) {
				//		cout << nl3 << " " << j << " " << eggdist.one[i] << " " << k << "\n";
				//	}
				if (k>cut) {
					t[j]=t[j]+(eggdist.two[i]*k);
				} else {
					break;
				}
			}
			for (int j=m-1;j>0;j--) {
				double k=binomp(nl3,j,eggdist.one[i],logfactorial);
				if (k>cut) {
					t[j]=t[j]+(eggdist.two[i]*k);
				} else {
					break;
				}
			}
		}
	
		//Coarse grain
		int min=0;
		int max=0;
		double mean=0;
		for (int i=0;i<=nl3;i++) {
			if (t[i]>0) {
//				if (verbose==1) {
//					cout << i << " " << t[i] << "\n";
//				}
				if (min==0) {min=i;}
				max=i;
				mean=mean+(t[i]*i);
			}
		}
	
		if (verbose==1) {
			cout << "Mean is " << mean << "\n";
		}
		//cout << min << " " << max << "\n";
		int cse=(max-min)/grid;
		//cout << cse << "\n";
		if (cse%2==0) {
			cse++;
		}
		int mean_int=mean;
		for (int i=min;i<=max;i=i+cse) {
			double id=0;
			for (int j=i;j<i+cse;j++) {
				id=id+t[j];
			}
			L3num[i+(cse/2)]=id;
		}
	

		if (verbose==1) {
			cout << "L3 resistant numbers\n";
			double sum=0;
			mean=0;
			for (int i=0;i<L3num.size();i++) {
				if (L3num[i]>0) {
					cout << i << " " << L3num[i] << "\n";
					sum=sum+L3num[i];
					mean=mean+(i*L3num[i]);
				}
			}
//			cout << sum << "\n";
			cout << "Mean is " << mean << " sum is " << sum << "\n";
		}
	}
}



void SheepSurvivalHetSus (int verbose, double sh, double s0, double cut, int nl3, map<int,double> L3num, dvec& resfreq, vector<double> logfactorial) {
	if (verbose==1) {
		cout << "SheepSurvival\n";
	}
	tvec temp;
	for (int i=0;i<=nl3;i++) {
		if (L3num[i]>0) {
			temp.one.push_back(i);
			temp.two.push_back(nl3-i);
			temp.three.push_back(L3num[i]);
		}
	}
	
	if (verbose==1) {
		cout << "Temp size " << temp.one.size() << "\n";
		for (unsigned int i=0;i<temp.one.size();i++) {
			cout << "Temp i " << i << " " << temp.one[i] << " " << temp.two[i] << " " << temp.three[i] << "\n";
		}
	}
	
	if (temp.one.size()==0) {
		resfreq.one.push_back(0);
		resfreq.two.push_back(1);
	} else {
	
		dvec het; //Number of heterozygous worms remaining
		dvec sus; //Number of homozygous susceptible worms remaining
		tvec temp2;
		for (unsigned int i=0;i<temp.one.size();i++) {
			//Number of hetero worms
			//	cout << "i= " << i << " " << temp.one[i] << " " << temp.two[i] << " " << temp.three[i] << "\n";
			//	cout << "Het\n";
			het.one.clear();
			het.two.clear();
			int im=floor(temp.one[i]*sh)+1;
			for (int j=im;j<=nl3;j++) {
				double k=binomp(temp.one[i],j,sh,logfactorial);
				if (k>cut) {
					het.one.push_back(j);
					het.two.push_back(k);
				} else {
					break;
				}
			}
		
			for (int j=im-1;j>=0;j--) {
				double k=binomp(temp.one[i],j,sh,logfactorial);
				if (k>cut) {
					het.one.push_back(j);
					het.two.push_back(k);
				} else {
					break;
				}
			}
		
	
			if (verbose==1) {
				cout << "Het\n";
				for (int j=0;j<het.one.size();j++) {
					cout << j << " " << het.one[j] << " " << het.two[j] << "\n";
				}
			}
			
			//Number of susceptible worms
			//	cout <<  "Sus \n";
			sus.one.clear();
			sus.two.clear();
			im=floor(temp.two[i]*s0)+1;
			for (int j=im;j<=nl3;j++) {
				double k=binomp(temp.two[i],j,s0,logfactorial);
				//cout << "Calc " << im-1 << " " << j << " " << s0 << " " << k << "\n";
				if (k>cut) {
					sus.one.push_back(j);
					sus.two.push_back(k);
				} else {
					break;
				}
			}
		
			for (int j=im-1;j>=0;j--) {
				double k=binomp(temp.two[i],j,s0,logfactorial);
				//cout << "Calc " << im-1 << " " << j << " " << s0 << " " << k << "\n";
				if (k>cut) {
					sus.one.push_back(j);
					sus.two.push_back(k);
				} else {
					break;
				}
			}
		
			if (verbose==1) {
				cout << "Sus\n";
				for (int j=0;j<sus.one.size();j++) {
					cout << j << " " << sus.one[j] << " " << sus.two[j] << "\n";
				}
			}
			
			for (int j=0;j<het.one.size();j++) {
				for (int k=0;k<sus.one.size();k++) {
					temp2.one.push_back(het.one[j]);
					temp2.two.push_back(sus.one[k]);
					temp2.three.push_back(het.two[j]*sus.two[k]*temp.three[i]);
				}
			}
		}
	
		for (int i=0;i<temp2.one.size();i++) {
			if (temp2.one[i]==0) {
				resfreq.one.push_back(0);
			} else {
				resfreq.one.push_back(temp2.one[i]/(temp2.one[i]+temp2.two[i]));
			}
			resfreq.two.push_back(temp2.three[i]);
		}
	
		if (verbose==1) {
			double sum=0;
			double mean=0;
			for (int i=0;i<resfreq.one.size();i++) {
				cout << resfreq.one[i] << " " << resfreq.two[i] << "\n";
				sum=sum+resfreq.two[i];
				mean=mean+resfreq.one[i]*resfreq.two[i];
			}
			cout << "Total " << sum << " Mean " << mean << "\n";
		}
	}
}

void L4fromL3 (int verbose, int nl4f, double wp, dvec resfreq, map<int,double>& nextgen, vector<double> logfactorial) {
	double tt=0;
	//if (verbose==1) {
	//	cout << "resfreq size " << resfreq.one.size() << "\n";
	//	for (int i=0;i<resfreq.one.size();i++) {
	//		cout << "RF " << i << " " << resfreq.one[i] << " " << resfreq.two[i] << "\n";
	//	}
	//}
	if (resfreq.one.size()==0) {
		nextgen[0]=nextgen[0]+wp;
	} else {
		for (int j=0;j<resfreq.one.size();j++) {
			//cout << j << " " << resfreq.one.size() << " " << nl4f << " " << resfreq.one[j] << " " << resfreq.two[j] << "\n";
			for (int i=0;i<=nl4f;i++) {
				nextgen[i]=nextgen[i]+(wp*(resfreq.two[j]*binomp(nl4f,i,resfreq.one[j],logfactorial)));
				if (verbose==1) {
				//	cout << i << " " << (wp*(resfreq.two[j]*binomp(nl4f,i,resfreq.one[j],logfactorial))) <<" " <<  nextgen[i] << "\n";
					//			tt=tt+(wp*(resfreq.two[j]*binomp(nl4f,i,resfreq.one[j],logfactorial)));
		//			cout << "Function " << i << " of " << nl4f << " " << resfreq.one[j] << " " << i << " " << tt << "\n";
				}
			}
		}
	}
//	cout << "TT " << tt << " " << wp << "\n";
	if (verbose==1) {
		double tot=0;
		for (int i=0;i<=nl4f;i++) {
			cout << "Pop " << i << " " << nextgen[i] << "\n";
			tot=tot+nextgen[i];
		}
		cout << "Check total probability " << tot << "\n";
	}
}



void L3fromEggs2 (int verbose, int nl3, int grid, tvec eggdist, map<int,double>& L3num, double cut, vector<double> logfactorial) {
	map<int,double> t;
	for (int i=0;i<=nl3;i++) {
		t[i]=0;
	}
	
	if (eggdist.one.size()==0) {
		L3num[0]=1;
	} else {
		
		if (verbose==1) {
			for (int i=0;i<eggdist.one.size();i++) {
				cout << "Here " << eggdist.one[i] << " " << eggdist.two[i] << "\n";
			}
		}
		
		for (unsigned int i=0;i<eggdist.one.size();i++) {
			int m=floor(eggdist.one[i]*nl3)+1;
			for (int j=m;j<=nl3;j++) {
				double k=binomp(nl3,j,eggdist.one[i],logfactorial);
				//	if (verbose==1) {
				//		cout << nl3 << " " << j << " " << eggdist.one[i] << " " << k << "\n";
				//	}
				if (k>cut) {
					t[j]=t[j]+(eggdist.two[i]*k);
				} else {
					break;
				}
			}
			for (int j=m-1;j>0;j--) {
				double k=binomp(nl3,j,eggdist.one[i],logfactorial);
				if (k>cut) {
					t[j]=t[j]+(eggdist.two[i]*k);
				} else {
					break;
				}
			}
		}
		
		//Coarse grain
		int min=0;
		int max=0;
		double mean=0;
		for (int i=0;i<=nl3;i++) {
			if (t[i]>0) {
				//				if (verbose==1) {
				//					cout << i << " " << t[i] << "\n";
				//				}
				if (min==0) {min=i;}
				max=i;
				mean=mean+(t[i]*i);
			}
		}
		
		if (verbose==1) {
			cout << "Mean is " << mean << "\n";
		}
		//cout << min << " " << max << "\n";
		int cse=(max-min)/grid;
		//cout << cse << "\n";
		if (cse%2==0) {
			cse++;
		}
		int mean_int=mean;
		for (int i=min;i<=max;i=i+cse) {
			double id=0;
			for (int j=i;j<i+cse;j++) {
				id=id+t[j];
			}
			L3num[i+(cse/2)]=id;
		}
		
		
		if (verbose==1) {
			cout << "L3 resistant numbers\n";
			double sum=0;
			mean=0;
			for (int i=0;i<L3num.size();i++) {
				if (L3num[i]>0) {
					cout << i << " " << L3num[i] << "\n";
					sum=sum+L3num[i];
					mean=mean+(i*L3num[i]);
				}
			}
			//			cout << sum << "\n";
			cout << "Mean is " << mean << " sum is " << sum << "\n";
		}
	}
}

void GetL3Struc (int n, double cut, L3struc& L3s, tvec eggdist, vector<double> logfactorial) {
	for (int i=0;i<eggdist.one.size();i++) {  //This does the full calculation...
		if (eggdist.three[i]>cut) {
	//		cout << "i= "<< i << "\n";
			int m1=floor(eggdist.one[i]*n)+1; //Mean hom res
			int m2=floor(eggdist.two[i]*n)+1; //Mean het
			vector<double> prb;
			vector<int> counts;
			prb.clear();
			prb.push_back(eggdist.one[i]);
			prb.push_back(eggdist.two[i]);
			prb.push_back(1-eggdist.one[i]-eggdist.two[i]);
			for (int j1=m1;j1<=n;j1++) {  //Calculate multinomial probs
				//cout << j1 << "\n";
				for (int j2=m2;j2<=n-j1;j2++) {
					counts.clear();
					counts.push_back(j1);
					counts.push_back(j2);
					counts.push_back(n-j1-j2);
					double k=eggdist.three[i]*multinomp(n,counts,prb,logfactorial);
					if (k>cut) {
						L3s.one.push_back(j1);
						L3s.two.push_back(j2);
						L3s.three.push_back(k);
					} else {
						break;
					}
				}
				for (int j2=m2-1;j2>=0;j2--) {
					counts.clear();
					counts.push_back(j1);
					counts.push_back(j2);
					counts.push_back(n-j1-j2);
					double k=eggdist.three[i]*multinomp(n,counts,prb,logfactorial);
					if (k>cut) {
						L3s.one.push_back(j1);
						L3s.two.push_back(j2);
						L3s.three.push_back(k);
					} else {
						break;
					}
				}
			}
			for (int j1=m1-1;j1>=0;j1--) {  //Calculate multinomial probs
				for (int j2=m2;j2<=n-j1;j2++) {
					counts.clear();
					counts.push_back(j1);
					counts.push_back(j2);
					counts.push_back(n-j1-j2);
					double k=eggdist.three[i]*multinomp(n,counts,prb,logfactorial);
					if (k>cut) {
						L3s.one.push_back(j1);
						L3s.two.push_back(j2);
						L3s.three.push_back(k);
					} else {
						break;
					}
				}
				for (int j2=m2-1;j2>=0;j2--) {
					counts.clear();
					counts.push_back(j1);
					counts.push_back(j2);
					counts.push_back(n-j1-j2);
					double k=eggdist.three[i]*multinomp(n,counts,prb,logfactorial);
					if (k>cut) {
						L3s.one.push_back(j1);
						L3s.two.push_back(j2);
						L3s.three.push_back(k);
					} else {
						break;
					}
				}
			}
		}
	}
}

void CalculatEggDist (int sel, int nl3x, double cut, double s1, double sh, double s0, L3struc L3s, tvec& eggdist) {
	eggdist.one.clear();
	eggdist.two.clear();
	eggdist.three.clear();
	for (int i=0;i<L3s.one.size();i++) {
		//cout << index << "\n";
		if (L3s.three[i]>cut) {
			double w1=(L3s.one[i]+0.)/(nl3x+0.);
			double wh=(L3s.two[i]+0.)/(nl3x+0.);
			double w0=(nl3x-w1-wh+0.)/(nl3x+0.);
			if (sel==1) {
				double sumfit=((s1*w1)+(sh*wh)+(s0*w0)+0.);
				w1=(w1*s1)/sumfit;
				wh=(wh*sh)/sumfit;
				w0=(w0*s0)/sumfit;
			}
			double wp=L3s.three[i];
			//Assume deterministic egg probabilities for sheep with L3 worm input
			double p1=(w1*w1)+(w1*wh)+(0.25*wh*wh);
			double ph=(w1*wh)+(2*w1*w0)+(0.5*wh*wh)+(wh*w0);
			double p0=(0.25*wh*wh)+(wh*w0)+(w0*w0);
			eggdist.one.push_back(p1);
			eggdist.two.push_back(ph);
			eggdist.three.push_back(wp);
			if (p1>1) {
				cout << "Error here " << w1 << " " << wh << " " << w0 << " " << p1 << " " << ph << " " << p0 << "\n";
			}
			
		}
	}
}


void GenerateBC4P1 (int nl3x, int nl4x, double cut, double s1, double sh, double s0, L3struc newL3, tvec eggdist, vector<double> logfactorial) {
	//Passage through naive sheep - get BC4P1
	CalculatEggDist (0,nl3x,cut,s1,sh,s0,newL3,eggdist);  //No selection
	
	cout << "Egg distribution\n";
	double tot=0;
	for (int i=0;i<eggdist.one.size();i++) {
		cout << i << " " << eggdist.one.size() << " " << eggdist.one[i] << " " << eggdist.two[i] << " " << eggdist.three[i] << "\n";
		tot=tot+eggdist.three[i];
	}
	cout << "Egg tot " << tot << "\n";
	
	//Collect 100 L4 worms; bottleneck
	
	cout << "Getting BC4P1 here\n";
	L3struc BC4P1;
	GetL3Struc(nl4x,cut,BC4P1,eggdist,logfactorial);
	tot=0;
	for (int i=0;i<BC4P1.one.size();i++) {
		//	cout << BC4P1.one[i] << " " << BC4P1.two[i] << " " << BC4P1.three[i] << "\n";
		tot=tot+BC4P1.three[i];
	}
	cout << "Total " << tot << "\n";
	
	cout << "Combine BC4P1\n";
	CombineL3Struc(BC4P1);
	
	ofstream dat1_file;
	dat1_file.open("BC4P1.dat");
	cout << "Here BC4P1 ";
	cout << BC4P1.one.size() << "\n";
	tot=0;
	for (int i=0;i<BC4P1.one.size();i++) {
		dat1_file << BC4P1.one[i] << " " << BC4P1.two[i] << " " << BC4P1.three[i] << "\n";
		tot=tot+BC4P1.three[i];
	}
	cout << "Total " << tot << "\n";
}

void GenerateBC4P1_IVM (int nl3x, int nl4x, double cut, double s1, double sh, double s0, L3struc newL3, tvec eggdist, vector<double> logfactorial) {
	//Passage through drugged sheep - get BC4P1_IVM
	
	CalculatEggDist (1,nl3x,cut,s1,sh,s0,newL3,eggdist); //With selection
	
	cout << "Egg distribution\n";
	double tot=0;
	for (int i=0;i<eggdist.one.size();i++) {
		cout << eggdist.one[i] << " " << eggdist.two[i] << " " << eggdist.three[i] << "\n";
		tot=tot+eggdist.three[i];
	}
	cout << "Egg tot " << tot << "\n";
	
	cout << "Getting BC4P1_IVM here\n";
	L3struc BC4P1_IVM;
	GetL3Struc(nl4x,cut,BC4P1_IVM,eggdist,logfactorial);
	tot=0;
	for (int i=0;i<BC4P1_IVM.one.size();i++) {
		//	cout << BC4P1.one[i] << " " << BC4P1.two[i] << " " << BC4P1.three[i] << "\n";
		tot=tot+BC4P1_IVM.three[i];
	}
	cout << "Total " << tot << "\n";
	
	
	cout << "Combine BC4P1_IVM\n";
	CombineL3Struc(BC4P1_IVM);
	
	
	ofstream dat2_file;
	dat2_file.open("BC4P1_IVM.dat");
	cout << "Here BC4P1_IVM ";
	cout << BC4P1_IVM.one.size() << "\n";
	tot=0;
	for (int i=0;i<BC4P1_IVM.one.size();i++) {
		dat2_file << BC4P1_IVM.one[i] << " " << BC4P1_IVM.two[i] << " " << BC4P1_IVM.three[i] << "\n";
		tot=tot+BC4P1_IVM.three[i];
	}
	cout << "Total " << tot << "\n";
}







double func_d (double x, double m1, double s1, double m2, double s2) {
    double a=func_a(x,s1,s2);
    double b=func_b(x,m1,s1,m2,s2);
    double c=func_c(m1,s1,m2,s2);
    double d=(pow(b,2)-(c*pow(a,2))/(2*pow(a,2)));
    d=exp(d);
    return d;
}

double emod (double x) {
    double e=0.5*(1+erf(x/sqrt(2)));
    return e;
}


int factorial(unsigned int n) {
    int f=1;
    if (n>1) {
        for (int i=1;i<=n;i++) {
            f=f*i;
        }
    }
    return f;
}

double logfact (unsigned int n) {
    double f=0;
    if (n>1) {
        for (int i=1;i<=n;i++) {
            f=f+log(i);
        }
    }
    return f;
}



double poisson(double rho, int x) {
    double p=exp(-rho)*pow(rho,x);
    p=p/(factorial(x)+0.);
    return p;
}

int imax (int a, int b) {
    if (a>b) {
        return a;
    } else {
        return b;
    }
}

double BetaBinomCalc(int depth, int r, float p) {
	double alpha=depth*p;
	double beta=depth*(1-p);
	double bin=1e10;
	if (alpha==0) {
		if (r==0) {
			bin=1;
		} else {
			bin=0;
		}
	} else if (beta==0) {
		if (r==1) {
			bin=1;
		} else {
			bin=0;
		}
	} else {
		bin=gsl_sf_lngamma(depth+1)+gsl_sf_lngamma(r+alpha)+gsl_sf_lngamma(depth-r+beta)+gsl_sf_lngamma(alpha+beta)-gsl_sf_lngamma(r+1)-gsl_sf_lngamma(depth-r+1)-gsl_sf_lngamma(depth+alpha+beta)-gsl_sf_lngamma(alpha)-gsl_sf_lngamma(beta);
		bin=exp(bin);
	}
	if (bin<1e-300) {
		bin=1e-300;
	}
	return(bin);
}

void GetDist (double my, double mf, double m1, dvec& ptint, int grid) {
    double sy=sqrt(my);
	double s1=sqrt(m1);
    //double nm=func_p(mf,m1,s1,my,sy);
    //Calibrate scale of calculation
    double scale=0.25;
    double next=func_p(mf+(grid*scale),m1,s1,my,sy);
    for (int i=1;i<=1000;i++) {
        if (next<1e-10) {
            scale=scale/1.1;
            next=func_p(mf+(grid*scale),m1,s1,my,sy);
        } else {
            break;
        }
    }
    dvec ptab;
    double x=mf-(grid*scale);
	
	//Check for probability at frequencies less than zero - exclude these...
	if (x<0) {
		x=0;
		scale=(mf+(grid*scale))/(2*(grid+1));
	}
	while (x<=mf+((grid+0.01)*scale)) {
		ptab.one.push_back(func_r(x));
		ptab.two.push_back(scale*func_p(x,m1,s1,my,sy));
		x=x+scale;
	}
	
	double tot=0;
	for (unsigned int i=0;i<ptab.one.size();i++) {
		tot=tot+ptab.two[i];
	//	cout << ptab.one[i] << " " << ptab.two[i] << "\n";
	}
	
	//cout << "Total " << tot << "\n";

	//Simpson's rule integration
//	double tot=0;
//    for (unsigned int i=0;i<=ptab.one.size()-2;i=i+2) {
//        ptint.one.push_back(ptab.one[i+1]);
//        double t=(scale/6)*(ptab.two[i]+(4*ptab.two[i+1])+ptab.two[i+2]);
//        ptint.two.push_back(t);
//		tot=tot+t;
//    }
	
//	cout << "Total " << tot << "\n";
	
	//Check for excluded probability from negative frequencies - above
//	if (tot<1) {
//		ptint.two[0]=ptint.two[0]+(1-tot);
//	}
	ptint.one=ptab.one;
	ptint.two=ptab.two;
	
}

void ConstructProbs (dvec ptint1, dvec ptint2, dvec ptint3, vvec& probs) {
	double p;
    double p1;
    double p2;
    double p3;
    double pt;
	probs.one.clear();
	probs.two.clear();
	probs.three.clear();
	probs.four.clear();
	probs.five.clear();
	for (unsigned int i=0;i<ptint1.one.size();i++) {
        for (unsigned int j=0;j<ptint2.one.size();j++) {
            for (unsigned int k=0;k<ptint3.one.size();k++) {
                p=1;
                p3=p*ptint3.one[k];
                p2=ptint2.one[j]*p*(ptint3.one[k]+1);
                p1=ptint1.one[i]*p*((ptint2.one[j]*(ptint3.one[k]+1))+ptint3.one[k]+1);
                pt=p+p1+p2+p3;
                probs.one.push_back(p1/pt);
                probs.two.push_back(p2/pt);
                probs.three.push_back(p3/pt);
                probs.four.push_back(p/pt);
                probs.five.push_back(ptint1.two[i]*ptint2.two[j]*ptint3.two[k]);
            }
        }
    }
	//cout << "Size of probs is " << probs.one.size() << "\n";

	//Reduce probs to those >10^-10
	vvec probs2;
	double tot2=0;
    for (unsigned int i=0;i<probs.one.size();i++) {
        if (probs.five[i]>1e-10) {
            probs2.one.push_back(probs.one[i]);
            probs2.two.push_back(probs.two[i]);
            probs2.three.push_back(probs.three[i]);
            probs2.four.push_back(probs.four[i]);
            probs2.five.push_back(probs.five[i]);
			tot2=tot2+probs.five[i];
        }
    }
	//cout << "Total here is " << tot2 << "\n";
    //cout << "New size of probs is " << probs2.one.size() << "\n";
    probs=probs2;
}

void calc_mprobs (int i, double cut, double pre_fac, vvec probs, map<long long,double>& mposs, vector<double> logfactorial) {
	//Find approx closest solution and search outwards...
    int a1=floor(50*probs.one[i]+0.5);
    int a2=floor(50*probs.two[i]+0.5);
    int a3=floor(50*probs.three[i]+0.5);
    int a4=50-a1-a2-a3;
    //cout << "Guess " << a1 << " " << a2 << " " << a3 << " " << a4 << "\n";
    vector<double> prb;
    prb.push_back(probs.one[i]);
    prb.push_back(probs.two[i]);
    prb.push_back(probs.three[i]);
    prb.push_back(probs.four[i]);
    double m;
    double chk1=-1;
    double chk2=-1;
	double chk3=-1;
	int d=0;
	long long s;
	vector<int> cts;
	int a;
	int b;
	int c;
    for (a=a1;a<=50;a++) {
        for (b=a2;b<=50-a;b++) {
            for (c=a3;c<=50-(a+b);c++) {
				//minternal3(1,a,b,c,chk1);
                d=50-a-b-c;
				cts.clear();
				cts.push_back(a);
                cts.push_back(b);
                cts.push_back(c);
                cts.push_back(d);
				//cout << "Prb " << prb[0] << " " << prb[1] << " " << prb[2] << " " << prb[3] << "\n";
				//cout << a << " " << b << " " << c << " " << d << " " << pre_fac << " " << probs.five[i] << " ";
                m=pre_fac*multinomp(50,cts,prb,logfactorial)*probs.five[i];
                //cout  << m << "\n";
				s=1e9*a+1e6*b+1000*c+d;
				mposs[s]=mposs[s]+m;
                if (m<cut) {
                    //cout << "Small\n";
                    chk1=c;
                    break;
                }
            }
            if (chk1!=a3&&a3>0) {
                for (c=a3-1;c>=0;c--) {
                    d=50-a-b-c;
					cts.clear();
                    cts.push_back(a);
                    cts.push_back(b);
                    cts.push_back(c);
                    cts.push_back(d);
                    m=pre_fac*multinomp(50,cts,prb,logfactorial)*probs.five[i];
					//cout << a << " " << b << " " << c << " " << d << " " << m << "\n";
					s=1e9*a+1e6*b+1000*c+d;
					mposs[s]=mposs[s]+m;
					if (m<cut) {
						//    cout << "Small\n";
						//    chk1=c;
                        break;
                    }
                }
            }
            if (chk1==a3&&m<cut) {
                //cout << "Small2\n";
                chk2=b;
                break;
            }
        }
		if (chk2!=a2&&a2>0) {
			//cout << "Repeat starts here\n";
			for (b=a2-1;b>=0;b--) {
				for (c=a3;c<=50-(a+b);c++) {
					d=50-a-b-c;
					cts.clear();
					cts.push_back(a);
					cts.push_back(b);
					cts.push_back(c);
					cts.push_back(d);
					m=pre_fac*multinomp(50,cts,prb,logfactorial)*probs.five[i];
					//cout << a << " " << b << " " << c << " " << d << " " << m << "\n";
					s=1e9*a+1e6*b+1000*c+d;
					mposs[s]=mposs[s]+m;
					if (m<cut) {
						//cout << "Small\n";
						chk1=c;
						break;
					}
				}
				if (chk1!=a3&&a3>0) {
					for (c=a3-1;c>=0;c--) {
						d=50-a-b-c;
						cts.clear();
						cts.push_back(a);
						cts.push_back(b);
						cts.push_back(c);
						cts.push_back(d);
						m=pre_fac*multinomp(50,cts,prb,logfactorial)*probs.five[i];
						//cout << a << " " << b << " " << c << " " << d << " " << m << "\n";
						s=1e9*a+1e6*b+1000*c+d;
						mposs[s]=mposs[s]+m;
						if (m<cut) {
							break;
						}
					}
				}
			}
			//cout << "Repeat done here\n";
		}
        if (chk2==a2&&m<cut) {
			//            cout << "Small3\n";
			chk3=a;
            break;
        }
    }
	//    cout << "chk3 " << chk3 << "\n";
	if (chk3!=a1&&a1>0) {
		for (a=a1-1;a>0;a--) {
			for (b=a2;b<=50-a;b++) {
				for (c=a3;c<=50-(a+b);c++) {
					//minternal3(1,a,b,c,chk1);
					d=50-a-b-c;
					cts.clear();
					cts.push_back(a);
					cts.push_back(b);
					cts.push_back(c);
					cts.push_back(d);
					m=pre_fac*multinomp(50,cts,prb,logfactorial)*probs.five[i];
					//cout << a << " " << b << " " << c << " " << d << " " << m << "\n";
					s=1e9*a+1e6*b+1000*c+d;
					mposs[s]=mposs[s]+m;
					if (m<cut) {
						//cout << "Small\n";
						chk1=c;
						break;
					}
				}
				if (chk1!=a3) {
					for (c=a3-1;c>=0;c--) {
						d=50-a-b-c;
						cts.clear();
						cts.push_back(a);
						cts.push_back(b);
						cts.push_back(c);
						cts.push_back(d);
						m=pre_fac*multinomp(50,cts,prb,logfactorial)*probs.five[i];
						//cout << a << " " << b << " " << c << " " << d << " " << m << "\n";
						s=1e9*a+1e6*b+1000*c+d;
						mposs[s]=mposs[s]+m;
						if (m<cut) {
							//    cout << "Small\n";
							//    chk1=c;
							break;
						}
					}
				}
				if (chk1==a3&&m<cut) {
					//cout << "Small2\n";
					chk2=b;
					break;
				}
			}
			if (chk2!=a2) {
				//cout << "Repeat starts here\n";
				for (b=a2-1;b>=0;b--) {
					for (c=a3;c<=50-(a+b);c++) {
						d=50-a-b-c;
						cts.clear();
						cts.push_back(a);
						cts.push_back(b);
						cts.push_back(c);
						cts.push_back(d);
						m=pre_fac*multinomp(50,cts,prb,logfactorial)*probs.five[i];
						//cout << a << " " << b << " " << c << " " << d << " " << m << "\n";
						s=1e9*a+1e6*b+1000*c+d;
						mposs[s]=mposs[s]+m;
						if (m<cut) {
							//cout << "Small\n";
							chk1=c;
							break;
						}
					}
					if (chk1!=a3) {
						for (c=a3-1;c>=0;c--) {
							d=50-a-b-c;
							cts.clear();
							cts.push_back(a);
							cts.push_back(b);
							cts.push_back(c);
							cts.push_back(d);
							m=pre_fac*multinomp(50,cts,prb,logfactorial)*probs.five[i];
							//cout << a << " " << b << " " << c << " " << d << " " << m << "\n";
							s=1e9*a+1e6*b+1000*c+d;
							mposs[s]=mposs[s]+m;
							if (m<cut) {
								break;
							}
						}
					}
				}
				//cout << "Repeat done here\n";
			}
			if (chk2==a2&&m<cut) {
				//cout << "Small3\n";
				chk3=a;
				break;
			}
		}
	}
}

/*
void GetParameters (run_params& p, int argc, const char **argv) {
	string p_switch;
	int x=1;
	p.backcross=0;
	p.N=50;
	p.L=101;
	p.tSim=4;
	p.sigma=0;
	p.rho=0.01;
	p.mu=1e-6;
    p.eggs=500;
	p.seed=(int) time(NULL);
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		cout << p_switch << "\n";
		if (p_switch.compare("--backcross")==0) {
			x++;
			p.backcross=atoi(argv[x]);
		}
		else if (p_switch.compare("--N")==0) {
			x++;
			p.N=atoi(argv[x]);
		}
		else if (p_switch.compare("--L")==0) {
			//cout << "Here\n";
			x++;
			p.L=atoi(argv[x]);
		}
		else if (p_switch.compare("--tSim")==0) {
			x++;
			p.tSim=atoi(argv[x]);
		}
		else if (p_switch.compare("--seed")==0) {
			x++;
			p.seed=atoi(argv[x]);
		}
		else if (p_switch.compare("--sigma")==0) {
			x++;
			p.sigma=atof(argv[x]);
		}
		else if (p_switch.compare("--rho")==0) {
			x++;
			p.rho=atof(argv[x]);
		}
		else if (p_switch.compare("--mu")==0) {
			x++;
			p.mu=atof(argv[x]);
		}
        else if (p_switch.compare("--eggs")==0) {
            x++;
            p.eggs=atof(argv[x]);
        }
		else {
			cout << "Incorrect usage\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
}
*/


