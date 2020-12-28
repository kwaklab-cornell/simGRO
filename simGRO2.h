#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <cmath>
#include <list>
#include "kernel.h"

using namespace std;
// Statefunction
#define MAX_STATE 20

// Copies of DNA templates
int POP_N = 1000;
#define MAXCOPY	20000
#define MAXOBJ 10000
#define DNALEN 1200
#define TSSPOS 200
#define PASPOS 800
#define PAUSESITE 240

// Promoter status bit
#define N_PRSTATUS	1
#define IC			0
// Pol2 object types
#define POL2		1
#define POL2PF		2
// Pol II definition
#define PROTECTION5	20
#define	PROTECTION3	15	
#define ICPOS	200       // Position for initiacion complex including TBP and TFIID

// Pol II status bit
#define N_POL2STATUS 6
#define PIC	        0
#define ELONGATING  1
#define PAUSED		2
#define POSTPOLYA   3
#define PFBOUND     4
#define POL2AGE		5

struct Pol2stat
{
	int Position;		// position on DNA template
	int Status[N_POL2STATUS];		// Pol II status
	void Init(int TSSpos);
	void Elongate(int delta_t,int *Factor_occupancy);		// Position change after delta_t millisecond
	int CheckPausing(int delta_t,int *Factor_occupancy);        // Returns 1 for pausing event
    int CheckEscape(int delta_t);      // Returns 1 for pause escape event
	int CheckCleavage(int delta_t);      // Returns 1 for post poly A cleavage  event
	int CheckTermination(int delta_t);		// Returns 1 for termination event
};

struct DNAstat
{
	int Length;
	int TSS;
	int PAS;
	int *Factor_occupancy;
	list<Pol2stat> Pol2list;
	int Promoter_status[N_PRSTATUS];
	void Init(int len, int TSSpos, int PASpos);
	void Progress(int delta_t);
	void ObjectList(int &n, int *objtype, int *objage, int *objpos);
	void Addmetagene(float *meta);
	void CheckIC(int delta_t);
};

struct StateFunction
{
	float State_mat[1000][MAX_STATE];	// Cummulative probability interval at x millisec  
	void Init(char *def_file);
	void InitPoisson(float average_occurrance);	// Initialize Poisson probability function per sec
	int StateChange(int delta_t);	// Probablistic assignment of state change
};

int Assay_res_type = 1;
DNAstat DNArep;     // representative DNA status
DNAstat DNAtem[MAXCOPY];    // population DNA status
StateFunction IC_on;        // initiation complex on rate
StateFunction IC_off;       // initiation complex off rate
StateFunction Recruitment;  // Pol II recruitment rate (including PIC formation and initial open complex formation)
StateFunction Elongation;   // Pol II elongation rate
StateFunction PF_on;        // Pause factor on rate
StateFunction PF_off;       // Pause factor off rate
StateFunction PA_Cleavage;   // PolyA cleavage rate
StateFunction Elongation_postPA;    // Post-polyA elongation rate
StateFunction Termination;; // termination rate
float Metagene[DNALEN];
kernel kern[3] = {kernel(1),kernel(5),kernel(25)};
int IC_on_hl;
int IC_off_hl;
int PF_off_hl;

// Function implementation


void ResetSimGRO()
{
	IC_on.InitPoisson(.2);
    IC_off.InitPoisson(.5);
    Recruitment.InitPoisson(1);
    Elongation.InitPoisson(30);
    PF_on.InitPoisson(10);
    PF_off.InitPoisson(.5);
    PA_Cleavage.InitPoisson(1);
    Elongation_postPA.InitPoisson(15);
    Termination.InitPoisson(0.2);
	IC_on_hl=5;
	IC_off_hl=2;
	PF_off_hl=2;
}
	
void InitSimGRO()
{
	DNArep.Init(DNALEN,TSSPOS,PASPOS);
	for(int i=0;i<MAXCOPY;i++) DNAtem[i].Init(DNALEN,TSSPOS,PASPOS);
	srand(time(NULL));
	ResetSimGRO();
}

int stcompare(char *s1, string s2)
{
	int len1, len2;
	len1=strlen(s1);
	len2=s2.length();
	
	if(len1<len2) return 1;
	
	for(int i=0;i<len2;i++)
	{
		if(!(s1[i]==s2[i]||s1[i]==s2[i]-32||s1[i]==s2[i]+32)) return 1;
	}
	
	return 0;
}

void ResetSimGRO(char *pref_file)
{
	ifstream pref;
	pref.open(pref_file);

	char buf[1024];
	float value;
	while(!pref.eof())
	{
		pref>>buf;
		if(stcompare(buf,"#")==0) pref.getline(buf,1024);
		else {
			pref>>value;
			if(stcompare(buf,"Recruitment")==0) Recruitment.InitPoisson(value);
			else if(stcompare(buf,"Elongation")==0) Elongation.InitPoisson(value);
			else if(stcompare(buf,"InitComplex_on")==0) { IC_on.InitPoisson(value); IC_on_hl=(int)(1/value+0.5); }
			else if(stcompare(buf,"InitComplex_off")==0) { IC_off.InitPoisson(value); IC_off_hl=(int)(1/value+0.5); }
			else if(stcompare(buf,"PauseFactor_on")==0) PF_on.InitPoisson(value);		
			else if(stcompare(buf,"PauseFactor_off")==0) { PF_off.InitPoisson(value); PF_off_hl=(int)(1/value+0.5); }
			else if(stcompare(buf,"PolyA_cleavage")==0) PA_Cleavage.InitPoisson(value);
			else if(stcompare(buf,"PostPA_Elongation")==0) Elongation_postPA.InitPoisson(value);
			else if(stcompare(buf,"Termination")==0) Termination.InitPoisson(value);
		}
	}
}

float tempmeta[DNALEN];

void UpdateSimGRO(int delta_t)
{
	DNArep.Progress(delta_t);
	for(int i=0;i<DNALEN;i++) tempmeta[i]=0;
	for(int i=0;i<POP_N;i++)
	{
		DNAtem[i].Progress(delta_t);
		DNAtem[i].Addmetagene(tempmeta);
	}

	float tempmodfactor[1];
	kern[Assay_res_type].addto(1,0,tempmodfactor,1);
	kern[Assay_res_type].addto(1,PROTECTION3+PROTECTION5,tempmodfactor,1);
	kern[Assay_res_type].addto(1,-PROTECTION3-PROTECTION5,tempmodfactor,1);
	
	for(int i=0;i<DNALEN;i++) Metagene[i]=0;
    for(int i=0;i<DNALEN;i++) kern[Assay_res_type].addto(tempmeta[i]/tempmodfactor[0],i,Metagene,DNALEN);		// Kernel smoothing
}

// Debugging code
int prev_time=0;
		
void PrintStatus(int time, int interval)
{
	cout<<"Elapsed time = "<<time/1000<<"\tFrames per second = "<<(float)1000/interval<<endl;
	float pp=0, gb=0;
	for(int i=DNArep.TSS;i<DNArep.TSS+100;i++) pp+=tempmeta[i]/POP_N;
	for(int i=DNArep.TSS+250;i<DNArep.PAS-250;i++) gb+=tempmeta[i]/POP_N/(DNArep.PAS-DNArep.TSS-500)*1000;
	cout<<"\t\t\tPromoter density = "<<pp<<" (molecules per promoter region from TSS to +100bp)"<<endl;
	cout<<"\t\t\tGene body density = "<<gb<<" (molecules per copy per kb)"<<endl;
	cout<<"\t\t\tPausing index = ";
	if(gb==0) cout<<"inf"<<endl;
	else cout<<pp/gb*10<<endl;
}

void Print(int time, int interval, int saveframe)
{
	if(time/1000>=prev_time+2)
	{
		PrintStatus(time, interval);
		prev_time = time/1000;
		if(saveframe>10000000) cout<<"Saving movie at frame "<<saveframe-10000000<<endl;
	}
}

void Pol2stat::Init(int TSSpos = TSSPOS)
{
	Position = TSSpos;
	for(int i=0;i<N_POL2STATUS;i++) Status[i] = 0;
	Status[ELONGATING]=1;
	Status[POL2AGE]=1;
}

void Pol2stat::Elongate(int delta_t, int *Factor_occupancy)
{
    if(Status[PIC]||Status[PAUSED]) return;
    else if(Status[ELONGATING])
	{
		int pr;     // progression length
		if(Status[POSTPOLYA]) pr=Elongation_postPA.StateChange(delta_t);
		else pr=Elongation.StateChange(delta_t);
		Factor_occupancy[Position]&=(~POL2);		// remove Pol II occupancy at present location;
		for(int i=Position+pr-PROTECTION5-PROTECTION3;i<Position+pr+PROTECTION3+PROTECTION5;i++)
        if(i>=0&&i<DNALEN) if(Factor_occupancy[i]&POL2)		// if new position is occupied, then don't advance
        {
            pr=0;
            break;
        }
		Position+=pr;
		if(Position>=DNALEN) Position=DNALEN-1;
		Factor_occupancy[Position]|=POL2;		// mark new positions occupied by pol2
		return;	
	}
}

int Pol2stat::CheckPausing(int delta_t,int *Factor_occupancy)
{
    if(Position<=PAUSESITE)             // check for pause factor binding
    {
        if((Status[PFBOUND]<=0)&&!Status[PAUSED]) 
        {
            if(PF_on.StateChange(delta_t)) Status[PFBOUND]=1;         // pause factor gets bound
        }
    }
    if(Position>PAUSESITE)          // force the pol2 to be at pause site
    {
        if(Status[PFBOUND]>0)
        {
            Factor_occupancy[Position]&=(~POL2);
			Position=PAUSESITE;
            Status[PAUSED]=1;
			Factor_occupancy[Position]|=POL2;
            return 1;
        }
    }
    return 0;
}

int Pol2stat::CheckEscape(int delta_t)
{
    if(Status[PAUSED])
    {
        if(PF_off.StateChange(delta_t))
        {
            Status[PFBOUND]=-8;
            Status[PAUSED]=0;
			Status[ELONGATING]=1;
            return 1;
        }
    }
    return 0;
}

int Pol2stat::CheckCleavage(int delta_t)
{
    if((Position>PASPOS)&&!Status[POSTPOLYA])
    {
        if(PA_Cleavage.StateChange(delta_t))
        {
            Status[POSTPOLYA]=1;
            return 1;
        }
    }
    return 0;
}


int Pol2stat::CheckTermination(int delta_t)
{
	if(Status[POSTPOLYA])
	{
		if(Termination.StateChange(delta_t)) return 1;
	}
	return 0;
}

void DNAstat::Init(int len, int TSSpos, int PASpos)
{
	Length = len;
	TSS = TSSpos;
	PAS = PASpos;
	
	Factor_occupancy = new int[Length];
	for(int i=0;i<Length;i++) Factor_occupancy[i]=0;
	for(int i=0;i<N_PRSTATUS;i++) Promoter_status[i]=0;
}

void DNAstat::CheckIC(int delta_t)
{
	// TBP bound state
	if(Promoter_status[IC]<=0)		// if IC is not bound
	{
		if(IC_on.StateChange(delta_t)) Promoter_status[IC]=1;
	}
	else		// IC bound
	{
		if(IC_off.StateChange(delta_t)) Promoter_status[IC]=-8;
	}
}

void DNAstat::Progress(int delta_t)
{
	CheckIC(delta_t);

	// Generate new Pol2
	int pol2_clear=1;
	for(int i=TSS;i<=PAUSESITE;i++) if(Factor_occupancy[i]&POL2) pol2_clear=0;
	if(Promoter_status[IC]>0&&pol2_clear)		// IC bound and No Pol2 at promoter
		if(Recruitment.StateChange(delta_t))
		{
			Pol2stat newpol2;
			newpol2.Init(TSS);
			Pol2list.push_front(newpol2);
			Factor_occupancy[TSS]|=POL2;		// mark TSS occupied by pol2
		}
	for(int i=0;i<N_PRSTATUS;i++) if(Promoter_status[i]!=0&&Promoter_status[i]!=8) Promoter_status[i]++;

	// Elongate Pol2 in the list
	for(list<Pol2stat>::iterator it = Pol2list.begin();it!=Pol2list.end();it++)
	{
		it->Elongate(delta_t,Factor_occupancy);
		it->CheckPausing(delta_t,Factor_occupancy);
		it->CheckEscape(delta_t);
		it->CheckCleavage(delta_t);
		for(int i=0;i<N_POL2STATUS;i++) if(it->Status[i]!=0&&it->Status[i]!=8) it->Status[i]++;

		if(it->CheckTermination(delta_t))		// Terminating Pol 2
		{
			it->Status[POL2AGE]=-8;			// mark as terminating when age reaches 0
		}
		if(it->Position>DNALEN-50||it->Status[POL2AGE]==0)
		{			
			Factor_occupancy[it->Position]&=(~POL2);		// terminate Pol2 and erase the mark
			Pol2list.erase(it);
		}
	}
}

void DNAstat::ObjectList(int &n, int *objtype, int *objage, int *objpos)
{
	list<Pol2stat>::iterator it;

	//Count number of objects
	n=0;
	if(Promoter_status[IC]!=0)
	{
		objtype[n] = IC;
		objage[n] = Promoter_status[IC];
		objpos[n] = TSS;
		n++;
	}
	// Pol 2 position list
	it = Pol2list.begin();
	while(it!=Pol2list.end()&&n<MAXOBJ)
	{
		objtype[n] = POL2;
		objage[n]=it->Status[POL2AGE];
		objpos[n] = it->Position;
		if(it->Status[PFBOUND]!=0)
		{
			n++;
			if(n==MAXOBJ) break;
			objtype[n] = POL2PF;
			objage[n] = it->Status[PFBOUND];
			objpos[n] = it->Position;
		}
		it++;
		n++;
	}
	
}

void DNAstat::Addmetagene(float *meta)
{
	for(list<Pol2stat>::iterator it = Pol2list.begin();it!=Pol2list.end();it++)
	{
		meta[it->Position]+=(float)1/POP_N;
	}
}

void StateFunction::Init(char *def_file)
{
	return;
}

void StateFunction::InitPoisson(float average_occurrance)
{
	float lambda;
	float f;
		
	for(int i=0;i<1000;i++)
	{
		lambda = average_occurrance*i/1000;
		for(int j=0;j<MAX_STATE;j++)
		{
			f = pow(lambda,j)*exp(-lambda);
			for(int k=1;k<=j;k++) f/=k;
			if(j==0) State_mat[i][j]=f;
			else State_mat[i][j]=State_mat[i][j-1]+f;
		}			
	}
	return;
}

int StateFunction::StateChange(int delta_t)
{
	float r;

	r = (float)rand()/RAND_MAX+(float)rand()/RAND_MAX/RAND_MAX;
	int state = 0;
	if(delta_t>=1000) delta_t=999;

	while(state<MAX_STATE)
	{
		if(r<State_mat[delta_t][state]) break;
		state++;
	}
	
	return state;
}
