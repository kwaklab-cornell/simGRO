#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <cmath>
#include <list>
#include <vector>

using namespace std;
// Statefunction
#define MAX_STATE 20
#define MAXRECORD 100

// Copies of DNA templates
int POP_N = 1000;
#define MAXCOPY	20000
#define MAXOBJ 50000
int DNALEN = 12000;
int TSSPOS = 200;
int PASPOS = 10200;
int PAUSESITE = 240;

// Promoter status bit
#define N_PRSTATUS	1
#define IC			0

// Pol2 object types
#define POL2		1
#define POL2PF		2

// Pol II definition
#define PROTECTION5	20
#define	PROTECTION3	15	
int ICPOS = 200;       // Position for initiation complex including TBP and TFIID

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
	int Activity;
	void Init();
	void Elongate(int delta_t,int *Factor_occupancy);		// Position change after delta_t millisecond
	int CheckPausing(int delta_t,int *Factor_occupancy);        // Returns 1 for pausing event
    int CheckEscape(int delta_t);      // Returns 1 for pause escape event
	int CheckCleavage(int delta_t);      // Returns 1 for post poly A cleavage  event
	int CheckTermination(int delta_t);		// Returns 1 for termination event
};

struct DNAstat
{
	int *Factor_occupancy;
	list<Pol2stat> Pol2list;
	int Promoter_status[N_PRSTATUS];
	void Init();
	void Progress(int delta_t);
	void ObjectList(int &n, int *objtype, int *objage, int *objpos);
	void Addmetagene(float *meta);
	void CheckIC(int delta_t);
};

// StateFunction: object that returns the state change based on a Poisson probability event
// defined by average event per second after a time interval delta_t
struct StateFunction
{
	float State_mat[1000][MAX_STATE];	// Cummulative probability interval at x millisec  
	void InitPoisson(float average_occurrance);	// Initialize Poisson probability function per sec
	int StateChange(int delta_t);	// Probablistic assignment of state change
};

// FunctionOf: two dimensional function, with respect to time and position to intrapolate
// dynamic event frequency
struct FunctionOf
{
	vector<vector<float> > data;
	vector<float> v1, v2;
	void Init(char *fileName);
	float value(float p1, float p2);
	float min, max;
};

// StateFunctionOf: scaffold of multiple StateFunctions to calculate state change efficiently
// given delta_t and average event frequency per second
struct StateFunctionOf
{
	float min,max;
	float step;
	vector<StateFunction> stateFunc;
	void Init(float st, char *fileName=NULL);
	int StateChange(int delta_t, float value1, float value2);
	FunctionOf rateFunc;
};

DNAstat DNArep;     // representative DNA status
DNAstat *DNAtem;    // population DNA status

StateFunctionOf IC_on;        // initiation complex on rate
StateFunctionOf IC_off;       // initiation complex off rate
StateFunctionOf Recruitment;  // Pol II recruitment rate (including PIC formation and initial open complex formation)
StateFunctionOf PF_on;        // Pause factor on rate
StateFunctionOf PF_off;        // Pause factor off rate
StateFunctionOf Elongation;   // Pol II elongation rate
StateFunctionOf PA_Cleavage;   // PolyA cleavage rate
StateFunctionOf Elongation_postPA;    // Post-polyA elongation rate
StateFunctionOf Termination; // termination rate
StateFunctionOf ActivityDistribution;	// Intrinsic Pol II activity determined at initiation

float *Metagene, *tempmeta;
float RecordTime[MAXRECORD];
int currentRecord=0;
int countRecord=0;
ofstream recordFile;

int simTime = 0;

// Function implementation

void ResetSimGRO(char *pref_file);
void ResetSimGRO()
{
	IC_on.Init(0.2);				// TF on every 5 s
	IC_off.Init(0.5);				// TF off every 2 s
    Recruitment.Init(1);			// Pol II can fire every 1 s
	PF_on.Init(10);					// Pause factor in 0.1 s
	PF_off.Init(0.5);				// Pause factor off every 2 s
	Elongation.Init(30);			// Elongation speed at 30 bp/s
    PA_Cleavage.Init(1);			// Poly(A) cleavage once every 1 s	
    Elongation_postPA.Init(15);		// Post-poly(A) elongation 15 bp/s
	Termination.Init(0.2);			// Termination rate once every 5 s
	ActivityDistribution.Init(100);	// Pol II acitivity 
}
	
void InitSimGRO(char *pref_file, char *recordFileName)
{
	//cerr<<"Initializing"<<endl;
	ResetSimGRO();
	//cerr<<"Constants initialized"<<endl;
	ResetSimGRO(pref_file);
	//cerr<<"Constants loaded"<<endl;
	DNArep.Init();
	//cerr<<"Representative DNA template initialized"<<endl;
	DNAtem = new DNAstat[MAXCOPY];
	for(int i=0;i<MAXCOPY;i++) DNAtem[i].Init();
	//cerr<<"DNA templates initialized"<<endl;
	srand(time(NULL));
	Metagene = new float[DNALEN];
	tempmeta = new float[DNALEN];
	//cerr<<"Output array initialized"<<endl;
	recordFile.open(recordFileName);
	//cerr<<"Record file opened"<<endl;
}

int stcompare(char *s1, string s2)
{
	int len1, len2;
	len1=strlen(s1);
	len2=s2.length();
	
	if(len1<len2) return 1;
	
	for(int i=0;i<len2;i++)
	{
		if(tolower(s1[i])!=tolower(s2[i])) return 1;
	}
	
	return 0;
}


void ResetSimGRO(char *pref_file)
{
	ifstream pref;
	pref.open(pref_file);
	//cerr<<"Opened preference file"<<endl;
	char buf[1024];
	char buf2[1024];
	float value;
	countRecord=0;
	int geneLen = DNALEN - TSSPOS - PASPOS;
	int postpaLen = DNALEN - PASPOS;
	int prmLen = TSSPOS;
	int psDist = PAUSESITE - TSSPOS;
	while(!pref.eof())
	{
		buf[0]='\0';
		pref>>buf;
		//cerr<<"Reading preference file line: "<<buf<<endl;
		if(stcompare(buf,"#")==0) {
			//cerr<<"Annotation line : "<<buf<<endl;
			pref.getline(buf,1024); }
		else {
			pref>>buf2;
			//cerr<<"Read value "<<buf2<<endl;
			if(stcompare(buf,"TF_on")==0) {
				//cerr<<"Initialozing TF on rates"<<endl;
				IC_on.Init(200, buf2); }
			else if(stcompare(buf,"TF_off")==0) IC_off.Init(200, buf2);
			else if(stcompare(buf,"Recruitment")==0) Recruitment.Init(200, buf2);
			else if(stcompare(buf,"PF_on")==0) PF_on.Init(200, buf2);
			else if(stcompare(buf,"PF_off")==0) PF_off.Init(200, buf2);
			else if(stcompare(buf,"Elongation")==0) Elongation.Init(200, buf2);
			else if(stcompare(buf,"Cleavage")==0) PA_Cleavage.Init(200, buf2);
			else if(stcompare(buf,"PostPA")==0) Elongation_postPA.Init(200, buf2);
			else if(stcompare(buf,"Termination")==0) Termination.Init(200, buf2);
			else if(stcompare(buf,"Activity")==0) ActivityDistribution.Init(200, buf2);
			else if(stcompare(buf,"Record")==0)
			{
				ifstream in(buf2);
				while(true)
				{
					float v=-2e37;
					in>>v;
					if(v<-1e37) break;
					RecordTime[countRecord]=v;
					countRecord++;
				}
				in.close();
				//cerr<<"Read recording times from "<<buf2<<endl;
			}
			else if(stcompare(buf,"GeneLength")==0) {
				stringstream in(buf2);
				in >> geneLen;
				//cerr<<"GeneLength = "<<geneLen<<endl;
			}
			else if(stcompare(buf,"PostPolyALength")==0) {
				//cerr<<"Assigning postPA length"<<endl;
				stringstream in(buf2);
				in >> postpaLen;
				//cerr<<"PostPALength = "<<postpaLen<<endl;
			}
			else if(stcompare(buf,"PromoterLength")==0) {
				stringstream in(buf2);
				in >> prmLen;
				//cerr<<"PromoterLength = "<<prmLen<<endl;
			}
			else if(stcompare(buf,"PauseSite")==0) {
				stringstream in(buf2);
				in >> psDist;
				//cerr<<"PauseSite = "<<psDist<<endl;
			}
			else if(stcompare(buf,"Copies")==0) {
				//cerr<<"Assigning simulation scale"<<endl;
				stringstream in(buf2);
				in >> POP_N;
				//cerr<<"Copies = "<<POP_N<<endl;
			}
		}
	}
	DNALEN = geneLen + postpaLen + prmLen;
	TSSPOS = prmLen;
	PASPOS = geneLen - postpaLen;
	PAUSESITE = prmLen + psDist;
	currentRecord=0;
}


void UpdateSimGRO(int delta_t)
{
	simTime += delta_t;
	DNArep.Progress(delta_t);
	for(int i=0;i<DNALEN;i++) tempmeta[i]=0;
	for(int i=0;i<POP_N;i++)
	{
		DNAtem[i].Progress(delta_t);
		DNAtem[i].Addmetagene(tempmeta);
	}

	for(int i=0;i<DNALEN;i++) Metagene[i]=tempmeta[i];
}

double drand()
{
    return (double)rand()/RAND_MAX+(double)rand()/RAND_MAX/RAND_MAX;
}

// Debugging code
int prev_time=0;
		
void PrintStatus(int time, int interval)
{
	cout<<"Elapsed time = "<<time/1000<<"\tFrames per second = "<<(float)1000/interval<<endl;
	cout<<"Simulated time = "<<simTime/1000<<" seconds"<<endl;
	float pp=0, gb=0;
	for(int i=TSSPOS;i<TSSPOS+100;i++) pp+=tempmeta[i]/POP_N;
	for(int i=TSSPOS+250;i<PASPOS-250;i++) gb+=tempmeta[i]/POP_N/(PASPOS-TSSPOS-500)*1000;
	cout<<"\t\t\tPromoter density = "<<pp<<" (molecules per promoter region from TSS to +100bp)"<<endl;
	cout<<"\t\t\tGene body density = "<<gb<<" (molecules per copy per kb)"<<endl;
	cout<<"\t\t\tPausing index = ";
	if(gb==0) cout<<"inf"<<endl;
	else cout<<pp/gb*10<<endl;
}

void Print(int time, int interval, int bin)
{
	if(time/1000>=prev_time+2)
	{
		PrintStatus(time, interval);
		prev_time = time/1000;
	}
	if(currentRecord<countRecord) 
	if(simTime/1000>RecordTime[currentRecord])
	{
		float binsum = 0;
		for(int i=0;i<DNALEN;i++) {
			binsum += Metagene[i];
			if(i%bin == bin-1) {
				recordFile<<RecordTime[currentRecord]<<"\t"<<i-bin+1<<"\t"<<i<<"\t"<<binsum<<endl;
				binsum = 0;
			}
		}
		currentRecord++;
	}
}

void Pol2stat::Init()
{
	Position = TSSPOS;
	for(int i=0;i<N_POL2STATUS;i++) Status[i] = 0;
	Status[ELONGATING]=1;
	Status[POL2AGE]=1;
	Activity = ActivityDistribution.StateChange(999, 0, simTime);
}

void Pol2stat::Elongate(int delta_t, int *Factor_occupancy)
{
    if(Status[PIC]||Status[PAUSED]) return;
    else if(Status[ELONGATING])
	{
		int pr;     // progression length
		if(Status[POSTPOLYA])
			pr=Elongation_postPA.StateChange(delta_t, Position, Activity);
		else pr=Elongation.StateChange(delta_t, Position, Activity);
		Factor_occupancy[Position]&=(~POL2);		// remove Pol II occupancy at present location;
		for(int i=Position+pr-PROTECTION5-PROTECTION3;i<Position+pr+PROTECTION3+PROTECTION5;i++)
        	if(i>=0&&i<DNALEN) if(Factor_occupancy[i]&POL2)		// if new position is occupied, then don't advance
        	{
            	pr=0;
            	break;
			}
		Position+=pr;
		if(Position>=DNALEN)
		{
			Position=DNALEN-1;
			Status[POL2AGE]=-8;			// mark as terminating when age reaches 0
		}
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
            if(PF_on.StateChange(delta_t, 0, (int)(simTime/1000))) Status[PFBOUND]=1;         // pause factor gets bound
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
        if(PF_off.StateChange(delta_t, 0, (int)(simTime/1000)))
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
    if(!Status[POSTPOLYA])
    {
        if(PA_Cleavage.StateChange(delta_t, Position, Activity))
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
		if(Termination.StateChange(delta_t, Position, Activity)) return 1;
	}
	return 0;
}

void DNAstat::Init()
{
	Factor_occupancy = new int[DNALEN];
	for(int i=0;i<DNALEN;i++) Factor_occupancy[i]=0;
	for(int i=0;i<N_PRSTATUS;i++) Promoter_status[i]=0;
}

void DNAstat::CheckIC(int delta_t)
{
	// TBP bound state
	if(Promoter_status[IC]<=0)		// if IC is not bound
	{
		if(IC_on.StateChange(delta_t, 0, simTime)) Promoter_status[IC]=1;
	}
	else		// IC bound
	{
		if(IC_off.StateChange(delta_t, 0, simTime)) Promoter_status[IC]=-8;
	}
}

void DNAstat::Progress(int delta_t)
{
	CheckIC(delta_t);

	// Generate new Pol2
	int pol2_clear=1;
	for(int i=TSSPOS;i<=PAUSESITE;i++) if(Factor_occupancy[i]&POL2) pol2_clear=0;
	if(Promoter_status[IC]>0&&pol2_clear)		// IC bound and No Pol2 at promoter
		if(Recruitment.StateChange(delta_t, 0, simTime))
		{
			Pol2stat newpol2;
			newpol2.Init();
			Pol2list.push_front(newpol2);
			Factor_occupancy[TSSPOS]|=POL2;		// mark TSS occupied by pol2
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
			Factor_occupancy[it->Position]&=(~POL2);		// force terminate Pol2 at the end and erase the mark
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
		objpos[n] = TSSPOS;
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

	r = drand();
	int state = 0;
	if(delta_t>=1000) delta_t=999;

	while(state<MAX_STATE)
	{
		if(r<State_mat[delta_t][state]) break;
		state++;
	}
	
	return state;
}


// Define StateFunctionOf

void StateFunctionOf::Init(float st, char *filename)
{
	if(filename) {
		step = st;
		rateFunc.Init(filename);
		min = rateFunc.min;
		max = rateFunc.max;
		int n=(int)((max-min)/step)+1;
		stateFunc.resize(n);
		for(int i=0;i<n;i++)
		{
			float x=min+step*i;
			stateFunc[i].InitPoisson(x);
		}
	}
	else {
		min = st/5;
		max = st*2;
		step = 100;
		int n=(int)((max-min)/step)+1;
		stateFunc.resize(n);
		for(int i=0;i<n;i++)
		{
			float x=min+step*i;
			stateFunc[i].InitPoisson(x);
		}
	}
}


int StateFunctionOf::StateChange(int delta_t, float v1, float v2)
{
	int n=(rateFunc.value(v1, v2)-min)/step;
	if(n<0) n=0;
	if(n>=stateFunc.size()) n=stateFunc.size()-1;
	return stateFunc[n].StateChange(delta_t);
}


void FunctionOf::Init(char *filename)
{
	// cerr<<"Reading function from "<<filename<<endl;
	data.clear();
	v1.clear();
	v2.clear();
	min = 2e37;
	max = -2e37;
	ifstream in(filename);
	string firstline;
	getline(in,firstline);
	// cerr<<"Reading firstline: "<<firstline<<endl;
	stringstream firstss(firstline);
	string description;
	firstss>>description;
	while(true)
	{
		float v=-2e37;
		firstss>>v;
		if(v<-1e37) break;
		v1.push_back(v);
		// cerr<<v<<" saved as v1"<<endl;
	}

	while(true)
	{
		string line;
		getline(in,line);
		stringstream ss(line);
		float vv=-2e37;
		ss>>vv;
		if(vv<-1e37) break;
		v2.push_back(vv);
		vector<float> r;
		while(true)
		{
			float vvv=-2e37;
			ss>>vvv;
			if(vvv<-1e37) break;
			r.push_back(vvv);
			if(vvv<min) min = vvv;
			if(vvv>max) max = vvv;
		}
		data.push_back(r);
	}
}

float FunctionOf::value(float p1, float p2)
{
	int i1,i2,j1,j2;
	if(p1<=v1.front()) i1=i2=0;
	else if(p1>=v1.back()) i1=i2=v1.size()-1;
	else
	{
		i1=lower_bound(v1.begin(),v1.end(),p1)-v1.begin()-1;
		i2=upper_bound(v1.begin(),v1.end(),p1)-v1.begin();
		if(i2-i1>1) i1=i2=i1+1;
	}
	if(p2<=v2.front()) j1=j2=0;
	else if(p2>=v2.back()) j1=j2=v2.size()-1;
	else
	{
		j1=lower_bound(v2.begin(),v2.end(),p2)-v2.begin()-1;
		j2=upper_bound(v2.begin(),v2.end(),p2)-v2.begin();
		if(j2-j1>1) j1=j2=j1+1;
	}
	float v11=data[j1][i1],v12=data[j2][i1],v21=data[j1][i2],v22=data[j2][i2];
	
	float x21=v1[i2]-v1[i1];
	float x1=p1-v1[i1];
	float x2=v1[i2]-p1;
	float y21=v2[j2]-v2[j1];
	float y1=p2-v2[j1];
	float y2=v2[j2]-p2;
	float r1=(x21==0)?v11:x2/x21*v11+x1/x21*v21;
	float r2=(x21==0)?v12:x2/x21*v12+x1/x21*v22;
	float r=(y21==0)?r1:y2/y21*r1+y1/y21*r2;

	return r;
}

