#include <cmath>

using namespace std;

struct kernel
{
    kernel(float bandwidth);
    void addto(float value, int position, float *destination, int size);

    float bw;
    int n;
    float *val;
};

kernel::kernel(float bandwidth)
{
    bw=bandwidth;
    n=2*(int)(5*bw+0.5)+1;
    val=new float[n];

    float sum=0;
    for(int i=0;i<n;i++)
    {
        float x=(float)i/bw-5;
        val[i]=exp(-x*x/2);
        sum+=val[i];
    }
    for(int i=0;i<n;i++) val[i]/=sum;
}

void kernel::addto(float value, int position, float *destination, int size)
{
    for(int i=0;i<n;i++)
    {
        int pos=position+i-5*bw;
        if(pos>=0&&pos<size) destination[pos]+=value*val[i];
    }
}


