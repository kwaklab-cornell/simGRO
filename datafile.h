#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <map>

using namespace std;

struct datafile
{
    int n_skipped_line;
	 
    std::vector <std::string> skipped_line;
    std::vector <std::vector<std::string> > data;

    void load(char *fn, int skip=1)
    {
        n_skipped_line=skip;
        std::ifstream in(fn);
        std::string tempstr;

        // read skipped line
        for(int i=0;i<n_skipped_line;i++)
        {
            std::getline(in,tempstr);
            skipped_line.push_back(tempstr);
        }

        // load data
        while(true)
        {
			std::string templine;
			std::getline(in,templine);
			if(in.fail()) break;
			std::stringstream ss(templine);
			std::vector<std::string> line;	
            while(true)
            {
				std::string buf;
				ss>>buf;
				if(ss.fail()) break;
				line.push_back(buf);
            }
			data.push_back(line);
        }
        in.close();
    };

    void get_data(std::vector<float> &v, int column)
    {
		int line=0;
		v.clear();
		while(true)
		{
			if(line>=data.size()) break;
			if(column>=data[line].size()) break;
			std::stringstream ss(data[line][column]);
			float a;
			ss>>a;
			if(ss.fail()) break;
			v.push_back(a);
			line++;
		}
	};

    void get_data(std::vector<int> &v, int column)
    {
		int line=0;
		v.clear();
		while(true)
		{
			if(line>=data.size()) break;
			if(column>=data[line].size()) break;
			std::stringstream ss(data[line][column]);
			int a;
			ss>>a;
			if(ss.fail()) break;
			v.push_back(a);
			line++;
		}
	};

    void get_data(std::vector<std::string> &v, int column)
    {
		int line=0;
		v.clear();
		while(true)
		{
			if(line>=data.size()) break;
			if(column>=data[line].size()) break;
			v.push_back(data[line][column]);
			line++;
		}
    };

    template <class T>
    void set_map(std::map<std::string,T> &m, int key_column, int data_column)
    {
        std::vector<T> d;
        get_data(d,data_column);
        int size=data[key_column].size();
        for(int i=0;i<size;i++) m[data[key_column][i]]=d[i];
    };
};


