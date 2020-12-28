#ifndef _ARGUMENTS_H
#define _ARGUMENTS_H

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <utility>

namespace arg
{
    std::vector <std::string>id_string;
    std::vector <std::string>description;

    std::vector <int *> p_int;
    std::vector <char *> p_char;
    std::vector <float *> p_float;
    std::vector <double *> p_double;
    std::vector <std::string *> p_string;
    std::vector <bool *> p_bool;
    std::vector <std::vector<std::string> *> p_string_vec;
    std::vector <std::vector<int> *> p_int_vec;
    std::vector <std::vector<float> *> p_float_vec;
    std::vector <std::vector<double> *> p_double_vec;
    std::vector <char **> p_char_ptr;

    std::vector <std::pair<int,int> > arg_type;
    std::vector <bool> is_optional;
    void push_p(int *p)
    {
        arg_type.push_back(std::make_pair(1,p_int.size()));
        p_int.push_back(p);
    };

    void push_p(char *p)
    {
        arg_type.push_back(std::make_pair(2,p_char.size()));
        p_char.push_back(p);
    };

    void push_p(float *p)
    {
        arg_type.push_back(std::make_pair(3,p_float.size()));
        p_float.push_back(p);
    };

    void push_p(double *p)
    {
        arg_type.push_back(std::make_pair(4,p_double.size()));
        p_double.push_back(p);
    };

    void push_p(std::string *p)
    {
        arg_type.push_back(std::make_pair(5,p_string.size()));
        p_string.push_back(p);
    };

    void push_p(bool *p)
    {
        arg_type.push_back(std::make_pair(6,p_bool.size()));
        p_bool.push_back(p);
    };

    void push_p(std::vector<std::string> *p)
    {
        arg_type.push_back(std::make_pair(7,p_string_vec.size()));
        p_string_vec.push_back(p);
    };

    void push_p(std::vector<int> *p)
    {
        arg_type.push_back(std::make_pair(8,p_int_vec.size()));
        p_int_vec.push_back(p);
    };
    
    void push_p(std::vector<float> *p)
    {
        arg_type.push_back(std::make_pair(9,p_float_vec.size()));
        p_float_vec.push_back(p);
    };
    
    void push_p(std::vector<double> *p)
    {
        arg_type.push_back(std::make_pair(10,p_double_vec.size()));
        p_double_vec.push_back(p);
    };

    void push_p(char **p)
    {
        arg_type.push_back(std::make_pair(11,p_char_ptr.size()));
        p_char_ptr.push_back(p);
    };

    template <class Type> void push(const char *id, const char *output, Type &vp, bool isoptional=false)
    {
        std::string t;
        t=id;
        id_string.push_back(t);
        t=output;
        description.push_back(t);
        push_p(&vp); 
        is_optional.push_back(isoptional);
    };

    bool get(int argc, char *argv[])
    {
        int arit=1,n=id_string.size();
        std::vector <bool>is_passed;
        is_passed.resize(n);
        std::fill(is_passed.begin(),is_passed.end(),false);
        while(arit<argc)
        {
            for(int i=0;i<n;i++)
            {
                if(id_string[i]==argv[arit])
                {
                    is_passed[i]=true;
                    int p_it=arg_type[i].second;
                    switch(arg_type[i].first)
                    {
                        case 1:
                        {   
                            int *p=p_int[p_it];
                            *p=atoi(argv[arit+1]);
                            arit++;
                            break;
                        }
                        case 2:
                        {   
                            char *p=p_char[p_it];
                            *p=argv[arit+1][0];
                            arit++;
                            break;
                        }
                        case 3:
                        {   
                            float *p=p_float[p_it];
                            *p=atof(argv[arit+1]);
                            arit++;
                            break;
                        }
                        case 4:
                        {   
                            double *p=p_double[p_it];
                            *p=atof(argv[arit+1]);
                            arit++;
                            break;
                        }
                        case 5:
                        {   
                            std::string *p=p_string[p_it];
                            *p=argv[arit+1];
                            arit++;
                            break;
                        }
                        case 6:
                        {   
                            bool *p=p_bool[p_it];
                            *p=true;
                            break;
                        }
                        case 7:
                        {   
                            std::vector<std::string> *p=p_string_vec[p_it];
                            p->clear();
                            arit++;
                            while(arit<argc)
                            {
                                bool is_idstring=false;
                                for(int j=0;j<n;j++) if(id_string[j]==argv[arit]) is_idstring=true;
                                if(is_idstring)
                                {
                                    arit--;
                                    break;
                                }
                                else
                                {
                                    std::string tempstr=argv[arit];
                                    p->push_back(tempstr);
                                    arit++;
                                }
                            }
                            break;
                        }
                        case 8:
                        {   
                            std::vector<int> *p=p_int_vec[p_it];
                            p->clear();
                            arit++;
                            while(arit<argc)
                            {
                                bool is_idstring=false;
                                for(int j=0;j<n;j++) if(id_string[j]==argv[arit]) is_idstring=true;
                                if(is_idstring)
                                {
                                    arit--;
                                    break;
                                }
                                else
                                {
                                    int tempint=atoi(argv[arit]);
                                    p->push_back(tempint);
                                    arit++;
                                }
                            }
                            break;
                        }
                        case 9:
                        {   
                            std::vector<float> *p=p_float_vec[p_it];
                            p->clear();
                            arit++;
                            while(arit<argc)
                            {
                                bool is_idstring=false;
                                for(int j=0;j<n;j++) if(id_string[j]==argv[arit]) is_idstring=true;
                                if(is_idstring)
                                {
                                    arit--;
                                    break;
                                }
                                else
                                {
                                    float tempf=atof(argv[arit]);
                                    p->push_back(tempf);
                                    arit++;
                                }
                            }
                            break;
                        }
                        case 10:
                        {   
                            std::vector<double> *p=p_double_vec[p_it];
                            p->clear();
                            arit++;
                            while(arit<argc)
                            {
                                bool is_idstring=false;
                                for(int j=0;j<n;j++) if(id_string[j]==argv[arit]) is_idstring=true;
                                if(is_idstring)
                                {
                                    arit--;
                                    break;
                                }
                                else
                                {
                                    double tempf=atof(argv[arit]);
                                    p->push_back(tempf);
                                    arit++;
                                }
                            }
                            break;
                        }
                        case 11:
                        {   
                            char **p=p_char_ptr[p_it];
                            *p=argv[arit+1];
                            arit++;
                            break;
                        }
                    }
                    break;
                }
            }
            arit++;
        }
        for(int i=0;i<n;i++)
        {
            if(is_optional[i]==false&&is_passed[i]==false)
            {
                std::cout<<"Arguments"<<std::endl;
                for(int j=0;j<n;j++) std::cout<<id_string[j]<<"\t:\t"<<description[j]<<std::endl;
                return false;
            }
        }
        return true;
    };
};
#endif
