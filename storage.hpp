#ifndef STORAGE__HPP
#define STORAGE_HPP

#include "vecteur.hpp"
#include <string>
#include <vector>

class storage //This class in an interface
{
protected:
    int ncol_;
    int nrow_;
    string name_;
public:
    //Constructor-Desrtructor
    storage(int i=0, int j=0, string name="");
    ~storage();

    //Access opreators
    virtual double operator()(int i, int j) const=0; // i=0...nrow-1 ; j=0...ncol-1
    virtual double& operator()(int i, int j) =0; // i=0...nrow-1 ; j=0...ncol-1
    pair<int, int> shape() const;

};

class sparseStorage : public storage
{
protected:
    vector<double> data_;
    vector<int> indices_;
    vector<int> indptr_;
public:
    //constructor-Destructor
    sparseStorage(int n, int m);
    ~sparseStorage();

    //Access operators
    double operator()(int i, int j) const override;
    double& operator()(int i, int j) override;
    tuple<vector<double>&, vector<int>&, vector<int>&> get_data_ref();
    tuple<vector<double>, vector<int>, vector<int>> get_data() const;
    void remove(int i, int j);
};


class denseStorage : public storage
{
protected:
    vector<vecteur> mat_;
public:
    //constructor-Destructor

    denseStorage(int n=0, int m=0, double x=0);
    ~denseStorage();

    //Access operators
    vecteur& operator[](int i);
    vecteur operator[](int i) const;
    double operator()(int i, int j) const;
    double& operator()(int i, int j);
    vector<vecteur> get_data();






};

#endif // STORAGE__HPP
