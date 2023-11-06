#include "storage.hpp"
#include "vecteur.hpp"
#include<string>
#include <tuple>
#include<vector>


///Class Storage
storage::storage(int n_line, int n_col, string name):ncol_(n_col), nrow_(n_line), name_(name){}
storage::~storage() { ncol_=0; nrow_=0; name_ ="";}

pair<int, int> storage::shape() const{return std::make_pair(ncol_, nrow_);}

///Class sparseStorage
//Constructor-Destructor
sparseStorage::sparseStorage(int n_line, int n_col): storage(n_line, n_col, "Sparse"), data_(), indices_(), indptr_(n_col+1,0){}
sparseStorage::~sparseStorage(){ncol_=0; nrow_=0; name_ =""; data_.clear(); indices_.clear(); indptr_.clear();}

//Access Opreators
double sparseStorage::operator()(int line_index, int col_index) const
{
    if (data_.size()==0) return (0.0);
    for(int k=indptr_[line_index]; k< indptr_[line_index+1]; k++){
        if ((indices_[k]==col_index)&&  (k<(int)indices_.size())) return(data_[k]);
    }
    return (0.0);
}

double& sparseStorage::operator()(int line_index, int col_index)
{
    if (data_.size()==0)
    {
        data_.push_back(0);
        indices_.push_back(col_index);
        for (int k=line_index+1; k<nrow_+1; k++) indptr_[k]++;
        return(data_.at(0));
    }

    int idx = indptr_[line_index];
    while ((idx < indptr_[line_index+1]) && (col_index >indices_[idx])) idx++;
    if ((indices_[idx]==col_index)&&  (idx<(int)indices_.size()) && (idx < indptr_[line_index+1])) return(data_[idx]);

    data_.insert(data_.begin()+(idx),0.0);
    indices_.insert(indices_.begin() +(idx),col_index);
    for (int k=line_index+1; k<nrow_+1; k++) indptr_[k]++;
    return(data_[idx]);
}
void sparseStorage::remove(int line_index, int col_index)
{
    int idx = indptr_[line_index];
    while ((idx < indptr_[line_index+1]) && (col_index >indices_[idx])) idx++;
    if ((indices_[idx]==col_index)&&  (idx<(int)indices_.size()) ){
        data_.erase(data_.begin()+idx);
        indices_.erase(indices_.begin()+idx);
        for (int k=line_index+1; k<nrow_+1; k++) indptr_[k]--;
    }
    else{
        cerr<< "element to remove does not exist";
        exit(-1);
    }
    return;
}

tuple<vector<double>&, vector<int>&, vector<int>&> sparseStorage::get_data_ref()
{
    return std::tie(data_, indptr_, indices_);
}

tuple<vector<double>, vector<int>, vector<int>> sparseStorage::get_data() const
{
    return std::make_tuple(data_, indptr_, indices_);
}

///Class denseStorage
denseStorage::denseStorage(int n_line, int n_col, double x):storage(n_line, n_col, "Dense"), mat_(vector<vecteur>(n_line,vecteur(n_col,x))){}
denseStorage::~denseStorage(){mat_.clear();}

vector<vecteur> denseStorage::get_data(){return mat_;}

//Access Opreators
vecteur& denseStorage::operator[](int i) {return(mat_.at(i));}
vecteur denseStorage::operator[](int i) const {return(mat_.at(i));}
double denseStorage::operator()(int i, int j) const {return ( (*this)[i][j]);}
double& denseStorage::operator()(int i, int j) {return ( (*this)[i][j]);}
