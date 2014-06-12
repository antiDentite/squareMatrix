#include "LinAlg.h"
#include <math.h>

///vct Class

template <class T,unsigned int dim> vct<T,dim>::~vct(){};
template <class T,unsigned int dim> vct<T,dim>::vct(const vct<T,dim> &v):std::vector<T>(v){};
template <class T,unsigned int dim> vct<T,dim>::vct():std::vector<T>(dim){};
template <class T,unsigned int dim> vct<T,dim>::vct(const T & x):std::vector<T>(dim,x){};
template <class T, unsigned int dim> T & vct<T,dim>::operator()(const unsigned int & pos){return (*this)[pos];};
template <class T, unsigned int dim> vct<T,dim> & vct<T,dim>:: operator +=(const vct<T,dim> & v){

std::transform(this->begin(),this->end(),v.begin(),this->begin(),std::plus<T>());
return *(this);
};

template <class T,unsigned int dim> vct<T,dim> & vct<T,dim>:: operator -=(const vct<T,dim> & v){

(*this)+=-v;

return *(this);

};

template <class T,unsigned int dim> vct<T,dim> vct<T,dim>:: operator+(const vct<T,dim> & v_1)const{vct v_2(*(this));v_2+=v_1;return v_2;};
template <class T,unsigned int dim> vct<T,dim> vct<T,dim>:: operator-(const vct<T,dim> &v_1)const{vct  v_2(*(this));v_2-=v_1;return v_2;};

template <class T,unsigned int dim> T vct<T,dim>::operator *(const vct<T,dim> & v_1)const{

std::vector<T> v(dim);
std::transform(v_1.begin(),v_1.end(),this->begin(),v.begin(),std::multiplies<T>());
T s(0);

for(int i=0;i<dim;i++)s+=v[i];

return s;

};


template <class T, unsigned int dim> vct<T,dim> vct<T,dim>::Random(const unsigned int & s,
const double & a, const double & b){
vct<T,dim> v;

///change the seed value

srand(s);

///load up the vector

for(int i=0;i<dim;i++)v[i]=rand()*((b-a)/RAND_MAX) + a;
return v;
};



///***************

///SQUARE class

template <class T,class S> const S & SQUARE<T,S>:: return_to_derived_2() const{

return static_cast<const S &>(*this);

};

template <class T,class S> T SQUARE<T,S>:: entry(unsigned int i,unsigned int j) const{
return static_cast<const S &>(*this).entry(i,j);};

template <class T,class S>  template<class Z> std::vector< std::vector<T> > SQUARE<T,S>:: slowAdd(const SQUARE<T,Z> & M_1)const{
unsigned int dim=M_1.return_to_derived_2().M.size();
std::vector< std::vector<T> > M_2(dim);

for(int i=0;i<dim; i++){

for(int j=0;j<dim;j++)M_2[i].push_back(entry(i,j)+
(M_1.return_to_derived_2()).entry(i,j));

};

return M_2;
};

template <class T,class S>   template<class Z> std::vector< std::vector<T> > SQUARE<T,S>:: slowMinus(const SQUARE<T,Z> & M_1)const{
unsigned int dim=M_1.return_to_derived_2().M.size();
std::vector< std::vector<T> > M_2(dim);

for(int i=0;i<dim; i++){

for(int j=0;j<dim;j++)M_2[i].push_back(entry(i,j)-M_1.return_to_derived_2().entry(i,j));

};


return M_2;

};

template <class T,class S>   template<class Z> std::vector< std::vector<T> > SQUARE<T,S>:: slowMult(const SQUARE<T,Z> & M_1)const{
unsigned int dim=M_1.return_to_derived_2().M.size();
std::vector< std::vector<T> > M_2(dim);

for(int i=0;i<dim; i++){

for(int j=0;j<dim;j++){
M_2[i].push_back(0);
for(int k=0;k<dim;k++)M_2[i][j]+=entry(i,k)*M_1.return_to_derived_2().entry(k,j);

};

};

return M_2;
};

template<class T, class S> std::vector< std::vector<T> > SQUARE<T,S>::slowTranspose()const{

unsigned int dim=return_to_derived_2().M.size();

std::vector< std::vector<T> > M_2(dim);

///load up the container

for(int i=0;i<dim; i++){

for(int j=0;j<dim;j++)M_2[i].push_back(entry(j,i));

};

return M_2;

};

template<class T, class S> std::vector< std::vector<T> > SQUARE<T,S>::negation()const{
unsigned int dim=return_to_derived_2().M.size();
std::vector< std::vector<T> > M_2(dim);

for(int i=0;i<dim;i++){
M_2[i].resize(return_to_derived_2().M[i].size());
std::transform(return_to_derived_2().M[i].begin(),return_to_derived_2().M[i].end(),M_2[i].begin(),std::negate<T>());

};

return M_2;
};

template<class T, class S> std::vector< std::vector<T> > SQUARE<T,S>::scalarMult(const T & t)const{
std::vector< std::vector<T> > M_2;
M_2=return_to_derived_2().M;
std::vector<T> v_1(M_2.size(),t);

for(int i=0;i<M_2.size();i++){

std::transform(M_2[i].begin(),M_2[i].end(),v_1.begin(),M_2[i].begin(),std::multiplies<T>());

};

return M_2;

};

///***************

///FULL<T> class

template <class T,unsigned int dim> FULL<T,dim>::FULL(){

M.resize(dim);

for(int i=0;i<dim;i++){M[i].resize(dim);};

};

template <class T,unsigned int dim> FULL<T,dim>::~FULL(){}; ///do nothing

template <class T,unsigned int dim> T & FULL<T,dim>::operator()(const unsigned int & i, const unsigned int & j){

assert(i<dim && j<dim);

return M[i][j];

};

template <class T,unsigned int dim> T  FULL<T,dim>::entry(const unsigned int & i, const unsigned int & j)const{

assert(i<dim && j<dim);

return M[i][j];

};

template <class T,unsigned int dim> FULL<T,dim> & FULL<T,dim>:: operator =(const BAND<T,dim> & M_1){
int adj;
for(int i=0;i<dim;i++){

adj =std::max<int>(0,i-M_1.p);
M[i].assign(dim,0);

for(int j=adj;j<=std::min<int>(dim-1,i+M_1.q);j++)M[i][j]=M_1.M[i][j-adj];


};

return *this;
};

template <class T,unsigned int dim> FULL<T,dim> & FULL<T,dim>:: operator =(const TRIDAG<T,dim> & M_1){

int adj;

for(int i=0;i<dim;i++){

adj =std::max<int>(0,i-1);
M[i].assign(dim,0);

for(int j=adj;j<=std::min<int>(dim-1,i+1);j++)M[i][j]=M_1.M[i][j-adj];

};

return *this;

};

///addition
template<class T,unsigned int dim> FULL<T,dim>  FULL<T,dim>:: operator + ( const FULL<T,dim> & M_1 ) const {
FULL<T,dim> M_2;
for(int i=0;i<dim;i++){

std::transform(M[i].begin(),M[i].end(),M_1.M[i].begin(),M_2.M[i].begin(),std::plus<T>());

};

return M_2;

};

template <class T,unsigned int dim> FULL<T,dim>  FULL<T,dim>:: operator + (const BAND<T,dim> & M_2) const{

return M_2 + *this;};
template <class T,unsigned int dim> FULL<T,dim>  FULL<T,dim>:: operator + (const TRIDAG<T,dim> & M_2)const {

return M_2+*this;

};

///substraction

template <class T,unsigned int dim> FULL<T,dim>  FULL<T,dim>:: operator - (const FULL<T,dim> & M_1 )const{

FULL<T,dim> M_2;
for(int i=0;i<dim;i++){

std::transform(M[i].begin(),M[i].end(),M_1.M[i].begin(),M_2.M[i].begin(),std::minus<T>());

};

return M_2;

};
template <class T, unsigned int dim> FULL<T,dim>  FULL<T,dim>:: operator - (const BAND<T,dim> & M_2)const{

BAND<T,dim> M_3(-M_2);

return M_3+*this;};

template <class T, unsigned int dim> FULL<T,dim>  FULL<T,dim>:: operator - (const TRIDAG<T,dim> & M_2)const{
TRIDAG<T,dim> M_3(-M_2);

return M_3+*this;

};

///multiplication
template<class T, unsigned int dim> vct<T, dim>  FULL<T,dim>:: operator * (const vct<T,dim> & v)const{

vct<T,dim> v_1(0);
std::vector< std::vector<T> > v_2(dim);

for(int i=0;i<dim;i++){
v_2[i].resize(dim);
std::transform(M[i].begin(),M[i].end(),v.begin(),v_2[i].begin(),std::multiplies<T>());

for(int k=0;k<dim;k++)v_1[i]+=v_2[i][k];
};

return v_1;

};

template<class T, unsigned int dim>  FULL<T,dim>   FULL<T,dim>::operator *(const FULL<T,dim> & M_1 )const{

FULL<T,dim> M_2;

for(int i=0;i<dim;i++){

for(int j=0;j<dim;j++){
M_2.M[i][j]=0;
for(int k=0;k<dim;k++){

M_2.M[i][j]+=M[i][k]*M_1.M[k][j];};

};

};

return M_2;

};

template<class T, unsigned int dim> FULL<T,dim>  FULL<T,dim>:: operator * (const BAND<T, dim> & M_1)const{

FULL<T,dim> M_2(transpose(M_1)*transpose(*this));

return transpose(M_2);

};
template<class T, unsigned int dim> FULL<T,dim>  FULL<T,dim>:: operator * (const TRIDAG<T, dim> & M_1)const{

FULL<T,dim> M_2(transpose(M_1)*transpose(*this));

return transpose(M_2);

};

///solver

template<class T, unsigned int dim> vct<T,dim>  FULL<T,dim>::fastSolve(const vct<T, dim> & b)const{

std::vector<T> v(dim);int k;

FULL<T,dim> M_1(*this);

std::vector< std::vector<T>  > L(dim),U(dim);

T s(0), s_1(0);

for(int i=0;i<dim;i++){
L[i].resize(i+1);
L[i][i]=1;
U[i].resize(dim-i);};

for(int i=0;i<dim;i++){

for(int j=i;j<dim;j++){

s=M_1.M[i][j];

for(int k=0;k<i;k++){s-=L[i][k]*U[k][j-k];}
U[i][j-i]=s;

};

for(int j=i+1;j<dim;j++){

s_1=M_1.M[j][i];

for(int k=0;k<i;k++){s_1-=L[j][k]*U[k][i-k];}

L[j][i]=s_1/U[i][0];

};};



std::vector<T> y_1(dim); vct<T,dim> x;

s=L[1][0]*b[0];y_1[0]=b[0];

for(int i=1;i<M.size()-1;i++){

for(int j=1;j<i;j++){

s+=L[i][j]*y_1[j];
}

y_1[i]=(b[i]-s);

s=L[i+1][0]*b[0];

};


for(int j=1;j<M.size()-1;j++){

T b_10=L[M.size()-1][j]*y_1[j];
s+=b_10;
};

y_1[M.size()-1]=b[M.size()-1]-s;

s=0;


x[M.size()-1]=y_1[M.size()-1]/U[M.size()-1][0];

for(int i=M.size()-2;i>=0;i--){

for(int j=i+1;j<M.size();j++){

s+=U[i][j-i]*x[j];

};

x[i]=(y_1[i]-s)/U[i][0];

s=0;
};

return x;

};

///Matrix Generators

template<class T, unsigned int dim> FULL<T,dim>  FULL<T,dim>::Random(unsigned int s,double a, double b){
FULL<T,dim> M_1;
///change the seed value

srand(s);

///load up the containers:

for(int i=0;i<dim;i++){

for(int k=0;k<dim;k++){

M_1.M[i][k]=rand()*((b-a)/RAND_MAX) + a;

};

};

return M_1;
};


///****************

///TRIDAG Class

template<class T,unsigned int dim> TRIDAG<T,dim>::TRIDAG():M(dim){

M[0].resize(2);

for(int i=1;i<dim-1;i++)M[i].resize(3);

M[dim-1].resize(2);


};
template<class T,unsigned int dim>  T & TRIDAG<T,dim>::operator()(const unsigned int & i, const unsigned int & j){
assert(i<dim && j<dim);
T t(0); T & r=t;
if(i>1+j || j>1+i){return r;} else{

if(i==0){if(j==0){return M[0][0];} else{return M[0][1];};}

else if(i==dim-1){if(j==dim-2){return M[i][0];} else{return M[dim-1][1];}}

else{if(j==i-1){return M[i][0];} else if(j==i){return M[i][1];}
else{return M[i][2];}

};

};};

template<class T,unsigned int dim>  T  TRIDAG<T,dim>::entry(const unsigned int & i, const unsigned int & j)const{
assert(i<dim && j<dim);
T t(0); T & r=t;
if(i>1+j || j>1+i){return r;} else{

if(i==0){if(j==0){return M[0][0];} else{return M[0][1];};}

else if(i==dim-1){if(j==dim-2){return M[i][0];} else{return M[dim-1][1];}}

else{if(j==i-1){return M[i][0];} else if(j==i){return M[i][1];}
else{return M[i][2];}

};

};};

template <class T,unsigned int dim> vct<T,dim>  TRIDAG<T,dim>::operator *(const vct<T,dim> & v)const{

vct<T,dim> v_1(dim);

v_1[0]=M[0][0]*v[0]+M[0][1]*v[1];

for(int i=1;i<dim-1;i++){

v_1[i]=M[i][0]*v[i-1] + M[i][1]*v[i] + M[i][2]*v[i+1];

};

v_1[dim-1]=M[dim-1][0]*v[dim-2]+M[dim-1][1]*v[dim-1];

vct<T,dim> & R=v_1;
return R;
};

template <class T, unsigned int dim> FULL<T,dim>  TRIDAG<T,dim>::operator *(const FULL<T,dim> & M_1)const{

FULL<T,dim> M_2;
int adj;
for(int i=0;i<dim;i++){
adj=std::max<int>(0,i-1);
for(int j=0;j<dim;j++){
M_2.M[i][j]=0;

for(int k=0;k<M[i].size();k++){
M_2.M[i][j]+=M[i][k]*M_1.M[k+adj][j];};

};

};
return M_2;
};

template<class T,unsigned int dim> BAND<T,dim>  TRIDAG<T,dim>::operator *(const TRIDAG<T,dim> & M_1)const{

BAND<T,dim> M_2(2,2);
int adj,adj_2,adj_3;
for(int i=0;i<=dim-1;i++){

adj=std::max<int>(0,i-2);
adj_2=std::max<int>(0,i-1);

for(int j=std::max<int>(0,i-2);j<=std::min<int>(i+2,dim-1);j++){
M_2.M[i][j-adj]=0;

for(int k=std::max<int>(0,j-1-adj_2);k<=std::min<int>(M[i].size()-1,j+1-adj_2);k++)
{adj_3=std::max<int>(0,k+adj_2-1);
M_2.M[i][j-adj]+=M[i][k]*M_1.M[k+adj_2][j-adj_3];};

};

};

return M_2;

};

template<class T,unsigned int dim> BAND<T,dim>  TRIDAG<T,dim>::operator *(const BAND<T,dim> & M_1)const{

int p_2=std::min<int>(dim-1,1+M_1.p);
int q_2=std::min<int>(dim-1,1+M_1.q);

BAND<T,dim> M_2(p_2,q_2);

int adj,adj_2,adj_3;
for(int i=0;i<dim;i++){

adj=std::max<int>(0,i-p_2);
adj_2=std::max<int>(0,i-1);

for(int j=std::max<int>(i-p_2,0);j<=std::min<int>(dim-1,i+q_2);j++){

M_2.M[i][j-adj]=0;


for(int k=std::max<int>(std::max<int>(0,i-1),j-M_1.q);k<=std::min<int>(std::min<int>(dim-1,i+1),j+M_1.p);k++){
adj_3=std::max<int>(0,k-M_1.p);
M_2.M[i][j-adj]+=M[i][k-adj_2]*M_1.M[k][j-adj_3];

};

};


};



return M_2;
};

template<class T,unsigned int dim> FULL<T,dim>  TRIDAG<T,dim>::operator + (const FULL<T,dim> & M_1)const{

FULL<T,dim> M_2;

M_2.M=M_1.M;
int adj;
for(int i=0;i<dim;i++){
adj=std::max<int>(0,i-1);

for(int k=0;k<M[i].size();k++){
M_2.M[i][k+adj]+=M[i][k];
};

};

return M_2;

};

template<class T, unsigned int dim> BAND<T,dim> TRIDAG<T,dim>::operator + (const BAND<T,dim> & M_1)const{

int p_2=std::min<int>(1+M_1.p,dim-1);
int q_2=std::min<int>(1+M_1.q,dim-1);
BAND<T,dim> M_2(p_2,q_2);

int adj,adj_2;
for(int i=0;i<dim;i++){
M_2.M[i].assign(0,M_2.M[i].size());
adj=std::max<int>(i-1,0);
adj_2=std::max<int>(i-p_2,0);

for(int k=0;k<M[i].size();k++){

M_2.M[i][k+adj-adj_2]=M[i][k];

};

adj=std::max<int>(i-M_1.p,0);
adj_2=std::max<int>(i-p_2,0);
for(int k=0;k<M_1.M[i].size();k++){

M_2.M[i][k+adj-adj_2]+=M_1.M[i][k];

};

};

return M_2;

};

template<class T,unsigned int dim> TRIDAG<T,dim>  TRIDAG<T,dim>::operator + (const TRIDAG<T,dim> & M_1)const{

TRIDAG<T,dim> M_2;

for(int i=0;i<dim;i++)transform(M[i].begin(),M[i].end(),M_1.M[i].begin(),M_2.M[i].begin(),std::plus<T>());
return M_2;
};

template<class T, unsigned int dim> FULL<T,dim>  TRIDAG<T,dim>::operator - ( const FULL<T,dim> &  M_1)const{

FULL<T,dim> M_2(-M_1);

return *this + M_2;

};
template<class T,unsigned int dim> BAND <T,dim> TRIDAG<T,dim>::operator - (const BAND<T,dim> & M_1)const{

BAND<T,dim> M_2(-M_1);

return *this + (M_2);

};
template<class T,unsigned int dim> TRIDAG<T,dim>  TRIDAG<T,dim>::operator - (const TRIDAG<T,dim> & M_1)const{
TRIDAG<T,dim> M_2(-M_1);
return *this + M_2;};

///Fast Solver

template <class T, unsigned int dim> vct<T,dim> TRIDAG<T,dim>:: fastSolve(const vct<T,dim> & b)const{

if(dim==2){
vct<T,dim> x;

x[0]=(b[1]-M[1][1]*b[0]/M[0][1])/(M[1][0]-M[0][0]*M[1][1]/M[0][1]);

x[1]=(b[0]-M[0][0]*x[0])/M[0][1];

return x;
}

else{
std::vector<T> v(dim);int k;

std::vector< std::vector<T> > L(dim),U(dim);

U[0].assign(1,M[0][1]);U[0].insert(U[0].begin(),M[0][0]);

U[1].assign(1,M[1][2]);

U[1].insert(U[1].begin(),M[1][1]-(M[1][0]/U[0][0])*M[0][1]);

for(int i=2;i<dim-1;i++){

U[i].assign(1,M[i][2]);

U[i].insert(U[i].begin(),M[i][1]-(M[i][0]/U[i-1][0])*M[i-1][2]);

}

U[dim-1].insert(U[dim-1].begin(),M[dim-1][1]-(M[dim-1][0]/U[dim-2][0])*M[dim-2][2]);



L[0].assign(1,1);


for(int i=1;i<dim;i++){

L[i].assign(1,L[0][0]);
L[i].insert(L[i].begin(),M[i][0]/U[i-1][0]);

}

std::vector<T> y_1(dim); vct<T,dim> x(dim);

y_1[0]=b[0];

for(int i=1;i<dim;i++){

y_1[i]=b[i]-(L[i][0]*y_1[i-1]);

}

x[dim-1]=y_1[dim-1]/U[dim-1][0];

for(int i=dim-2;i>=0;i--){

x[i]=(y_1[i]-(U[i][1]*x[i+1]))/U[i][0];

}


return x;};

};


///Matrix Generator

template<class T, unsigned int dim> TRIDAG<T,dim>  TRIDAG<T,dim>::Random(unsigned int s,double a, double b){

TRIDAG<T,dim> M_1;

///change the seed value

srand(s);

///load up the containers:

for(int i=0;i<dim;i++){

for(int k=0;k<M_1.M[i].size();k++){

M_1.M[i][k]=rand()*((b-a)/RAND_MAX) + a;

};

};


return M_1;
};


///*****************

///BAND Class

template <class T,unsigned int dim> BAND<T,dim>::BAND(const unsigned int & p_1,const unsigned int & q_1){

M.resize(dim);
p=std::min<int>(dim-1,p_1);
q=std::min<int>(dim-1,q_1);

for(int i=0;i<dim;i++){
M[i].resize(std::min<int>(i+1,p+1)+std::min<int>(dim-i-1,q));
};


};

template <class T, unsigned int dim> BAND<T,dim>::BAND():p(0),q(0){
M.resize(dim);
for(int i=0;i<dim;i++)M[i].resize(1);
};

template<class T, unsigned int dim> BAND<T,dim>::~BAND(){};

template<class T,unsigned int dim> T & BAND<T,dim>::operator()(const unsigned int & i, const unsigned int & j){

assert(i<dim && j<dim);
T t(0); T & r=t;
if(i>p+j || j>q+i){return r;} else{

if(i>p){return M[i][p+(j-i)];}
else{
return M[i][j];};

};

};

template<class T,unsigned int dim> T BAND<T,dim>::entry(const unsigned int & i, const unsigned int & j)const{

assert(i<dim && j<dim);
T t(0); T & r=t;
if(i>p+j || j>q+i){return r;} else{

if(i>p){return M[i][p+(j-i)];}
else{
return M[i][j];};

};

};

template<class T,unsigned int dim> BAND<T,dim> & BAND<T,dim>::operator =(const TRIDAG<T,dim> & M_1){

M.clear();

M=M_1.M;

p=1;q=1;

return *this;

};

///addition

template<class T,unsigned int dim>  BAND<T,dim>  BAND<T,dim>::operator +(const BAND<T,dim> & M_1)const{

int p_2=std::max<int>(p,M_1.p);
int q_2=std::max<int>(q,M_1.q);

BAND<T,dim> M_2(p_2,q_2);


if(p<p_2){
for(int i=0;i<dim;i++){

for(int j=std::max<int>(0,i-p_2);j<std::max<int>(0,i-p);j++)M_2.M[i][j-(std::abs<int>(i-p_2)+(i-p_2))/2]=M_1.M[i][j-(std::abs<int>(i-M_1.p)+(i-M_1.p))/2];

};

if(q<q_2){
for(int i=0;i<dim;i++){

for(int j=std::max<int>(0,i-p);j<std::min<int>(dim,i+q);j++)M_2.M[i][j-(std::abs<int>(i-p_2)+(i-p_2))/2]=M_1.M[i][j-(std::abs<int>(i-M_1.p)+(i-M_1.p))/2] + M[i][j-(std::abs<int>(i-p)+(i-p))/2];

};

for(int i=0;i<M.size();i++){

for(int j=std::min<int>(M.size(),i+q+1);j<=std::min<int>(M.size()-1,i+q_2);j++)M_2.M[i][j-(std::abs<int>(i-p_2)+(i-p_2))/2]=M_1.M[i][j-(std::abs<int>(i-M_1.p)+(i-M_1.p))/2];

};

return M_2;

} else{

for(int i=0;i<dim;i++){

for(int j=std::max<int>(0,i-p);j<=std::min<int>(dim-1,i+M_1.q);j++)M_2.M[i][j-(std::abs<int>(i-p_2)+(i-p_2))/2]=M_1.M[i][j-(std::abs<int>(i-M_1.p)+(i-M_1.p))/2] + M[i][j-(std::abs<int>(i-p)+(i-p))/2];

};

for(int i=0;i<dim;i++){

for(int j=std::min<int>(dim,i+M_1.q+1);j<=std::min<int>(dim-1,i+q_2);j++)M_2.M[i][j-(std::abs<int>(i-p_2)+(i-p_2))/2]=M[i][j-(std::abs<int>(i-p)+(i-p))/2];

};

return M_2;

};


}

else{

for(int i=0;i<dim;i++){

for(int j=std::max<int>(0,i-p_2);j<std::max<int>(0,i-M_1.p);j++)M_2.M[i][j-(std::abs<int>(i-p_2)+(i-p_2))/2]=M[i][j-(std::abs<int>(i-p)+(i-p))/2];

};

if(q<q_2){
for(int i=0;i<dim;i++){

for(int j=std::max<int>(0,i-M_1.p);j<std::min<int>(dim,i+q);j++)M_2.M[i][j-(std::abs<int>(i-p_2)+(i-p_2))/2]=M_1.M[i][j-(std::abs<int>(i-M_1.p)+(i-M_1.p))/2] + M[i][j-(std::abs<int>(i-p)+(i-p))/2];

};

for(int i=0;i<dim;i++){

for(int j=std::min<int>(dim,i+q+1);j<=std::min<int>(dim-1,i+q_2);j++)M_2.M[i][j-(std::abs<int>(i-p_2)+(i-p_2))/2]=M_1.M[i][j-(std::abs<int>(i-M_1.p)+(i-M_1.p))/2];

};

return M_2;

} else{

for(int i=0;i<dim;i++){

for(int j=std::max<int>(0,i-M_1.p);j<=std::min<int>(dim-1,i+M_1.q);j++)M_2.M[i][j-(std::abs<int>(i-p_2)+(i-p_2))/2]=M_1.M[i][j-(std::abs<int>(i-M_1.p)+(i-M_1.p))/2] + M[i][j-(std::abs<int>(i-p)+(i-p))/2];

};

for(int i=0;i<dim;i++){

for(int j=std::min<int>(dim,i+M_1.q+1);j<=std::min<int>(dim-1,i+q_2);j++)M_2.M[i][M_2.M[i][j-(std::abs<int>(i-p_2)+(i-p_2))/2]]=M[i][j-(std::abs<int>(i-p)+(i-p))/2];

};

return M_2;

};


};


return M_2;

};

template<class T, unsigned int dim> FULL<T,dim>   BAND<T,dim>::operator +(const FULL<T,dim> & M_1)const{

FULL<T,dim> M_2;
M_2.M=M_1.M;
int adj;
for(int i=0;i<dim;i++){
adj=std::max<int>(0,i-p);
for(int k=0;k<M[i].size();k++){

M_2.M[i][k+adj]+=M[i][k];};

};


return M_2;

};

template<class T,unsigned int dim> BAND<T,dim>  BAND<T,dim>:: operator +(const TRIDAG<T,dim> & M_1)const{
return M_1+*this;
};

///Substraction

template<class T, unsigned int dim> FULL<T,dim>  BAND<T,dim>::operator -(const FULL<T,dim> & M_1)const{

FULL<T,dim> M_2(-M_1);
return *this + M_2;};

template<class T, unsigned int dim> BAND <T,dim>  BAND<T,dim>::operator -(const BAND<T,dim> & M_1)const{

BAND<T,dim> M_2(-M_1);
return *this + M_2;

};

template<class T,unsigned int dim> BAND <T,dim>  BAND<T,dim>::operator -(const TRIDAG<T,dim> & M_1)const{

TRIDAG<T,dim> M_2(-M_1);
return *this + M_2;};

///multiplication

template <class T,unsigned int dim> vct<T,dim>  BAND<T,dim>:: operator *(const vct<T,dim> & v)const{
vct<T,dim> v_1;

for(int i=0;i<dim;i++){
v_1[i]=0;
for(int k=0;k<M[i].size();k++){

v_1[i]+=M[i][k]*v[k+std::max<int>(0,i-p)];

};

};


return v_1;
};

template<class T,unsigned int dim> FULL<T,dim>   BAND<T,dim> ::operator *(const FULL<T,dim> & M_1)const{

FULL<T,dim> M_2;
int adj;
for(int i=0;i<dim;i++){

adj=std::max<int>(0,i-p);

for(int j=0;j<dim;j++){

M_2.M[i][j]=0;
for(int k=0;k<M[i].size();k++){

M_2.M[i][j]+=M[i][k]*M_1.M[k+adj][j];

};


};

};
return M_2;
};

template<class T,unsigned int dim> BAND<T,dim>  BAND<T,dim> ::operator *(const BAND<T,dim> & M_1)const{

int p_2=std::min<int>(dim-1,p+M_1.p);
int q_2=std::min<int>(dim-1,q+M_1.q);

BAND<T,dim> M_2(p_2,q_2);

int adj,adj_2,adj_3;
for(int i=0;i<dim;i++){

adj=std::max<int>(0,i-p_2);
adj_2=std::max<int>(0,i-p);

for(int j=std::max<int>(i-p_2,0);j<=std::min<int>(dim-1,i+q_2);j++){

M_2.M[i][j-adj]=0;

for(int k=std::max<int>(std::max<int>(0,i-p),j-M_1.q);k<=std::min<int>(std::min<int>(dim-1,i+q),j+M_1.p);k++){
adj_3=std::max<int>(0,k-M_1.p);
M_2.M[i][j-adj]+=M[i][k-adj_2]*M_1.M[k][j-adj_3];

};

};


};

return M_2;

};

template<class T,unsigned int dim> BAND<T,dim>  BAND<T,dim>:: operator *(const TRIDAG<T,dim> & M_1)const{

BAND<T,dim> M_2(transpose(M_1)*transpose(*this));

return transpose(M_2);

};
///Matrix Generator

template<class T, unsigned int dim> BAND<T,dim>  BAND<T,dim>::Random(unsigned int p, unsigned int q,unsigned int s,double a, double b){

BAND<T,dim> M_1(p,q);

///change the seed value

srand(s);

///load up the containers:

for(int i=0;i<dim;i++){

for(int k=0;k<M_1.M[i].size();k++){

M_1.M[i][k]=rand()*((b-a)/RAND_MAX) + a;

};

};


return M_1;
};

template<class T, unsigned int dim> vct<T,dim> BAND<T,dim>::fastSolve(const vct<T,dim> & b)const{

std::vector<T> v; int k;

std::vector< std::vector<T> > L(dim),U(dim);

for(int i=0;i<dim;i++){
if(p<=i){L[i].resize(p+1);L[i][p]=1;} else{L[i].resize(i+1);L[i][i]=1;};
if(q<=(dim-(i+1))){U[i].resize(q+1);} else{U[i].resize(dim-i);};
};

T t;

///build the first row of U
for(int k=0;k<U[0].size();k++){
U[0][k]=M[0][k];
};

///build the first column of L

for(int k=1;k<=std::min<int>(dim-1,p);k++){
L[k][0]=M[k][0]/U[0][0];
};
int adj;

for(int i=1;i<dim;i++){

adj=std::max<int>(i-p,0);

///the next row of U

for(int j=0;j<U[i].size();j++){

t=M[i][(j+i)-adj];

for(int k=0;k<L[i].size()-1;k++){

t-=L[i][k]*U[k+adj][i+j-(k+adj)];};

U[i][j]=t;};

///now the next column of L

for(int j=i+1;j<=std::min<int>(dim-1,i+p);j++){
adj=std::max<int>(0,j-p);
t=M[j][i-adj];
for(int k=0;k<i-adj;k++)t-=L[j][k]*U[k][i-k];

L[j][i-adj]=t/U[i][0];

};

};

vct<T,dim> y_1;
vct<T,dim> x;

/// we have LU(x)=b, solve first for Ly_1=b then solve for
/// U(x)=y_1;

for(int i=0;i<dim;i++){
t=b[i];
for(int k=0;k<L[i].size()-1;k++)t-=L[i][k]*y_1[k+std::max<int>(0,i-p)];
y_1[i]=t;
};

///now solve for U(x)=y_1:

for(int i=dim-1;i>=0;i--){
t=y_1[i];
for(int k=U[i].size()-1;k>=1;k--)t-=U[i][k]*x[k+i];

x[i]=t/U[i][0];};


return x;

};
