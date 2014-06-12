#ifndef LINALG_H_INCLUDED
#define LINALG_H_INCLUDED

#include <functional>
#include <algorithm>
#include <iostream>
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>

template <class T, unsigned int dim> class vct:public std::vector<T>{

public:

vct();
~vct();
vct(const vct<T, dim> &);
vct(const T  &);

T & operator()(const unsigned int &);

friend vct<T, dim>  operator -(const vct<T,dim> &v){

vct v_1;
std::transform(v.begin(),v.end(),v_1.begin(),std::negate<T>());
return v_1;

};
vct<T,dim> & operator +=(const vct<T,dim> &);
vct<T,dim> & operator -=(const vct<T,dim> &);

vct<T,dim> operator+(const vct<T,dim> & ) const;
vct<T,dim> operator-(const vct<T,dim>  &)const;

T operator *(const vct<T,dim> &)const;

friend vct<T,dim> operator *(const  T & s, const vct<T,dim> & v){

vct v_1(s);

std::transform(v.begin(),v.end(),v_1.begin(),v_1.begin(),std::multiplies<T>());

return v_1;
};

static vct<T,dim> Random(const unsigned int &,const double & a, const double & b);
friend std::ostream & operator << (std::ostream & os,const vct<T,dim> & v){

for(int i=0;i<dim-1;i++){os<<v[i]<<',';};
os<<v[dim-1];

return os;

};

};

template<class T,class S> class SQUARE{

public:

const S & return_to_derived_2()const;
T entry(unsigned int, unsigned int) const;
///Arithmetic checks

template <class Z> std::vector< std::vector<T> > slowAdd (const SQUARE<T,Z> &)const;
template <class Z> std::vector< std::vector<T> > slowMinus (const SQUARE<T,Z> &)const;
template <class Z> std::vector< std::vector<T> > slowMult (const SQUARE<T,Z> &)const;

/// slow transpose
 std::vector< std::vector<T> > slowTranspose()const;

///negation

std::vector< std::vector<T> > negation()const;

///scalar multiplication

std::vector< std::vector<T> > scalarMult(const T & t)const;

///Streamer

friend std::ostream & operator <<(std::ostream & os,const SQUARE<T,S> & M_1){
unsigned int dim=M_1.return_to_derived_2().M.size();
for(int i=0;i<dim;i++){
for(int j=0;j<dim-1;j++)os<<M_1.entry(i,j)<<',';
os<<M_1.entry(i,dim-1);
os<<'\n';

};

return os;

};};

///forward declaration of some classes:
template<class T, unsigned int dim> class FULL;
template<class T,unsigned int dim> class BAND;
template<class T, unsigned int dim> class TRIDAG;

template <class T,unsigned int dim> class FULL :public SQUARE<T,FULL<T,dim> >{

public:

std::vector< std::vector<T> > M;

FULL();
~FULL();
T & operator()(const unsigned int & i, const unsigned int & j);
T entry(const unsigned int & i, const unsigned int & j)const;

///equality

FULL<T,dim>  & operator =(const BAND<T,dim> &);
FULL<T,dim> & operator =(const TRIDAG<T,dim> &);

///negation

friend FULL<T,dim> operator -(const FULL<T,dim> & M_1){

FULL<T,dim> M_2; M_2.M=const_cast<FULL<T,dim> &>(M_1).negation();

return M_2;

};

///Addition

FULL<T,dim>  operator + (const FULL<T,dim> & ) const;
FULL<T,dim>  operator + (const BAND<T,dim> & )const ;
FULL<T,dim>  operator + (const TRIDAG<T,dim> & ) const ;

///substraction

FULL<T,dim>  operator - (const FULL<T,dim> &)const ;
FULL<T,dim>  operator - (const BAND <T,dim> &)const ;
FULL<T,dim>  operator - (const TRIDAG<T,dim> &)const ;

///Multiplication operators

vct<T,dim>  operator * (const vct<T,dim> &)const;
FULL<T,dim>  operator *(const FULL<T,dim> &)const;
FULL<T,dim>  operator *(const BAND <T,dim> &)const;
FULL<T,dim>  operator *(const TRIDAG <T,dim> &)const;

friend FULL<T,dim> operator *(const FULL<T,dim> & M_1,const T & t){

FULL<T,dim> M_2;

M_2.M=const_cast<FULL<T,dim> &>(M_1).scalarMult(t);

return M_2;

};

friend FULL<T,dim> operator *(const T & t, const FULL<T,dim> & M_1){
return M_1*t;
};

///transpose

friend FULL<T,dim> transpose(const FULL<T,dim> & M_1){

FULL<T,dim> M_2;

for(int i=0;i<dim;i++){

for(int j=0;j<dim;j++){M_2.M[i][j]=M_1.M[j][i];};

};

return M_2;
}

///Solvers
vct<T,dim>  fastSolve(const vct<T,dim> &)const ;
vct<T,dim>  slowSolve(const vct<T,dim> &)const; ///not supported yet

///Generating Matrices with Random Coefficients

static FULL<T,dim> Random(unsigned int,double, double);

///Streamer

};

template <class T,unsigned int dim>  class BAND: public SQUARE<T,BAND<T,dim> >{

public:

std::vector< std::vector<T> > M;
int p,q;

BAND(const unsigned int &,const unsigned int &);
BAND();
~BAND();

T & operator()(const unsigned int &, const unsigned int &);
T entry(const unsigned int &, const unsigned int &)const;

///equality

BAND<T,dim> & operator =(const TRIDAG<T,dim> &);

///negation

friend BAND<T,dim> operator -(const BAND<T,dim> & M_1){

BAND <T,dim> M_2(M_1.p,M_1.q); M_2.M=const_cast<BAND<T,dim> &>(M_1).negation();

return M_2;

};

///Addition

FULL<T,dim>  operator +(const FULL<T,dim> &)const ;
BAND<T,dim>  operator +(const BAND<T,dim> &)const ;
BAND<T,dim>  operator +(const TRIDAG<T,dim> &)const ;

///substraction

FULL<T,dim> operator -(const FULL<T,dim> &)const;
BAND<T,dim>  operator -(const BAND<T,dim> &)const;
BAND<T,dim>  operator-(const TRIDAG<T,dim> &)const;

///multiplication

vct<T,dim> operator *(const vct<T,dim> &) const;
FULL<T,dim>  operator *(const FULL<T,dim> &)const;
BAND<T,dim>  operator *(const BAND<T,dim> &)const;
BAND<T,dim>  operator *(const TRIDAG<T,dim> &)const;

friend BAND<T,dim> operator *(const BAND<T,dim> & M_1,const T & t){

BAND<T,dim> M_2(M_1.p,M_1.q);

M_2.M=const_cast< BAND<T,dim> > (M_1).scalarMult(t);

return M_2;

};

friend BAND<T,dim> operator *(const T &  t, const BAND<T,dim> &  M_1){
return M_1*t;
};

///transpose
friend BAND<T,dim>  transpose(const BAND<T,dim> & M_1){

const unsigned int q(M_1.q);
const unsigned int p(M_1.p);

BAND M_2(q,p);

for(int i=0;i<dim;i++){
int adj=std::max<int>(0,i-q);
for(int j=adj;j<=std::min<int>(dim-1,i+p);j++){

M_2.M[i][j-adj]=M_1.M[j][i-std::max<int>(0,j-p)];

};

};


return M_2;

};


///generating matrices with random coefficients

static BAND<T,dim> Random(unsigned int, unsigned int,unsigned int, double, double);

vct<T,dim> fastSolve(const vct<T,dim> &)const;

};

template <class T,unsigned int dim> class TRIDAG: public SQUARE<T,TRIDAG<T,dim> >{

public:

std::vector< std::vector<T> > M;

TRIDAG();
~TRIDAG(){};

T & operator()(const unsigned int &, const unsigned int &);
T entry(const unsigned int &, const unsigned int &)const;

///negation

friend TRIDAG<T,dim> operator -(const TRIDAG<T,dim> & M_1){

TRIDAG <T,dim> M_2; M_2.M=const_cast<TRIDAG<T,dim> &>(M_1).negation();

return M_2;

};

///addition

TRIDAG<T,dim>  operator +(const TRIDAG<T,dim> &)const;
BAND<T,dim>  operator +(const BAND<T,dim> &)const;
FULL<T,dim>  operator +(const FULL<T,dim> &)const;

///substraction

TRIDAG<T,dim>  operator -(const TRIDAG<T,dim> &)const;
BAND<T,dim>  operator -(const BAND<T,dim> &)const;
FULL<T,dim>  operator -(const FULL<T,dim> &)const;

///multiplication
vct<T,dim>  operator *(const vct<T,dim> &)const;
FULL<T,dim>  operator *(const FULL<T,dim> &)const;
BAND<T,dim>  operator *(const BAND<T,dim> &)const;
BAND<T,dim>  operator *(const TRIDAG<T,dim> & )const;

///scalar multiplication

friend TRIDAG<T,dim> operator *(const TRIDAG<T,dim> M_1,const T t){

TRIDAG<T,dim> M_2;

M_2.M=const_cast<TRIDAG<T,dim>  &>(M_1).scalarMult(t);

return M_2;

};

friend TRIDAG<T,dim> operator *(const T  t, const TRIDAG<T,dim>  M_1){
return M_1*t;
};

///transpose

friend TRIDAG<T,dim>  transpose(const TRIDAG<T,dim>  M_1){

TRIDAG M_2;

M_2.M[1][0]=M_1.M[0][1];
M_2.M[0][0]=M_1.M[0][0];
M_2.M[0][1]=M_1.M[1][0];

for(int i=1;i<dim-1;i++){

M_2.M[i][1]=M_1.M[i][1];
M_2.M[i][2]=M_1.M[i+1][0];
M_2.M[i+1][0]=M_1.M[i][2];

};

M_2.M[dim-1][1]=M_1.M[dim-1][1];

return M_2;

};

///solving a system:

vct<T,dim>  fastSolve(const vct<T,dim> &)const;

///Generating Matrices with random coefficients:

static TRIDAG<T,dim> Random(unsigned int, double, double);

};


#endif // LINALG_H_INCLUDED
