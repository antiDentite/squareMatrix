#include "SQUARE.h"
#include "SQUARE.cpp"

using namespace std;

int main(void){

/// 1. vct<T,dim> template class

///Vector of type double and dimension 5
vct<double,5> v_1;
///Vector of type int and dimension 10
vct<int,10> v_2;
///Change the first coordinate of v_1:
v_1(0)=3.23;

///Define vector of type double and dimension 4
///with all entries equal to 2.5:

vct<double,4> v_3(2.5);

///Define a vector of type double with 4 randomly
///chosen coordinate values between -2 and 2:

vct<double,4> v_4=vct<double,4>::Random(time(NULL),-2,2); ///the first value is the seed.
cout<<v_4;
///Add, subtract, dot product, equate, scalar multiply:

v_4*v_4; v_3+v_4; v_3-v_4;v_3=v_4;3.27*v_3;

/// 2. The Matrix classes

///2.1 The FULL<T,dim> template class

///Define square full matrix of type int and
///dimension 5

FULL<int,5> M_1;

///Define square full matrix of type double and
///dimension 7 with random entries between
///1.5 and 2.75

FULL<double,7> M_2=FULL<double,7>::Random(time(NULL),1.5,2.75);
cout<<"\n \n";
cout<<M_2;

///Take the transpose
cout<<"\n \n";
cout<<transpose(M_2);

///Change the first entry of the first row of M_2:
M_2(0,0)=2.3;
///Change the second entry of the fifth row:
M_2(4,1)=3.1;

/// 2.2 The Band<T,dim> template class

///Recall that a matrix B is a band matrix with
///parameters (p,q) whenever B(i,j)=0 for
/// j-i>q or i-j<p.

///Declaring a square band matrix of type double
///dimension 4 and parameters p=1,q=2;

BAND<double,4> B_1(1,2);

///...with random coefficients between -1.5 and 2.75

BAND<double,4> B_2=BAND<double,4>::Random(1,2,time(NULL),-1.5,2.75);
cout<<"\n \n";
cout<<B_2;

///2.3 The TRIDAG<T,dim> template class

///Recall that a matrix T is tridiagonal whenever
///T(i,j)=0 for |i-j|>1.

///Declaring a square tridiagonal matrix of type
///double and dimension 5:

TRIDAG<double,5> T_1;

///...with random coefficients between -5 and 5:

TRIDAG<double,5> T_2=TRIDAG<double,5>::Random(time(NULL),-5,5);
cout<<"\n \n";
cout<<T_2;

///3. Linear Algebra

FULL<double,5> M_3=FULL<double,5>::Random(time(NULL),-5,5);
BAND<double,5> B_3=BAND<double,5>::Random(1,2,time(NULL),-5,5);
TRIDAG<double,5> T_3=TRIDAG<double,5>::Random(time(NULL),-5,5);
vct<double,5> v_5=vct<double,5>::Random(time(NULL),-5,5);

///Multiplication:

B_3*B_3;B_3*M_3;T_3*M_3;T_3*B_3;M_3*M_3;
B_3*M_3*T_3;

///Addition

B_3+B_3;B_3+M_3;B_3+T_3;T_3+B_3;M_3+M_3;

///Substraction

B_3-B_3;B_3-M_3;M_3-B_3;M_3-T_3;T_3-B_3;

///Scalar and vector multiplication

B_3*v_5;M_3*v_5;T_3*v_5;v_5*v_5;
2.53*B_3; 2.1*M_3; -M_3; -3.2*T_3;
T_3*1.243;M_3*1.3; B_3*2.1;

///4. LU Decomposition and solving linear systems

///4.1 decompBin struct

///The struct decompBin can be used to
///attempt to find an LU decomposition for
///a square matrix of class FULL, BAND or TRIDAG.
///It only works when the underlying matrix is
///found to be invertible. The struct has a member
///function "solve" that can be used to solve
///a linear system. The following examples should
///clarify the process:

///constructors
FULL<double,5>::decompBin bin_1(M_3);
BAND<double,5>::decompBin bin_2(B_3);
TRIDAG<double,5>::decompBin  bin_3(T_3);

///The member variable bool valid will indicate
///whether we could find an LU decomposition:
cout<<"\n \n";
cout<<bin_1.valid<<','<<bin_2.valid<<','<<bin_3.valid;
cout<<"\n \n";
///If the  member variable "valid" is 1 then
///an LU decomposition has been found.
/// The constructor tries to build L and U and
///then compares the values of the diagonal
///of U with a threshold value that is
///a static member variable of a template class
/// template <class T> class THRESH.

///In this case, we have the following declaration
///at the end of the LinAlg.cpp file:
/// template<> double THRESH<double>::value=1.0e-20;

///When one of the diagonal elements of U is less
/// than 1.0e-20, then valid is set to 0 and the
/// decomposition "failed". One can simply change that
///value for double scalar type or even for custom
/// multi-precision scalar types.

/// 4.2 solving linear systems

///Once the decompBin struct has been defined
///we can attempt to solve any linear system of
///the form Mx=b with M a square matrix.
///If valid is 0 and you attempt to use "solve"
///, then the process will exit and an error message
/// will alert you. You should do
///something like this:

vct<double,5> v_6;

if(bin_1.valid){
v_6=bin_1.solve(v_5);
cout<<"Error is ="<<v_5-M_3*v_6<<'\n';
};

if(bin_2.valid){
v_6=bin_2.solve(v_5);
cout<<"Error is ="<<v_5-B_3*v_6<<'\n';
};

if(bin_3.valid){

v_6=bin_3.solve(v_5);
cout<<"Error is ="<<v_5-T_3*v_6<<'\n';

};

///Generating an error when calling solve:
FULL<double,5> M_4;

FULL<double,5>::decompBin bin_4(M_4);

bin_4.solve(v_5);



return 0;}
