This is a cleaned up version of a project that I completed last year when I was learning C and C++ for scientific computing. The are three main matrix template classes that one can use and a vector template class. The matrix classes are for square full matrices, square band matrices and square tridiagonal matrices. The template parameters enable one to define the dimension and scalar type for both the vector class and the matrix classes. All standard linear algebra operations are supported as well as the transpose operator and a linear system solver. The linear system solver uses the LU decomposition and currently only works when the underlying matrix is assumed to be invertible (the diagonal elements of U are not zero). Have a look at the sample file, sample.cpp, for an overview of how to use the classes. 


## **Installation**

Simply insert the header file LinAlg.h in your main source file. This header file contains declarations and full definitions for the friend functions of the various classes. The full definitions for the member functions are in the LinAlg.cpp file which is included in the LinAlg.h. 

## **Ongoing Work**

In the near future, I will add the following features:

* A symmetric matrix class that will make use of the Cholesky factorization.

* A slowSolve method that impliments the general LU factorization. I already have the code written for this feature but need to debug it before making a new commit. 

* Add a section in the sample file on using multi-precision floating numbers such as mpfr. 


## **Contact Information**

If you have any questions about this project, you can contact Marc Fortier at reitrof.cram@gmail.com



