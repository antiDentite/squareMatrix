This is a cleaned up version of a project that I completed last year when I was learning C and C++ for scientific computing. There are three main matrix template classes that one can use and a vector template class. The matrix classes are for square full matrices, square band matrices and square tridiagonal matrices. The template parameters enable one to define the dimension and scalar type for both the vector class and the matrix classes. All standard linear algebra operations are supported as well as the transpose operator and a linear system solver. The linear system solver uses the LU decomposition and currently only works when the underlying matrix is deemed to be invertible. Have a look at the sample file, sample.cpp, for an overview of how to use these classes. 


## **Usage**

Simply insert the header file SQUARE.h and the source file SQUARE.cpp into your main source file. Remember to include the source file after the header file. 

## **Ongoing Work**

In the near future, I wish to add the following features:

* A symmetric matrix class that will make use of the Cholesky factorization.

* A slowSolve method that impliments the general LU factorization. I already have the code written for this feature but need to debug it before making a new commit. 

* Add a section in the sample file on using multi-precision floating numbers such as mpfr. 


## **Contact Information**

If you have any questions or comments about this project, you can contact Marc Fortier at reitrof.cram@gmail.com



