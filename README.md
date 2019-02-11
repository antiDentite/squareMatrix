This is a small project that I completed while learning scientific computing in C++. My hope is that others find it useful in their journey into C++. There are three main matrix template classes and a vector template class. The matrix classes are for square full matrices, square band matrices and square tridiagonal matrices. The template parameters enable one to define the dimension and scalar type for both the vector class and the matrix classes. All standard linear algebra operations are supported as well as the transpose operator and a linear system solver. The linear system solver uses the LU decomposition and currently only works when the underlying matrix is deemed to be invertible. The sample file (sample.cpp) gives an overview of how to use these classes. 


## **Usage**

Simply insert the header file SQUARE.h and the source file SQUARE.cpp into your main source file. Remember to include the source file after the header file. 

## **Contact Information**

If you have any questions or comments about this project, you can contact Marc Fortier at reitrof.cram@gmail.com



