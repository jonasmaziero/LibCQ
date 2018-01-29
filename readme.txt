LibCQ: C Library for Quantum Information Science

In this set of functions, I'm translating LibForQ and LibForro to C. I'm also adding 
some new functionalities that are more convenient to implement in C. Further explanations and
references can be found in the comments in the original Fortran code.


Disclaimer:
Anyone can use this code, respecting the licencing terms of course. But I should observe that I
make the code for my own use and for the use of my research group. So, I do not offer any warranty 
for it. Also, I usually make some tests before releasing the code, but nothing too significant. 
So I shall thank any feedback.


Remarks:
- For using functions needing Lapack and Blas use the command:
$ gcc *.c -llapacke -lblas -lm
The -lm is used to link the math functions of C.
- To compile the CUDA code use:
$ nvcc name.cu 