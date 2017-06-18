1. Compile the LU decomposition example

gfortran -O3 hs_lu_decomp.f90 -llapack -o hs.x

2. Execute code with some n = (100, 1000) check the times

./hs.x > plot.dat

3. Plot data file

gnuplot plot.gnu

4. Check the execution time in case we use tridiagonal solver


4.1 Compile example

gfortran -O3 hs_tridiag.f90 -llapack -o hs.x

4.2 Execute code for same set of n values (100, 1000, 10000) compare the times with Step 2.

5. Loop alighment example. In source folder there is two example codes which implement Jacobi iterative scheme

5.1 hs_jacobi.f90 - loop nest with poor cache use

5.2 hs_jacobi_loop_optimized.f90 - loop nest with sride-one access

6. Review the above examples, compile and compare the times for n (100, 200, 400)

gfortran -O3 hs_jacobi.f90 -o hs.x
./hs.x > plot.dat
 
gfortran -O3 hs_jacobi_loop_optimized.f90 -o hs.x
./hs.x > plot.dat

7. Comment the results

8. Check Mathematica example for problem size n = 5000 and compare the results with Gauss and Tridiagonal solvers

module load mathematica
MathKernel < hs_mathematica.m

9. Which code perform worst and wchich is fastest? Why? 
