/*----------- JACOBI'S ITERATION METHOD TO SOLVE LINEAR EQUATIONS -----*/

/*	THE PROGRAM SOLVES THE SYSTEM OF LINEAR EQUATIONS USING

	JACOBI'S ITERATION METHOD.

	INPUTS :  1) Number of variables in the equation.

		  2) Coefficient's of linear equations.

	OUTPUTS : Results of every iteration till 'q' is pressed.               */

	/*---------------------------    PROGRAM  -----------------------------*/
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <conio.h>
#include <omp.h>
#include <algorithm> 
#include <cmath>

using namespace std;

int main()
{

	std::cout << "\n        Computational Techniques - J. S. CHITODE" <<
		"\n   JACOBI'S ITERATION METHOD TO SOLVE LINEAR EQUATIONS" <<
		"\n\n             The form of equations is as follows\n\n" <<
		"                 a11x1 + a12x2 + ... + a1nxn = b1\n" <<
		"                 a21x1 + a22x2 + ... + a2nxn = b2\n" <<
		"                 a31x1 + a32x2 + ... + a3nxn = b3\n" <<
		"                 ................................\n" <<
		"                 an1x1 + an2x2 + ... + annxn = bn\n";
	//Ініціалізація системи рівнянь


	double a[11][11], x[10], y[10];
	double coef[3][3] = { {2,5,7},{6,3,4},{5,-2,-3} };
	double rmatrix[10] = { 1,0,0,0,1,0,0,0,1,0 };
	//double e = 0.00000000001;
	double e = 0.000008;

	double bufferarray[10];

	double priorX[10];

	double maxValue = -1000.0;
	double count = 10000;
	/* ARRAY OF a[n][n] STORING COEFFICIENTS OF EQUATIONS */

	int n;
	/* ENTER THE NUMBER OF VARIABLES IN THE EQUATION   */
	std::cout << "\n\nEnter the number of variables (max 9) = ";
	std::cin >> n;
	//Введення елементів матриці
	for (int i = 0; i < sqrt(n); i++)
	{
		for (int j = 0; j < sqrt(n); j++)
		{
			std::cout << "a" << i + 1 << j + 1 << ": " << coef[i][j]<<endl;
			//cin >> coef[i][j];
		}

	}
	std::cout << endl << "Введена матриця: " << endl;
	for (int i = 0; i < sqrt(n); i++)
	{
		for (int j = 0; j < sqrt(n); j++)
		{
			std::cout << coef[i][j] << " ";

		}
		std::cout << endl;
	}

	std::cout << "Введіть одиничну матрицю: " << endl;
	for (int i = 0; i < n; i++)
	{
		printf("b%d = ", i);
		cout << rmatrix[i] << endl;
		//cin >> rmatrix[i];
	}

	int buffer = 1;
	int z = 1;
	int m = 0;
	std::cout << "Кількість потоків: ";
	int numberOfThreads = 0;
	cin >> numberOfThreads;
	omp_set_num_threads(numberOfThreads);
	double start = omp_get_wtime();

#pragma omp parallel sections 
	{
		//#pragma omp parallel for shared(a) private(i, j)
#pragma omp section
#pragma omp parallel for collapse(2)
		for (int i = 0; i < 11; i++)
		{
			for (int j = 0; j < 11; j++)
			{
				a[i][j] = 0;
			}

		}
#pragma omp section
#pragma omp parallel for
		for (int i = 0; i < 10; i++)
		{
			bufferarray[i] = 0;
			priorX[i] = 0;
		}
	}
	//Перезапис у систему рівнянь
   // #pragma omp parallel for
	for (int i = 1; i <= n; i++)
	{
		buffer = z;
		// #pragma omp critical
		{
			for (int j = 0; j < sqrt(n); j++)
			{

				a[i][z] = coef[m][j];
				z += sqrt(n);
			}
		}
#pragma omp parallel sections 
		{
#pragma omp section 
			{
				int r = sqrt(n);
				if (i % r == 0)
				{
					m++;
					z = 1;
				}
				else {
					z = buffer;
					z++;
				}
			}
			//printf("b%d = ", i);
			//cin >> a[i][n + 1];
#pragma omp section 
			{
				a[i][n + 1] = rmatrix[i - 1];
				x[i] = y[i] = 0;
			}
		}
	}
#pragma omp parallel sections
	{
#pragma omp section
		{
			std::cout << endl << "Отримана система рівнянь: " << endl;
			for (int i = 0; i < 11; i++)
			{
				for (int j = 0; j < 11; j++)
				{
					std::cout << a[i][j] << " ";
				}
				std::cout << endl;
			}
		}

#pragma omp section
		{
			//Знаходження коренів системи (елементів оберненої матриці)
			//int i;
			printf("\n\nThe results are as follows....\n\n");
		}
#pragma omp section
		{

			while (abs(maxValue) > e)
			{
#pragma omp parallel for 
				for (int i = 1; i <= n; i++)
				{
					/*  LOOP TO CALCULATE VALUES OF x1,x2,...,xn etc  */
					double root = 0;
					int j = 0;
					//std::cout << "root" << root << endl;

				   // #pragma omp parallel for 
					for (j = 1; j <= n; j++)
					{

						if (i != j) {
							//double sum = a[i][j] * y[j];
							//#pragma omp critical
#pragma omp critical
							{

								root -= a[i][j] * y[j];
							}
							//std::cout << root << " " << a[i][j] << " " << y[j] << endl;
							//std::cout << "i: " << i << "j: " << j<<endl;
							//std::cout << a[i][j] << " " << y[j]<<endl;


						}
					}
					//cout << "final: " << root;
#pragma omp critical
					{
						x[i] = root + a[i][j];
						x[i] = x[i] / a[i][i];
					}
				}
#pragma omp parallel for 
				for (int i = 1; i <= n; i++)
				{
#pragma omp parallel sections
					{
#pragma omp section
						{
							/* LOOP TO PRINT VALUES OF x1,x2,...xn etc */
							bufferarray[i - 1] = abs(x[i] - priorX[i - 1]);
							y[i] = x[i];
							priorX[i - 1] = x[i];
							x[i] = 0;
						}
#pragma omp section
						{
							//#pragma omp task
														//{
							printf("x%d = %.8lf  \n", i, y[i]);
							//}
						}

					}

				}
				std::cout << "-----------------------\n";
				//maxValueOld = maxValue;

				maxValue = bufferarray[0];
#pragma omp parallel for
				for (int i = 0; i < n; i++)
				{
					//#pragma omp critical
					{
						if (bufferarray[i] > maxValue)
							maxValue = bufferarray[i];
					}
				}

			}
		}
	}
	double end = omp_get_wtime();
	double time = end - start;
	std::cout << endl << "Час: " << time << endl;
#pragma omp parallel
	{
		std::cout << "Кількість процесів: " << omp_get_num_threads() << endl;
	}
	std::cout << "Обернена матриця: " << endl;
	for (int i = 0; i < n; i++)
	{
		std::cout << priorX[i] << " ";
		if ((i + 1) % int(sqrt(n)) == 0)
		{
			std::cout << endl;
		}
	}

	return 0;

}
/*-------------------------------- END OF PROGRAM -----------------------*/

/*----------- JACOBI'S ITERATION METHOD TO SOLVE LINEAR EQUATIONS -----*/

/*	THE PROGRAM SOLVES THE SYSTEM OF LINEAR EQUATIONS USING

	JACOBI'S ITERATION METHOD.

	INPUTS :  1) Number of variables in the equation.

		  2) Coefficient's of linear equations.

	OUTPUTS : Results of every iteration till 'q' is pressed.               */

	/*---------------------------    PROGRAM  -----------------------------*/
