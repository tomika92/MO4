// zad3.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

void process(int N, double q, double r, double s, double **A, double *b)
{
    for (int i = 0; i <= N; i++)
    {
        for (int j = 0; j <= N; j++)
        {
            A[i][j] = 0;
        }
        b[i] = 0;
    }

    double h = 5. / N;
    A[0][0] = 1;
    A[N][N] = 1;
    for (int i = 1; i < N; i++)
    {
        A[i][i - 1] = (1. / (h * h)) - (q / (2 * h));
        A[i][i] = r - (2. / (h * h));
        A[i][i + 1] = (1. / (h * h)) + (q / (2 * h));
        b[i] = s;
    }

}
void Gauss(double **A, double *b, int N)
{
    for (int k = 0; k < N - 1; k++) //dekompozycja LU
    {
        for (int i = k + 1; i < N; i++)
        {
            double aux = A[i][k] / A[k][k];
            for (int j = k + 1; j <= N; j++)
            {
                A[i][j] = A[i][j] - A[k][j] * aux;
            }
            A[i][k] = aux;
        }
    }
    for (int k = 0; k < N - 1; k++)//eliminacja w przod
    {
        for (int i = k + 1; i < N; i++)
        {
            b[i] = b[i] - b[k] * A[i][k];
        }
    }
    b[N - 1] = b[N - 1] / A[N - 1][N - 1];//podstawienie wstecz
    for (int i = N - 2; i >= 0; i--)
    {
        double s = 0;
        for (int j = i + 1; j < N; j++)
        {
            s = s + A[i][j] * b[j];
        }
        b[i] = (b[i] - s) / A[i][i];
    }
    cout << setprecision(6) << scientific;
    for (int i = 0; i < N; i++)
    {
        cout << "Krok " << i << " = " << b[i] << "\t\t" << endl;
    }
}

int main()
{
    double q = 0.41;
    double r = -0.49;
    double s = 1.88;
    int N = 19;
    double** A;
    A = new double* [N + 1];
    double* b;
    b = new double[N + 1];
    for (int i = 0; i <= N; i++)
    {
        A[i] = new double[N + 1];
    }
    cout << "Rozwiazanie dla siatki z N odcinkow" << endl;
    process(N, q, r, s, A, b);
    Gauss(A, b, N+1);
    
    N = 2*N;
    double** A2;
    A2 = new double* [N + 1];
    double* b2;
    b2 = new double[N + 1];
    for (int i = 0; i <= N; i++)
    {
        A2[i] = new double[N + 1];
    }
     cout << "\n\nRozwiazanie dla siatki z 2N odcinkow" << endl;
     process(N, q, r, s, A2, b2);
     Gauss(A2, b2, N + 1);

     N = 2 * N;
     double** A4;
     A4 = new double* [N + 1];
     double* b4;
     b4 = new double[N + 1];
     for (int i = 0; i <= N; i++)
     {
         A4[i] = new double[N + 1];
     }
     cout << "\n\nRozwiazanie dla siatki z 4N odcinkow" << endl;
     process(N, q, r, s ,A4, b4);
     Gauss(A4, b4, N + 1);
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
