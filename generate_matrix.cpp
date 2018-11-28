
#include <cstdlib>


double**  gen_matrix2(int n, int m)
{
	double **d;

	d = (double**)malloc(sizeof(double*)*n);

	for (int i3 = 0; i3 < n; i3++)
	{
		d[i3] = (double*)malloc(sizeof(double)*m);
	}

	return d;
}
double*** gen_matrix3(int n, int m, int m2)
{
	double ***d;

	d = (double***)malloc(sizeof(double**)*n);


	for (int i3 = 0; i3 < n; i3++)
	{
		d[i3] = gen_matrix2(m, m2);
	}

	return d;
}

double**** gen_matrix4(int n, int m, int m2, int m3)
{
	double ****d;

	d = (double****)malloc(sizeof(double***)*n);

	for (int i3 = 0; i3 < n; i3++)
	{
		d[i3] = gen_matrix3(m, m2, m3);
	}

	return d;
}

double***** gen_matrix5(int n, int m, int m2, int m3, int m4)
{
	double *****d;

	d = (double*****)malloc(sizeof(double****)*n);

	for (int i4 = 0; i4 < n; i4++)
	{
		d[i4] = gen_matrix4(m, m2, m3, m4);
	}

	return d;
}

double****** gen_matrix6(int n, int m, int m2, int m3, int m4, int m5)
{
	double ******d;

	d = (double******)malloc(sizeof(double*****)*n);

	for (int i5 = 0; i5 < n; i5++)
	{
		d[i5] = gen_matrix5(m, m2, m3, m4, m5);
	}

	return d;
}

double******* gen_matrix7(int n, int m, int m2, int m3, int m4, int m5, int m6)
{
	double *******d;

	d = (double*******)malloc(sizeof(double******)*n);

	for (int i6 = 0; i6 < n; i6++)
	{
		d[i6] = gen_matrix6(m, m2, m3, m4, m5, m6);
	}

	return d;
}