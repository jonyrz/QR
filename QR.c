#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
typedef struct {
	int m, n;
	double ** v;
} mat_t, *mat;
 
mat matriz_nueva(int m, int n)
{
	mat x = malloc(sizeof(mat_t));
	x->v = malloc(sizeof(double*) * m);
	x->v[0] = calloc(sizeof(double), m * n);
	for(int i=0; i<m; i++)
		x->v[i] = x->v[0] + n * i;
	x->m = m;
	x->n = n;
	return x;
}
 
void matriz_borrar(mat m)
{
	free(m->v[0]);
	free(m->v);
	free(m);
}
 
void matriz_transpuesta(mat m)
{
	for (int i = 0; i < m->m; i++) {
		for (int j = 0; j < i; j++) {
			double t = m->v[i][j];
			m->v[i][j] = m->v[j][i];
			m->v[j][i] = t;
		}
	}
}
 
mat matriz_copiar(int n, double a[][n], int m)
{
	mat x = matriz_nueva(m, n);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			x->v[i][j] = a[i][j];
	return x;
}
 
mat matriz_multi(mat x, mat y)
{
	if (x->n != y->m) return 0;
	mat r = matriz_nueva(x->m, y->n);
	for (int i = 0; i < x->m; i++)
		for (int j = 0; j < y->n; j++)
			for (int k = 0; k < x->n; k++)
				r->v[i][j] += x->v[i][k] * y->v[k][j];
	return r;
}
 
mat matriz_menor(mat x, int d)
{
	mat m = matriz_nueva(x->m, x->n);
	for (int i = 0; i < d; i++)
		m->v[i][i] = 1;
	for (int i = d; i < x->m; i++)
		for (int j = d; j < x->n; j++)
			m->v[i][j] = x->v[i][j];
	return m;
}
 
/* c = a + b * s */
double *vmadd(double a[], double b[], double s, double c[], int n)
{
	for (int i = 0; i < n; i++)
		c[i] = a[i] + s * b[i];
	return c;
}
 
/* m = I - v v^T */
mat vmul(double v[], int n)
{
	mat x = matriz_nueva(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			x->v[i][j] = -2 *  v[i] * v[j];
	for (int i = 0; i < n; i++)
		x->v[i][i] += 1;
 
	return x;
}
 
/* ||x|| */
double vnorm(double x[], int n)
{
	double sum = 0;
	for (int i = 0; i < n; i++) sum += x[i] * x[i];
	return sqrt(sum);
}
 
/* y = x / d */
double* vdiv(double x[], double d, double y[], int n)
{
	for (int i = 0; i < n; i++) y[i] = x[i] / d;
	return y;
}
 
/* take c-th column of m, put A v */
double* mcol(mat m, double *v, int c)
{
	for (int i = 0; i < m->m; i++)
		v[i] = m->v[i][c];
	return v;
}
 
void matriz_mostrar(mat m)
{
	for(int i = 0; i < m->m; i++) {
		for (int j = 0; j < m->n; j++) {
			printf(" %8.3f", m->v[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}
 
void algoritmo(mat m, mat *R, mat *Q)
{
	mat q[m->m];
	mat z = m, z1;
	for (int k = 0; k < m->n && k < m->m - 1; k++) {
		double e[m->m], x[m->m], a;
		z1 = matriz_menor(z, k);
		if (z != m) matriz_borrar(z);
		z = z1;
 
		mcol(z, x, k);
		a = vnorm(x, m->m);
		if (m->v[k][k] > 0) a = -a;
 
		for (int i = 0; i < m->m; i++)
			e[i] = (i == k) ? 1 : 0;
 
		vmadd(x, e, a, e, m->m);
		vdiv(e, vnorm(e, m->m), e, m->m);
		q[k] = vmul(e, m->m);
		z1 = matriz_multi(q[k], z);
		if (z != m) matriz_borrar(z);
		z = z1;
	}
	matriz_borrar(z);
	*Q = q[0];
	*R = matriz_multi(q[0], m);
	for (int i = 1; i < m->n && i < m->m - 1; i++) {
		z1 = matriz_multi(q[i], *Q);
		if (i > 1) matriz_borrar(*Q);
		*Q = z1;
		matriz_borrar(q[i]);
	}
	matriz_borrar(q[0]);
	z = matriz_multi(*Q, m);
	matriz_borrar(*R);
	*R = z;
	matriz_transpuesta(*Q);
}
 
double A[][3] = {
	{ 12, -51,   4},
	{  6, 167, -68},
	{ -4,  24, -41},
	{ -1, 1, 0},
	{ 2, 0, 3},
};
 
int main()
{
	mat R, Q;
	mat x = matriz_copiar(3, A, 5);
	algoritmo(x, &R, &Q);
 
	puts("Q"); matriz_mostrar(Q);
	puts("R"); matriz_mostrar(R);
 
	// to show their product is the input matrix
	mat m = matriz_multi(Q, R);
	puts("Q * R"); matriz_mostrar(m);
 
	matriz_borrar(x);
	matriz_borrar(R);
	matriz_borrar(Q);
	matriz_borrar(m);
	return 0;
}