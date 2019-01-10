#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
	int m, n;
	double ** v;
} mat_t, *mat;

behavior matrix_transpose(mat m) 
{
	void main(void)
	{
		for (int i = 0; i < m->m; i++) {
			for (int j = 0; j < i; j++) {
				double t = m->v[i][j];
				m->v[i][j] = m->v[j][i];
				m->v[j][i] = t;
			}
		}
	}
}

behavior matrix_copy(int n, double a[][n], int m)
{
	mat main(void)
	{
		mat x;
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				x->v[i][j] = a[i][j];
		return x;
	}
}

behavior matrix_mul(mat x, mat y)
{
	mat main(void)
	{
		if (x->n != y->m) return 0;
		mat r;
		for (int i = 0; i < x->m; i++)
			for (int j = 0; j < y->n; j++)
                for (int k = 0; k < x->n; k++)
                        r->v[i][j] += x->v[i][k] * y->v[k][j];
		return r;
	}
}

behavior matrix_minor(mat x,int d)
{
    mat main(mat x, int d)
    {
        mat m = matrix_new(x->m, x->n);
        for (int i = 0; i < d; i++)
            m->v[i][i] = 1;
        for (int i = d; i < x->m; i++)
            for (int j = d; j < x->n; j++)
                m->v[i][j] = x->v[i][j];
        return m;
    }
}

/* c = a + b * s */
behavior vmadd(double a[], double b[], double s, double c[], int n)
{
    double* main(void)
    {
        for (int i = 0; i < n; i++)
            c[i] = a[i] + s * b[i];
        return c;
    }
}

/* m = I - v v^T */
behavior vmul(double v[], int n)
{
    mat main(void)
    {
        mat x = matrix_new(n, n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                x->v[i][j] = -2 * v[i] * v[j];
        for (int i = 0; i < n; i++)
            x->v[i][i] += 1;

        return x;
    }
}

/* ||x|| */
behavior vnorm(double x[], int n)
{
    double main(void)
    {
        double sum = 0;
        for (int i = 0; i < n; i++) sum += x[i] * x[i];
        return sqrt(sum);
    }
}

/* y = x / d */
behavior vdiv(double x[], double d, double y[], int n)
{
    double* main(void)
    {
	    for (int i = 0; i < n; i++) y[i] = x[i] / d;
	    return y;
    }
}

/* take c-th column of m, put in v */
behavior mcol(mat m, double *v, int c)
{
    double* main(void)
    {
        for (int i = 0; i < m->m; i++)
            v[i] = m->v[i][c];
        return v;
    }
}

behavior householder(mat m, mat *R, mat *Q)
{
    void main(void)
    {
        const int idx = m->m;
        mat q[idx];
        mat z = m, z1;
        for (int k = 0; k < m->n && k < m->m - 1; k++) {
            double e[m->m], x[m->m], a;
            z1 = matrix_minor(z, k);
            z = z1;

            mcol(z, x, k);
            a = vnorm(x, m->m);
            if (m->v[k][k] > 0) a = -a;

            for (int i = 0; i < m->m; i++)
                e[i] = (i == k) ? 1 : 0;

            vmadd(x, e, a, e, m->m);
            vdiv(e, vnorm(e, m->m), e, m->m);
            q[k] = vmul(e, m->m);
            z1 = matrix_mul(q[k], z);
            z = z1;
        }
        *Q = q[0];
        *R = matrix_mul(q[0], m);
        for (int i = 1; i < m->n && i < m->m - 1; i++) {
            z1 = matrix_mul(q[i], *Q);
            *Q = z1;
        }
        z = matrix_mul(*Q, m);
        *R = z;
        matrix_transpose(*Q);
    }
}


behavior Main {
    double in[][3] = {
        { 12, -51,   4},
        {  6, 167, -68},
        { -4,  24, -41},
        { -1, 1, 0},
        { 2, 0, 3},
    };

    int main()
    {
        mat R, Q;
        mat x = matrix_copy(3, in, 5);
        householder(x, &R, &Q);

        puts("Q"); matrix_show(Q);
        puts("R"); matrix_show(R);

        // to show their product is the input matrix
        mat m = matrix_mul(Q, R);
        puts("Q * R"); matrix_show(m);

        return 0;
    }
}