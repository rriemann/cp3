
#include <string>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <TRandom3.h>
#include "geom_pbc.h"


using std::string;
using std::cout;
using std::cerr;
using std::endl;

typedef double field_type;
typedef field_type **field;
typedef field_type  *vector;
void free_vector(vector vektor);
void free_field(field feld);
vector field_vector(field feld, vector vektor, int L);
double dot_product(vector u, vector v, int L);
void cg(int L, vector x, vector b, void (func)(vector, double, vector), double tolerance, int max_iterations, bool flag);
void print_vector(vector v, int L);
void print_matrix(field f, int L);

const int L = 8;
int ndim = 3;
int *lsize;
int nvol, **nn;

int main(int argc, char *argv[]) {

    lsize = (int*)malloc(ndim*sizeof(int));
    for(int i = 0; i < ndim; i++)
      lsize[i] = L;
    geom_pbc(); // assign nn[][] and nvol

    field_type eta[nvol];
    field_type phi[nvol];

    // assign random numbers between 0 and 1
    TRandom3 *ran = new TRandom3();
    for (int i = 0; i < nvol; i++) {
        eta[i] = ran->Uniform();
    }
    delete ran;

    result = cg(nvol, x, b, laplace, tolerance, max_iterations, false);

    /*
    print_vector(vektor, L);
    print_matrix(feld, L);

    vector result = field_vector(feld, vektor);
    print_vector(result, L);

    free_vector(result);

    int max_iterations = L*3;
    double tolerance = 1E-10;
    vector result = NULL;
    result = cg(L, feld2, x,      b, field_vector, tolerance, max_iterations, false);
    print_vector(result, L);
    result = cg(L, feld2, result, b, field_vector, tolerance, max_iterations, true);
    print_vector(result, L);

    free_vector(b);
    free_vector(x);
    */

    free(lsize);
    return 0;
}

void laplace(vector x, double m2, vector result){
  for(int i = 0; i < nvol; ++i){
    result[i] = (2*ndim + m2)*x[i];
    for(int j = 1; j <= ndim; ++j){
      result[i] -= x[nn[j][i]] + x[nn[ndim+j][i]];
    }
  }
}

void cg(int nvol, vector x, vector b, void (func)(vector, double m2, vector), double tolerance, int max_iterations, bool flag) {
    field_type r[nvol], p[nvol], s[nvol];
    double bb = dot_product(b, b, L);
    double rr_tol = bb * tolerance*tolerance;
    double rr;

    if (flag) {
        s = func(A,x,L);
        rr = 0;
        for (int i = 0; i < L; i++)
            r[i] = b[i] - s[i];
        rr = dot_product(r, r, L);
    } else {
        for (int i = 0; i < L; i++) {
            x[i] = 0;
            r[i] = b[i];
        }
        rr = bb;
    }

    cout << "iter = " << 0 << "   rr = " << rr << endl;
    if (rr > rr_tol) {

        for(int i = 0; i < L; i++) {
            p[i] = r[i];
        }

        double alpha;
        for (int iteration = 1; iteration < max_iterations; iteration++) {
            func(p, m2, s);
            alpha = rr/dot_product(p, s, L);
            for (int i = 0; i < L; i++) {
                x[i] += alpha*p[i];
                r[i] -= alpha*s[i];
            }
            double rr_old = rr;
            rr = dot_product(r, r, L);
            cout << "iter = " << iteration << "   rr = " << rr << endl;
            if (rr <= tolerance)
                break;

            double beta = rr/rr_old;
            for (int i = 0; i < L; i++) {
                p[i] = r[i]+beta*p[i];
            }
        }
        if(rr > tolerance)
            cerr << "error in cg: no convergence" << endl;
    }
}

void print_vector(vector v, int L) {
    for (int i = 0; i < L; i++) {
        cout << v[i] << endl;
    }
}

void print_matrix(field f, int L) {
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            cout << f[i][j] << " ";
        }
        cout << endl;
    }
}

double dot_product(vector u, vector v, int L) {
    double sum = 0;
    for (int i = 0; i < L; i++) {
        sum += u[i]*v[i];
    }
    return sum;
}
