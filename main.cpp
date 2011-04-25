
#include <string>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <TRandom3.h>

using std::string;
using std::cout;
using std::cerr;
using std::endl;

typedef double field_type;
typedef field_type **field;
typedef field_type  *vector;
field malloc_field(int L);
inline vector malloc_vector(int L) {
    return (vector)malloc(L * sizeof(field_type));
}
void free_vector(vector vektor);
void free_field(field feld);
vector field_vector(field feld, vector vektor, int L);
double dot_product(vector u, vector v, int L);
vector cg(int L, field A, vector x, vector b, vector (func)(field, vector, int), double tolerance, int max_iterations, bool flag);
void print_vector(vector v, int L);
void print_matrix(field f, int L);

const int L = 2;

int main(int argc, char *argv[]) {

    field feld;
    feld = malloc_field(L);

    // assign random numbers between 0 and 1 to matrix feld
    const int L2 = L*L;
    TRandom3 *ran = new TRandom3();
    for (int i = 0; i < L2; i++) {
        feld[0][i] = ran->Uniform();
    }
    feld[0][0] = 0.211322;
    feld[0][1] = 0.876808;
    feld[1][0] = 0.931371;
    feld[1][1] = 0.463767;

    // make matrix feld non-singular: M = (M+M^t)/2
    field feld2 = malloc_field(L);
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            feld2[i][j] = 0;
            for(int k = 0; k < L; k++)
                feld2[i][j] += (feld[i][k]*feld[j][k]);
        }
    }

    vector x = malloc_vector(L);
    vector b = malloc_vector(L);
    for (int i = 0; i < L; i++) {
        x[i] = ran->Uniform();
        b[i] = ran->Uniform();
    }
    delete ran;

    b[0] = 0.716600;
    b[1] = 0.086319;

    /*
    print_vector(vektor, L);
    print_matrix(feld, L);

    vector result = field_vector(feld, vektor);
    print_vector(result, L);

    free_vector(result);
    */

    int max_iterations = L*3;
    double tolerance = 1E-10;
    vector result = NULL;
    result = cg(L, feld2, x,      b, field_vector, tolerance, max_iterations, false);
    print_vector(result, L);
    result = cg(L, feld2, result, b, field_vector, tolerance, max_iterations, true);
    print_vector(result, L);

    free_vector(b);
    free_field(feld);
    free_field(feld2);
    free_vector(x);
    return 0;
}

void free_field(field feld) {
    free_vector(feld[0]);
    free(feld);
    feld = NULL;
}
void free_vector(vector vektor) {
    free(vektor);
    vektor = NULL;
}

vector cg(int L, field A, vector x, vector b, vector (func)(field, vector, int), double tolerance, int max_iterations, bool flag) {
    vector r = malloc_vector(L);
    vector p = malloc_vector(L);
    vector s = NULL;
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
            free_vector(s);
            s = field_vector(A, p, L);
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
        cerr << "error in cg: no convergence" << endl;
    }

    free_vector(p);
    free_vector(r);
    free_vector(s);

    return x;
}

vector field_vector(field feld, vector vektor, int L) {
    vector result = malloc_vector(L);

    for (int i = 0; i < L; i++) {
        result[i] = 0;
        for (int j = 0; j < L; j++) {
            result[i] += feld[i][j]*vektor[j];
        }
    }

    return result;
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

field malloc_field(int L)
{
    field feld = (field)malloc(L * sizeof(field_type*));
    feld[0] = malloc_vector(L*L);
    if (feld[0] == NULL) {
        cerr << "Not enough memory for allocating field" << endl;
        exit(1);
    }
    for (int i = 1; i < L; i++)
        feld[i] = feld[0] + L * i;
    return feld;
}
