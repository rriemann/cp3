
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

const int L = 4;

int main(int argc, char *argv[]) {

    field feld;
    feld = malloc_field(L);

    // assign random numbers between 0 and 1 to matrix feld
    const int L2 = L*L;
    TRandom3 *ran = new TRandom3(0);
    for (int i = 0; i < L2; i++) {
        feld[0][i] = ran->Uniform();
    }

    // make matrix feld non-singular: M = (M+M^t)/2
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            feld[i][j] = (feld[i][j] + feld[j][i])/2.0;
        }
    }

    vector vektor = malloc_vector(L);
    vector b = malloc_vector(L);
    for (int i = 0; i < L; i++) {
        vektor[i] = ran->Uniform();
        b[i] = ran->Uniform();
    }
    delete ran;

    /*
    print_vector(vektor, L);
    print_matrix(feld, L);

    vector result = field_vector(feld, vektor);
    print_vector(result, L);

    free_vector(result);
    */

    vector result = cg(L, feld, vektor, b, field_vector, 10E-8, 1000, false);
    print_vector(result, L);

    free_vector(b);
    free_field(feld);
    free_vector(vektor);
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
    vector r0 = malloc_vector(L);
    vector x0 = malloc_vector(L);
    if (flag) {
        for (int i = 0; i < L; i++)
            x0[i] = 0;
        r0 = b;
    } else {
        r0 = func(A, x, L);
        for (int i = 0; i < L; i++) {
            r0[i] = b[i] - r0[i];
        }
    }

    double r02 = dot_product(r0, r0, L);
    if (r02 < tolerance)
        return r0;

    vector p0 = malloc_vector(L);
    for (int i = 0; i < L; i++) {
        p0[i] = r0[i];
    }

    vector s;
    vector x1 = malloc_vector(L);
    vector r1 = malloc_vector(L);
    double alpha;
    for (int iteration = 0; iteration < max_iterations; iteration++) {
        s = field_vector(A, p0, L);
        alpha = dot_product(p0, r0, L)/dot_product(p0, s, L);
        for (int i = 0; i < L; i++) {
            x1[i] = x0[i]+alpha*p0[i];
            r1[i] = r0[i]-alpha*s[i];
        }
        double r12 = dot_product(r1, r1, L);
        if (r12 < tolerance)
            return r1;

        double beta = r12/r02;
        for (int i = 0; i < L; i++) {
            p0[i] = r1[i]+beta*p0[i];
            r0[i] = r1[i];
        }
        free_vector(s);
    }

    free_vector(x1);
    free_vector(r1);

    if (!flag) free_vector(r0);
    return r1;
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
