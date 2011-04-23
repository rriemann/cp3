
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
inline vector malloc_vector(int L) {return (vector)malloc(L * sizeof(field_type));}
vector field_vector(field feld, vector vektor);
double dot_product(vector u, vector v);
vector cg(int L, field A, vector x, vector b, vector (func)(field, vector), double tolerance, int max_iterations, bool flag);
void print_vector(vector v);
void print_matrix(field f);

const int L = 4;


int main(int argc, char *argv[]) {
  cout << "erster Test" << endl;

  // create matrix LxL
  field feld;
  feld = malloc_field(L);
  // assign random numbers between 0 and 1
  const int L2 = L*L;
  TRandom3 *ran = new TRandom3(0);
  for(int i = 0; i < L2; i++) {
    feld[0][i] = ran->Uniform();
  }

  field feld_non_singular = malloc_field(L);
  for(int i = 0; i < L; i++) {
    for(int j = 0; j < L; j++) {
      feld_non_singular[i][j] = (feld[i][j] + feld[j][i])/2.0;
    }
  }
  free(feld[0]);
  feld = feld_non_singular;
//   free(feld);
  // create vector
  vector vektor = malloc_vector(L);
  vector b = malloc_vector(L);
  for(int i = 0; i < L; i++) {
    vektor[i] = ran->Uniform();
    b[i] = ran->Uniform();
  }
  delete ran;

  print_vector(vektor);
  print_matrix(feld);

  vector result;
  result = field_vector(feld, vektor);
  print_vector(result);
  
  free(result);

  result = cg(L, feld, vektor, b, field_vector, 10E-8, 1000, false);
  print_vector(result);
  return 0;
}

vector cg(int L, field A, vector x, vector b, vector (func)(field, vector), double tolerance, int max_iterations, bool flag) {
  vector r0, x0;
  if(flag) {
    for(int i = 0; i < L; i++)
      x0[i] = 0;
    r0 = b;
  } else {
    r0 = func(A, x);
    for(int i = 0; i < L; i++) {
      r0[i] = b[i] - r0[i];
    }
  }

  double r02 = dot_product(r0,r0);
  if(r02 < tolerance)
    return r0;

  vector p0;
  for(int i = 0; i < L; i++) {
    p0[i] = r0[i];
  }

  vector s;
  vector x1 = malloc_vector(L);
  vector r1 = malloc_vector(L);
  double alpha;
  for(int iteration = 0; iteration < max_iterations; iteration++) {
    s = field_vector(A, p0);
    alpha = dot_product(p0,r0)/dot_product(p0,s);
    for(int i = 0; i < L; i++) {
      x1[i] = x0[i]+alpha*p0[i];
      r1[i] = r0[i]-alpha*s[i];
    }
    double r12 = dot_product(r1,r1);
    if(r12 < tolerance)
      return r1;

    double beta = r12/r02;
    for(int i = 0; i < L; i++) {
      p0[i] = r1[i]+beta*p0[i];
      r0[i] = r1[i];
    }
    free(s);
  }

  free(x1);
  free(r1);

  if(!flag) free(r0);
  return r1;
}

vector field_vector(field feld, vector vektor) {
  vector result = malloc_vector(L);

  for(int i = 0; i < L; i++) {
    result[i] = 0;
    for(int j = 0; j < L; j++) {
      result[i] += feld[i][j]*vektor[j];
    }
  }

  return result;
}

void print_vector(vector v) {
  for(int i = 0; i < L; i++) {
    cout << v[i] << endl;
  }
}

void print_matrix(field f) {
  for(int i = 0; i < L; i++) {
    for(int j = 0; j < L; j++) {
      cout << f[i][j] << " ";
    }
    cout << endl;
  }
}

double dot_product(vector u, vector v) {
  double sum = 0;
  for(int i = 0; i < L; i++) {
    sum += u[i]*v[i];
  }
  return sum;
}

field malloc_field(int L)
{
    field feld;
    feld = (field)malloc(L * sizeof(field_type*));
    feld[0] = (vector)malloc(L * L * sizeof(field_type));
    if(feld[0] == NULL) {
        cerr << "Not enough memory for allocating field" << endl;
        exit(1);
    }
    for(int x = 1; x < L; x++)
        feld[x] = feld[0] + L * x;
    return feld;
}
