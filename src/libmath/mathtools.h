#define PI 3.141593
#define inv_ln_10 1./log(10)
#define ln_10 log(10)
#define D_H 3000 

double average(double*, int);
double maximum(double*, int);
double minimum(double*, int);
int int_maximum(int*, int);

double* shellsort(double* , int);
int* int_shellsort(int* , int);

int *generate_random_subset(int, int, int*);

double* invert_array(double*, int);

double* log_stepper(double, double, int);
double* lin_stepper(double, double, int);

void lin_bin(double*, double*, int, int, int*);
void average_bin(double*, double*, double *, double*, double*, int, int);

double get_interpolated_value(double*, double*, int, double);

double solid_angle(double, void*);
double integrate_solid_angle(double,double,double,double);
