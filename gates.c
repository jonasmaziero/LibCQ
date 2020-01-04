#include <stdlib.h>
#include <math.h>

void id(int *d, double *idr) {
    int j,k;
    for(j = 0; j < (*d); j++){
        for(k = 0; k < (*d); k++){
            if(j == k){
                *(idr+j*(*d)+k) = 1.0;
            }
            else{
                *(idr+j*(*d)+k) = 0.0;
            }
        }
    }
}


void id_c(int *d, double _Complex *idc) {
    int j,k;
    for(j = 0; j < (*d); j++){
        for(k = 0; k < (*d); k++){
            if(j == k){
                *(idc+j*(*d)+k) = 1.0;
            }
            else{
                *(idc+j*(*d)+k) = 0.0;
            }
        }
    }
}

void zm(int *d, double *z) {
    int j,k;
    for(j = 0; j < (*d); j++){
        for(k = 0; k < (*d); k++){
            *(z+j*(*d)+k) = 0.0;
        }
    }
}

void zm_c(int *d, double _Complex *z) {
    int j,k;
    for(j = 0; j < (*d); j++){
        for(k = 0; k < (*d); k++){
            *(z+j*(*d)+k) = 0.0;
        }
    }
}

/*
int main(){
    double *A;
    int d;
    d = 2;
    A = (double *)malloc(d*d*sizeof(double));
    void array_display(int *, int *, double *);
    void id(int *, double *);
    id(&d, A);
    void zm(int *, double *);
    zm(&d, A);
    //array_display(&d, &d, A);
    double _Complex *B;
    B = (double _Complex *)malloc(d*d*sizeof(double _Complex));
    void id_c(int *, double _Complex *);
    id_c(&d, B);
    void zm_c(int *, double _Complex *);
    zm_c(&d, B);
    void array_display_c(int *, int *, double _Complex *);
    array_display_c(&d, &d, B);
    return 0;
}
*/
