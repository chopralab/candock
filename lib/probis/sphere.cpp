#include "sphere.h"

void Sphere::print_sphere(ostream& os) {
    char* buffer = new char[100];
    sprintf(buffer, "%5d%5.2f%8.3f%8.3f%8.3f", num, r, crd.x, crd.y, crd.z);
    os << buffer;
    delete[] buffer;
}
