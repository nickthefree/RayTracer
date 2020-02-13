#include <iostream>
#include <fstream>
#include <cmath>

class myVec {
    public:
  double x,y,z;
  myVec(double x, double y, double z) : x(x), y(y), z(z) {}
  myVec operator + (const myVec& v) const {
      return myVec(x+v.x, y+v.y, z+v.z);
  }
  myVec operator - (const myVec& v) const {
      return myVec(x-v.x, y-v.y, z-v.z);
  }
  myVec operator * (double d) const {
      return myVec(x*d, y*d, z*d);
  }
  myVec operator / (double d) const {
      return myVec(x/d, y/d, z/d);
  }
  myVec normalize() const {
    double mg = sqrt(x*x + y*y + z*z);
    return myVec(x/mg,y/mg,z/mg);
  }
};
double dot(const myVec& a, const myVec& b) {
  return (a.x*b.x + a.y*b.y + a.z*b.z);
}

class Ray {
    public:
  myVec o,d;
  Ray(const myVec& o, const myVec& d) : o(o), d(d) {}
};

void RGB(myVec& col) {
  col.x = (col.x > 255) ? 255 : (col.x < 0) ? 0 : col.x;
  col.y = (col.y > 255) ? 255 : (col.y < 0) ? 0 : col.y;
  col.z = (col.z > 255) ? 255 : (col.z < 0) ? 0 : col.z;
}

class Sphere {
    public:
  myVec c;
  double r;
  Sphere(const myVec& c, double r) : c(c), r(r) {}
  myVec getNormal(const myVec& pi) const { return (pi - c) / r; }
  bool intersect(const Ray& ray, double &t) const {
    const myVec o = ray.o;
    const myVec d = ray.d;
    const myVec oc = o - c;
    const double b = 2 * dot(oc, d);
    const double c = dot(oc, oc) - r*r;
    double disc = b*b - 4 * c;
    if (disc < .0001) return false;
    disc = sqrt(disc);
    const double t0 = -b - disc;
    const double t1 = -b + disc;
    t = (t0 < t1) ? t0 : t1;
    return true;
  }
};



int main() {

  //sphere
  const myVec black(0, 0, 0);
  const myVec white(255, 255, 255);
  const myVec blue(0,0,255);
  const myVec red(255, 0, 0);
  const myVec green(0,255,0);

  const int H = 500;
  const int W = 500;



  const Sphere sphere(myVec(W*0.5, H*0.5, 50), 100);
  const Sphere light(myVec(0, 0, 50), 1);

  std::ofstream out("picture.ppm");
  out << "P3\n" << W << ' ' << H << ' ' << "255\n";

  double t;
  myVec pix_col(black);

  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      pix_col = black;

      const Ray ray(myVec(x,y,0),myVec(0,0,1));
      if (sphere.intersect(ray, t)) {
        const myVec pi = ray.o + ray.d*t;
        const myVec L = light.c - pi;
        const myVec N = sphere.getNormal(pi);
        const double dt = dot(L.normalize(), N.normalize());

        pix_col = (blue + white*dt) * 0.5;
        RGB(pix_col);
      }
      out << (int)pix_col.x << ' '
          << (int)pix_col.y << ' '
          << (int)pix_col.z << '\n';
    }
  }
  //triangle
  const myVec p0(100,100,50);
  const myVec p1(200,100,50);
  const myVec p2(100,200,50);
  const myVec e0 = p1 - p0;
  const myVec e1 = p2 - p1;
  const myVec norm(0,0,-1);
  double triK = dot(p0, norm);
  // R = Camera position in world coordinates
  // P = pixel on screen world coordinates
  // D = R - P
  for(int i = 0; i < 100; i++){
    for(int j = 0; j < 100; j++){
      //for(int vertices = 0; vertices < 3; vertices++){
        pix_col = black;
        const Ray ray(myVec(i+100,j+100,0),myVec(0,0,1));

        if(i + 100 >= 100 && i + 100 <= 200){
          if(j + 100 >= 100 & j + 100 <= 200){
              myVec myP = ray.d;
              myVec myR = ray.o;
              myVec myD = myR - myP;
              double myT = (triK - dot(myP, norm)/(dot(myD,norm)));
              pix_col = (green + white*myT) * 0.5;
              RGB(pix_col);
      }
      out << (int)pix_col.x << ' '
          << (int)pix_col.y << ' '
          << (int)pix_col.z << '\n';
        }
      //}
    }
  }
  std::cout << triK << std::endl;


}
