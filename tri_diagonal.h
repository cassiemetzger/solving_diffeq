#pragma once

#include <cmath>
#include <iostream>
#include <vector>

template <typename T> class tri_diagonal {
public:
  tri_diagonal(const std::vector<T> a0, const std::vector<T> b0,
               const std::vector<T> c0) {
    a = a0;
    b = b0;
    c = c0;
    u.resize(a.size());
    l.resize(a.size());
    y.resize(a.size());
  }

  void compute_lu() {
    u[0] = b[0]; 
    for (int i = 1; i < a.size(); i++) {
      l[i] = a[i] / u[i - 1];
      u[i] = b[i] - l[i] * c[i - 1];
    }
  }

  std::vector<T> solve(const std::vector<T> &r) {
    if (r.size() != a.size()) {
      std::cout << "Error: r.size() != a.size()" << std::endl;
      return std::vector<T>();
    }
    //from notes 
    std::vector<T> x(r.size());
    //y0 = r0
    y[0] = r[0]; 
    for(int i = 1; i < r.size(); i++){
        y[i] = r[i] - l[i]*y[i-1];
        //yi = ri - li*yi-1
    }
    //xn-1 = yn-1/un-1 
    x[r.size()-1] = y[r.size()-1]/u[r.size()-1];
    //r.size = n 
    for(int i = r.size()-2; i >= 0; i--){
        x[i] = (y[i] - c[i]*x[i+1])/u[i];
        //xi = yi - ci*xi+1/un
    }

    return x;
  }

  std::vector<T> multiply(const std::vector<T> &x) {
    if (x.size() != a.size()) {
      std::cout << "Error: x.size() != a.size()" << std::endl;
      return std::vector<T>();
    }
    std::vector<T> r(x.size());
    //matrix multiplication of A*x 
    for(int i = 0; i < x.size(); i++){
        r[i] = b[i]*x[i];
        if(i > 0){
            r[i] += a[i]*x[i-1];
        }
        if(i < x.size()-1){
            r[i] += c[i]*x[i+1];
        }
    }
    return r;
  }

  std::vector<T> a, b, c;
  std::vector<T> l, u;
  std::vector<T> y;
};
