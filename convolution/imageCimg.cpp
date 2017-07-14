
// g++ image_cimg.cpp -O2 -L/usr/X11R6/lib -lm -lpthread -lX11



#include "CImg.h"
#include "convolution.h"
using namespace cimg_library;

int main() {

  CImg<double> image("img/Armadillo.png");

  CImgDisplay main_disp(image,"imagen"); 





  while (!main_disp.is_closed() ) {
    main_disp.wait();
  }

  return 0;
}