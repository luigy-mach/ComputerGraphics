
//      


#include <opencv2/imgproc/imgproc.hpp>

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <iostream>
#include <string>

#include "armadillo"

using namespace arma;
using namespace cv;
using namespace std;




/*template <typename T>
Mat<T> to_arma(const cv::Mat_<T> &src)
{
  Mat<T> dst(src.cols, src.rows);
  src.copyTo({src.rows, src.cols, dst.memptr()});
  return dst;
}
*/


template<class T3>
arma::Mat <T3> cvMat2armaMat(cv::Mat & cvMatIn) 
{ 
    return arma::Mat <T3> (cvMatIn.data, cvMatIn.rows, cvMatIn.cols,false); 
}

int main( int argc, char** argv )
{



    //cv::Mat greyMat, colorMat;
    //cv::cvtColor(colorMat, greyMat, CV_BGR2GRAY);

    String imageName( "/home/luigy/Desktop/topicos_ia/cnn/images.jpg" ); // by default
    if( argc > 1)
    {
        imageName = argv[1];
    }
    cv::Mat image, greyMat;
    image = imread( imageName, IMREAD_COLOR ); // Read the file
    if( image.empty() )                      // Check for invalid input
    {
        cout <<  "Could not open or find the image" << std::endl ;
        return -1;
    }
    cv::cvtColor(image, greyMat, CV_BGR2GRAY);
    
    


    for(int row = 0; row < greyMat.rows; ++row) {
        uchar* p = greyMat.ptr(row);
        for(int col = 0; col < greyMat.cols; ++col) {
             cout<<">: "<< row<<"+"<<col<<"-"<<*p++<<endl;  //points to each pixel value in turn assuming a CV_8UC1 greyscale image 

        }
    }

    //arma::mat arma_mat( reinterpret_cast<double*> greyMat.data, greyMat.rows, greyMat.cols )


    namedWindow( "Display window", WINDOW_AUTOSIZE ); // Create a window for display.
    //imshow( "Display window", image );                // Show our image inside it.
    imshow( "Display window", greyMat );                // Show our image inside it.
    waitKey(0); // Wait for a keystroke in the window
    return 0;
}