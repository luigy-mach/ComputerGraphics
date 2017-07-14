all: programCimgAndCuda

programCimgAndCuda: convolution.o
	g++ -o programCimgAndCuda convolution.o imageCimg.cpp -O2 -L/usr/X11R6/lib -lm -lpthread -lX11  

convolution.o:
    nvcc -c convolution.cu 
