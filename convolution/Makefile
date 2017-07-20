all: convolutionImgKernelRandom convolutionImgKernelEdgeDetection \
		convolutionImgKernelGaussianBlur convolutionImgKernelSharpen \
		convolutionMatrixComplete convolutionMatrixSinBordes


convolutionImgKernelRandom:
	nvcc -lcuda  convolutionImgKernelRandom.cu -arch=compute_20 -o convolutionImgKernelRandom

convolutionImgKernelEdgeDetection:
	nvcc -lcuda  convolutionImgKernelEdgeDetection.cu -arch=compute_20 -o convolutionImgKernelEdgeDetection

convolutionImgKernelGaussianBlur:
	nvcc -lcuda  convolutionImgKernelGaussianBlur.cu -arch=compute_20 -o convolutionImgKernelGaussianBlur

convolutionImgKernelSharpen:
	nvcc -lcuda  convolutionImgKernelSharpen.cu -arch=compute_20 -o convolutionImgKernelSharpen

convolutionMatrixComplete:
	nvcc -lcuda  convolutionMatrixComplete.cu -arch=compute_20 -o convolutionMatrixComplete

convolutionMatrixSinBordes:
	nvcc -lcuda  convolutionMatrixSinBordes.cu -arch=compute_20 -o convolutionMatrixSinBordes

