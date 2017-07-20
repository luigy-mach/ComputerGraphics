

edge=./convolutionImgKernelEdgeDetection
gauss=./convolutionImgKernelGaussianBlur
random=./convolutionImgKernelRandom
sharpen=./convolutionImgKernelSharpen

num_procs=('lena.ascii.pgm' 'saturn.ascii.pgm' 'mario.pgm')

echo "TEST convolucion sobre imagenes PGM"

for i in "${num_procs[@]}"
do
  $edge $i
  $gauss $i
  $random $i
  $sharpen $i
done


