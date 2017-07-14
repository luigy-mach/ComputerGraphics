# Luigy Machaca Arcana
# Computer science - Arequipa, Perú  2017


nvcc matrix_mult_memory_shared.cu -o app 

bin=./app

#num_procs=(500 1000 1500 2000 2500 3000)
num_procs=(	#32 64 128 256 384 512 640 768 896  1024 1280 1536 1792
			#32 64 128 256 512 1024 1280 1536 1792 
			#32 64 128 256 512 1024 1280 1536 
			#1100 1200 1300 1400 1500 1600 1700 1800 1900 2000
			2000 400 6000 8000 10000 
			)

for i in "${num_procs[@]}"
do
  mult[$i]=`$bin $i`	
  #echo $i ,${mult[$i]},${mult[$i]}
  echo $i ,${mult[$i]}
done


