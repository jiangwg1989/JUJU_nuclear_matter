INCLUDE = -I/home/slime/software/gsl/include

cc          =    g++  ${INCLUDE}
cflag       =    -g  
LOCAL_LIBS  =    -L/home/slime/software/gsl/lib -lgsl -lgslcblas 
all_objects =    CCD_neutronmatter_emulator.o CCD_general_eigvalue.o

CCD_neutronmatter_emulator.exe : ${all_objects}    
	${cc} -o CCD_neutronmatter_emulator.exe ${all_objects}  ${LOCAL_LIBS}

CCD_neutronmatter_emulator.o  : CCD_neutronmatter_emulator.cpp
	${cc} ${cflag} -c  CCD_neutronmatter_emulator.cpp
#gauss_Laguerre.o:        gauss_Laguerre.cpp
#        ${cc} ${cflag} -c gauss_Laguerre.cpp
#Gauss_Legendre.o:       Gauss_Legendre.cpp
#        ${cc} ${cflag} -c Gauss_Legendre.cpp
#moshinsky.o:    moshinsky.cpp
#        ${cc} ${cflag} -c moshinsky.cpp
#Track_generation.o:     Track_generation.cpp
#        ${cc} ${cflag} -c Track_generation.cpp

#%.o: %.cpp
#       ${cc} ${cflag} -c $<
#
#%.o: %.c
#       ${cc} -c $<



clean:  
	rm *.o *.exe

