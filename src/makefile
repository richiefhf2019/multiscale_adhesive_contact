#===============================================================================
# makefile, by Houfu Fan, 1.26.2014
#===============================================================================

#SUFFIXES DECLARATION
.SUFFIXES:
.SUFFIXES: .f90 .o

#===============================Macro Definition================================
#Compile Information
FC=gfortran
#FLAGS=-O3 -W -Wall -c
FLAGS=-O3 -c -g
#Source Files
SRCS= preprocess.f90\
initialize.f90\
write_tecplot.f90\
compute_mass.f90\
compute_force.f90\
dynamic_solve.f90\
formkf.f90\
output.f90\
write_vtu.f90\
readaline.f90

MODULE_SRC =  module_config.f90\
module_constant.f90\
module_math.f90\
module_constitutive.f90\
module_fem.f90\
module_elem01.f90\
module_elem02.f90\
module_elem03.f90\
module_elem04.f90


MAIN_SRC= main.f90

#Objective Files
OBJS= $(SRCS:.f90=.o)
MODULE_OBJ = $(MODULE_SRC:.f90=.o)
MAIN_OBJ= main.o

TARGET=dynamic2d.x


#=================================Main Target===================================
$(TARGET): $(OBJS) $(MAIN_OBJ)
	$(FC) -o $(TARGET) $(MODULE_OBJ) $(OBJS) $(MAIN_OBJ)
	mv $(TARGET) ../inp
	date
#====================================Clean======================================
clean:
	rm -f *.o *.mod *.bak

#===============================Inference Rules=================================
.f90.o:
	$(FC) $(FLAGS) $<

#===========================Define the dependency===============================
$(MAIN_OBJ): $(OBJS) $(MAIN_SRC)
$(OBJS):$(MODULE_OBJ)



 