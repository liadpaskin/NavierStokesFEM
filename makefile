FFLAGS= -g -O0  -fbacktrace  -fbounds-check -ffpe-trap=zero -fcray-pointer -ffree-line-length-none #-fdefault-integer-8
TARGET = NavierStokes_Liad2019
MODULES = bin/debug/MODULES.o
TECIOLIB= tecplot/lib/libtecio.a
OBJECTS =   $(MODULES) \
			bin/debug/IOMNGR.o  \
			bin/debug/CONTRL.o  \
			bin/debug/ALLOCATIONS.o \
			bin/debug/MEMORY.o \
			bin/debug/INELEM.o \
			bin/debug/INMESH.o \
			bin/debug/JOINTS.o \
			bin/debug/MATSET.o \
			bin/debug/profil.o \
			bin/debug/aponta.o \
			bin/debug/grafos.o \
			bin/debug/sort.o \
			bin/debug/genrcm.o \
			bin/debug/splot.o \
			bin/debug/matrixtrans2.o \
			bin/debug/fnroot.o \
			bin/debug/rcm.o \
			bin/debug/rootls.o \
			bin/debug/degree.o \
			bin/debug/LOADS.o \
			bin/debug/RHSV_READ.o \
			bin/debug/RHSV_MAKE.o \
			bin/debug/RHSV_LOAD.o \
			bin/debug/NavierStokes3D_Solver.o \
			bin/debug/HEXA8_NAVIER_STOKES_3D_CONVEC.o \
			bin/debug/HEXA8_NAVIER_STOKES_3D_DIFUS.o \
			bin/debug/HEXA8_NAVIER_STOKES_3D_MASSA.o \
			bin/debug/HEXA8_PRESSURE.o \
			bin/debug/HEXA8_CFLcalc.o \
			bin/debug/HEXA27_NAVIER_STOKES_3D_CONVEC.o \
			bin/debug/HEXA27_NAVIER_STOKES_3D_DIFUS.o \
			bin/debug/HEXA27_NAVIER_STOKES_3D_MASSA.o \
			bin/debug/HEXA27_PRESSURE.o \
			bin/debug/HEXA27_CFLcalc.o \
			bin/debug/TRIA6_NAVIER_STOKES_2D_CONVEC.o \
			bin/debug/TRIA6_NAVIER_STOKES_2D_DIFUS.o \
			bin/debug/TRIA6_NAVIER_STOKES_2D_MASSA.o \
			bin/debug/TRIA6_PRESSURE.o \
			bin/debug/TRIA6_CFLcalc.o \
			bin/debug/QUAD4_NAVIER_STOKES_2D_CONVEC.o \
			bin/debug/QUAD4_NAVIER_STOKES_2D_DIFUS.o \
			bin/debug/QUAD4_NAVIER_STOKES_2D_MASSA.o \
			bin/debug/QUAD4_PRESSURE.o \
			bin/debug/QUAD4_CFLcalc.o \
			bin/debug/SOLVER.o \
			bin/debug/getNALHS.o \
			bin/debug/CROUTFACT.o \
			bin/debug/oneStepLinearImplicit.o \
			bin/debug/matvec.o \
			bin/debug/local.o \
			bin/debug/ddot.o \
			bin/debug/addmatvec.o \
			bin/debug/FUNC.o \
			bin/debug/OUTENSIGHT.o \
			bin/debug/OUTTECPLT2D.o \
			bin/debug/OUTTECPLT3D.o \
			bin/debug/tec_vec_read.o \
			bin/debug/plotline2D.o \
			bin/debug/plotline3D.o \
			bin/debug/initialInterpolate.o \
			tecplot/include/TECIO.h
			
EXTRAINCLUDES=-I tecplot/include
LINKLIBS=-lpthread

makefile:
all: bin/debug/$(TARGET)

clean:
	-rm bin/debug/*.o bin/debug/*.mod bin/debug/$(TARGET)  #$(OBJECTS)
	
#NavierStokes_Liad2019 : 
#	gfortran src/NavierStokes_Liad2019.f90 -o NavierStokes_Liad2019.exe
	
bin/debug/$(TARGET) : $(OBJECTS)
	gfortran $(FFLAGS) src/$(TARGET).f90 $(OBJECTS) -I bin/debug $(EXTRAINCLUDES) $(TECIOLIB) -o bin/debug/$(TARGET)  $(LINKLIBS)  -lstdc++ 

bin/debug/ALLOCATIONS.o : $(MODULES) bin/debug/MEMORY.o
	gfortran $(FFLAGS) -c src/ALLOCATIONS/ALLOCATIONS.f90 -I bin/debug -o bin/debug/ALLOCATIONS.o
	
bin/debug/CONTRL.o : $(MODULES)
	gfortran $(FFLAGS) -c src/CONTRL/CONTRL.f90 -I bin/debug -o bin/debug/CONTRL.o
	
bin/debug/IOMNGR.o : $(MODULES)
	gfortran $(FFLAGS) -c src/IOMNGR/IOMNGR.f90 -I bin/debug -o bin/debug/IOMNGR.o
			
bin/debug/INMESH.o : $(MODULES)
	gfortran $(FFLAGS) -c src/INPUTMESH/INMESH.f90 -I bin/debug -o bin/debug/INMESH.o
	
bin/debug/INELEM.o : $(MODULES)
	gfortran $(FFLAGS) -c src/INPUTMESH/INELEM.f90 -I bin/debug -o bin/debug/INELEM.o
	
bin/debug/JOINTS.o : $(MODULES)
	gfortran $(FFLAGS) -c src/INPUTMESH/JOINTS.f90 -I bin/debug -o bin/debug/JOINTS.o
	
bin/debug/MATSET.o : $(MODULES)
	gfortran $(FFLAGS) -c src/INPUTMESH/MATSET.f90 -I bin/debug -o bin/debug/MATSET.o
	
bin/debug/LOADS.o : $(MODULES)
	gfortran $(FFLAGS) -c src/LOADS/LOADS.f90 -I bin/debug -o bin/debug/LOADS.o
	
bin/debug/oneStepLinearImplicit.o : $(MODULES)
	gfortran $(FFLAGS) -c src/SOLVER/oneStepLinearImplicit.f90 -I bin/debug -o bin/debug/oneStepLinearImplicit.o
bin/debug/getNALHS.o :
	gfortran $(FFLAGS) -c src/CROUTFACT/getNALHS.f90 -I bin/debug -o bin/debug/getNALHS.o
bin/debug/CROUTFACT.o :
	gfortran $(FFLAGS) -c src/CROUTFACT/CROUTFACT.f90 -I bin/debug -o bin/debug/CROUTFACT.o
	
bin/debug/HEXA8_CFLcalc.o  :
	gfortran $(FFLAGS) -c src/SOLVER/HEXA8/HEXA8_CFLcalc.f90 -I bin/debug -o bin/debug/HEXA8_CFLcalc.o
	
bin/debug/HEXA27_CFLcalc.o  :
	gfortran $(FFLAGS) -c src/SOLVER/HEXA27/HEXA27_CFLcalc.f90 -I bin/debug -o bin/debug/HEXA27_CFLcalc.o
	
bin/debug/TRIA6_CFLcalc.o  :
	gfortran $(FFLAGS) -c src/SOLVER/TRIA6/TRIA6_CFLcalc.f90 -I bin/debug -o bin/debug/TRIA6_CFLcalc.o
	
bin/debug/QUAD4_CFLcalc.o  :
	gfortran $(FFLAGS) -c src/SOLVER/QUAD4/QUAD4_CFLcalc.f90 -I bin/debug -o bin/debug/QUAD4_CFLcalc.o
	
bin/debug/tec_vec_read.o  :
	gfortran $(FFLAGS) -c src/OUTPUT/tec_vec_read.f90 -I bin/debug -o bin/debug/tec_vec_read.o
	
bin/debug/plotline2D.o  :
	gfortran $(FFLAGS) -c src/OUTPUT/plotline2D.f90 -I bin/debug -o bin/debug/plotline2D.o
bin/debug/plotline3D.o  :
	gfortran $(FFLAGS) -c src/OUTPUT/plotline3D.f90 -I bin/debug -o bin/debug/plotline3D.o
	
	
bin/debug/matvec.o : bin/debug/local.o
	gfortran $(FFLAGS) -c src/UTILITIES/matvec.f90 -I bin/debug -o bin/debug/matvec.o
bin/debug/local.o :
	gfortran $(FFLAGS) -c src/UTILITIES/local.f90 -I bin/debug -o bin/debug/local.o
bin/debug/ddot.o :
	gfortran $(FFLAGS) -c src/UTILITIES/ddot.f90 -I bin/debug -o bin/debug/ddot.o
bin/debug/addmatvec.o :
	gfortran $(FFLAGS) -c src/UTILITIES/addmatvec.f90 -I bin/debug -o bin/debug/addmatvec.o
bin/debug/FUNC.o :
	gfortran $(FFLAGS) -c src/UTILITIES/FUNC.f90 -I bin/debug -o bin/debug/FUNC.o
bin/debug/OUTENSIGHT.o :
	gfortran $(FFLAGS) -c src/OUTPUT/OUTENSIGHT.f90 -I bin/debug -o bin/debug/OUTENSIGHT.o
bin/debug/OUTTECPLT2D.o :
	gfortran $(FFLAGS) -c src/OUTPUT/OUTTECPLT2D.f90  $(FFLAGS) -I bin/debug $(EXTRAINCLUDES) $(TECIOLIB) -o bin/debug/OUTTECPLT2D.o $(LINKLIBS)  -lstdc++ 
bin/debug/OUTTECPLT3D.o :
	gfortran $(FFLAGS) -c src/OUTPUT/OUTTECPLT3D.f90  $(FFLAGS) -I bin/debug $(EXTRAINCLUDES) $(TECIOLIB) -o bin/debug/OUTTECPLT3D.o $(LINKLIBS)  -lstdc++ 
	
bin/debug/MODULES.o :
	gfortran $(FFLAGS) -c src/MODULES/MODULES.f90 -J bin/debug -o bin/debug/MODULES.o
	
bin/debug/MEMORY.o :
	gfortran $(FFLAGS) -c src/UTILITIES/MEMORY.f90 -I bin/debug -o bin/debug/MEMORY.o
	
bin/debug/SOLVER.o : $(MODULES) bin/debug/NavierStokes3D_Solver.o
	gfortran $(FFLAGS) -c src/SOLVER/SOLVER.f90 -I bin/debug -o bin/debug/SOLVER.o

bin/debug/NavierStokes3D_Solver.o : $(MODULES) bin/debug/getNALHS.o
	gfortran $(FFLAGS) -c src/SOLVER/NavierStokes3D_Solver.f90 -I bin/debug -o bin/debug/NavierStokes3D_Solver.o

bin/debug/HEXA8_NAVIER_STOKES_3D_CONVEC.o :
	gfortran $(FFLAGS) -c src/SOLVER/HEXA8/HEXA8_NAVIER_STOKES_3D_CONVEC.f90 -I bin/debug -o bin/debug/HEXA8_NAVIER_STOKES_3D_CONVEC.o
bin/debug/HEXA8_NAVIER_STOKES_3D_DIFUS.o : 
	gfortran $(FFLAGS) -c src/SOLVER/HEXA8/HEXA8_NAVIER_STOKES_3D_DIFUS.f90 -I bin/debug -o bin/debug/HEXA8_NAVIER_STOKES_3D_DIFUS.o
bin/debug/HEXA8_NAVIER_STOKES_3D_MASSA.o :
	gfortran $(FFLAGS) -c src/SOLVER/HEXA8/HEXA8_NAVIER_STOKES_3D_MASSA.f90 -I bin/debug -o bin/debug/HEXA8_NAVIER_STOKES_3D_MASSA.o
bin/debug/HEXA8_PRESSURE.o : $(MODULES)
	gfortran $(FFLAGS) -c src/SOLVER/HEXA8/HEXA8_PRESSURE.f90 -I bin/debug -o bin/debug/HEXA8_PRESSURE.o
	
bin/debug/HEXA27_NAVIER_STOKES_3D_CONVEC.o :
	gfortran $(FFLAGS) -c src/SOLVER/HEXA27/HEXA27_NAVIER_STOKES_3D_CONVEC.f90 -I bin/debug -o bin/debug/HEXA27_NAVIER_STOKES_3D_CONVEC.o
bin/debug/HEXA27_NAVIER_STOKES_3D_DIFUS.o : 
	gfortran $(FFLAGS) -c src/SOLVER/HEXA27/HEXA27_NAVIER_STOKES_3D_DIFUS.f90 -I bin/debug -o bin/debug/HEXA27_NAVIER_STOKES_3D_DIFUS.o
bin/debug/HEXA27_NAVIER_STOKES_3D_MASSA.o :
	gfortran $(FFLAGS) -c src/SOLVER/HEXA27/HEXA27_NAVIER_STOKES_3D_MASSA.f90 -I bin/debug -o bin/debug/HEXA27_NAVIER_STOKES_3D_MASSA.o
bin/debug/HEXA27_PRESSURE.o : $(MODULES)
	gfortran $(FFLAGS) -c src/SOLVER/HEXA27/HEXA27_PRESSURE.f90 -I bin/debug -o bin/debug/HEXA27_PRESSURE.o
	
bin/debug/TRIA6_NAVIER_STOKES_2D_CONVEC.o :
	gfortran $(FFLAGS) -c src/SOLVER/TRIA6/TRIA6_NAVIER_STOKES_2D_CONVEC.f90 -I bin/debug -o bin/debug/TRIA6_NAVIER_STOKES_2D_CONVEC.o
bin/debug/TRIA6_NAVIER_STOKES_2D_DIFUS.o : 
	gfortran $(FFLAGS) -c src/SOLVER/TRIA6/TRIA6_NAVIER_STOKES_2D_DIFUS.f90 -I bin/debug -o bin/debug/TRIA6_NAVIER_STOKES_2D_DIFUS.o
bin/debug/TRIA6_NAVIER_STOKES_2D_MASSA.o :
	gfortran $(FFLAGS) -c src/SOLVER/TRIA6/TRIA6_NAVIER_STOKES_2D_MASSA.f90 -I bin/debug -o bin/debug/TRIA6_NAVIER_STOKES_2D_MASSA.o
bin/debug/TRIA6_PRESSURE.o : $(MODULES)
	gfortran $(FFLAGS) -c src/SOLVER/TRIA6/TRIA6_PRESSURE.f90 -I bin/debug -o bin/debug/TRIA6_PRESSURE.o
	
bin/debug/QUAD4_NAVIER_STOKES_2D_CONVEC.o :
	gfortran $(FFLAGS) -c src/SOLVER/QUAD4/QUAD4_NAVIER_STOKES_2D_CONVEC.f90 -I bin/debug -o bin/debug/QUAD4_NAVIER_STOKES_2D_CONVEC.o
bin/debug/QUAD4_NAVIER_STOKES_2D_DIFUS.o : 
	gfortran $(FFLAGS) -c src/SOLVER/QUAD4/QUAD4_NAVIER_STOKES_2D_DIFUS.f90 -I bin/debug -o bin/debug/QUAD4_NAVIER_STOKES_2D_DIFUS.o
bin/debug/QUAD4_NAVIER_STOKES_2D_MASSA.o :
	gfortran $(FFLAGS) -c src/SOLVER/QUAD4/QUAD4_NAVIER_STOKES_2D_MASSA.f90 -I bin/debug -o bin/debug/QUAD4_NAVIER_STOKES_2D_MASSA.o
bin/debug/QUAD4_PRESSURE.o : $(MODULES)
	gfortran $(FFLAGS) -c src/SOLVER/QUAD4/QUAD4_PRESSURE.f90 -I bin/debug -o bin/debug/QUAD4_PRESSURE.o
	
bin/debug/matrixtrans2.o : 
	gfortran $(FFLAGS) -c src/INPUTMESH/CuthillMckee/matrixtrans2.f90 -I bin/debug -o bin/debug/matrixtrans2.o
bin/debug/profil.o : 
	gfortran $(FFLAGS) -c src/INPUTMESH/CuthillMckee/profil.f90 -I bin/debug -o bin/debug/profil.o
bin/debug/aponta.o : 
	gfortran $(FFLAGS) -c src/INPUTMESH/CuthillMckee/aponta.f90 -I bin/debug -o bin/debug/aponta.o
bin/debug/grafos.o : 
	gfortran $(FFLAGS) -c src/INPUTMESH/CuthillMckee/grafos.f90 -I bin/debug -o bin/debug/grafos.o
bin/debug/sort.o : 
	gfortran $(FFLAGS) -c src/INPUTMESH/CuthillMckee/sort.f90 -I bin/debug -o bin/debug/sort.o
bin/debug/genrcm.o : bin/debug/rcm.o bin/debug/fnroot.o
	gfortran $(FFLAGS) -c src/INPUTMESH/CuthillMckee/genrcm.f90 -I bin/debug -o bin/debug/genrcm.o	
bin/debug/splot.o : 
	gfortran $(FFLAGS) -c src/INPUTMESH/CuthillMckee/splot.f90 -I bin/debug -o bin/debug/splot.o
bin/debug/fnroot.o : bin/debug/rootls.o
	gfortran $(FFLAGS) -c src/INPUTMESH/CuthillMckee/fnroot.f90 -I bin/debug -o bin/debug/fnroot.o
bin/debug/rcm.o : bin/debug/degree.o
	gfortran $(FFLAGS) -c src/INPUTMESH/CuthillMckee/rcm.f90 -I bin/debug -o bin/debug/rcm.o
bin/debug/rootls.o : 
	gfortran $(FFLAGS) -c src/INPUTMESH/CuthillMckee/rootls.f90 -I bin/debug -o bin/debug/rootls.o
bin/debug/degree.o : 
	gfortran $(FFLAGS) -c src/INPUTMESH/CuthillMckee/degree.f90 -I bin/debug -o bin/debug/degree.o
	
bin/debug/RHSV_READ.o : 
	gfortran $(FFLAGS) -c src/LOADS/RHSV_READ.f90 -I bin/debug -o bin/debug/RHSV_READ.o
bin/debug/RHSV_MAKE.o : 
	gfortran $(FFLAGS) -c src/LOADS/RHSV_MAKE.f90 -I bin/debug -o bin/debug/RHSV_MAKE.o
bin/debug/RHSV_LOAD.o : 
	gfortran $(FFLAGS) -c src/LOADS/RHSV_LOAD.f90 -I bin/debug -o bin/debug/RHSV_LOAD.o
	
bin/debug/initialInterpolate.o : 
	gfortran $(FFLAGS) -c src/UTILITIES/initialInterpolate.f90 -I bin/debug -o bin/debug/initialInterpolate.o
