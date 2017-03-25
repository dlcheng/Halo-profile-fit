#------------------------------------------Always recommended 
OPT += -DLOG_SCALE         
OPT += -DIGNORE_SL_CUT                     #Ignore the cut of smoothing length
#------------------------------------------Analyze velocity profile
#OPT += -DVEL_PROFILE
#------------------------------------------If decay dark matter
#OPT += -DDECAY_DARK_MATTER              
#------------------------------------------Fitting profiles
OPT += -DGENERAL_NFW
#OPT += -DBURKERT
#OPT += -DBURKERT_NFW
#------------------------------------------Select target computer
SYSTYPE="ITSC"
#------------------------------------------Adjust settings for target computer

ifeq ($(SYSTYPE),"Workstation")
CC       =   gcc   
OPTIMIZE =   -O3 -Wall 
GSL_INCL =  -I/home/dlcheng/Install/gsl/include
GSL_LIBS =  -L/home/dlcheng/Install/gsl/lib
endif

ifeq ($(SYSTYPE),"dlcheng")
CC       =   gcc   
OPTIMIZE =   -O3 -Wall 
GSL_INCL =  -I/home/dalong/Install/gsl/include
GSL_LIBS =  -L/home/dalong/Install/gsl/lib
endif

ifeq ($(SYSTYPE),"ITSC")
CC       =   gcc   
OPTIMIZE =   -O3 -Wall
GSL_LIBS=   -L/users/s0902248/Lib/gsl-1.9/lib  -Wl,"-R /users/s0902248/Lib/gsl-1.9/lib"
GSL_INCL =  -I/users/s0902248/Lib/gsl-1.9/include
endif

OPTIONS =  $(OPTIMIZE) $(OPT) 

EXEC   = Halo_profile

OBJS   = allvars.o fit_data.o fit_model.o init.o main.o output.o process.o

INCL   = allvars.h proto.h define.h Makefile


CFLAGS = $(OPTIONS) $(GSL_INCL) 


LIBS   = $(GSL_LIBS) -lgsl -lgslcblas -lm 

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC) *.gch

