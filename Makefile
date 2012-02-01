###############################################################################
#
#  Makefile template for the sources
#
###############################################################################

# -----------------------------------------------------------------------------
#   Sources for all modules
# -----------------------------------------------------------------------------
BINNAME = pd_sdp
CPPSRC	= pd_sdp.cpp\
            TPM.cpp\
            Matrix.cpp\
            SPM.cpp\
            PHM.cpp\
            SUP.cpp\
            EIG.cpp\
	    DPM.cpp\
	    PPHM.cpp\
	    T2PM.cpp\
       GutMat.cpp\
       Math.cpp

OBJ	= $(CPPSRC:.cpp=.o)

# -----------------------------------------------------------------------------
#   These are the standard libraries, include paths and compiler settings
# -----------------------------------------------------------------------------

BRIGHT_ROOT= .

INCLUDE = -Iinclude

LIBS= -llapack -lblas

CC	= gcc
CXX	= g++

# -----------------------------------------------------------------------------
#   Compiler & Linker flags
# -----------------------------------------------------------------------------
CFLAGS	= $(INCLUDE) -g -Wall -O2
LDFLAGS	=  -g -Wall -O2


# =============================================================================
#   Targets & Rules
# =============================================================================
all:
	@echo
	@echo '  +++ Building $(BINNAME)...'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-DPQG"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built successfully!'; \
	   echo; \
	 fi

#------------------------------------------------------------------------------
#  Compile with only P and Q conditions activated
#------------------------------------------------------------------------------

PQ:
	@echo
	@echo '  +++ Building $(BINNAME) with P and Q conditions'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-DPQ"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with P and Q conditions successfully!'; \
	   echo; \
	 fi

PQG:
	@echo
	@echo '  +++ Building $(BINNAME) with P, Q and G conditions active'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-DPQG"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with P, Q and G conditions successfully!'; \
	   echo; \
	 fi

PQGT1:
	@echo
	@echo '  +++ Building $(BINNAME) with P, Q, G and T_1 conditions'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-DPQGT1"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with P, Q, G and T_1 conditions successfully!'; \
	   echo; \
	 fi

PQGT:
	@echo
	@echo '  +++ Building $(BINNAME) with P, Q, G, T_1 and T_2 conditions'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-DPQGT"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with P, Q, G, T_1 and T_2 conditions successfully!'; \
	   echo; \
	 fi

PQGT2P:
	@echo
	@echo '  +++ Building $(BINNAME) with P, Q, G and T_2'\'' conditions'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-DPQGT2P"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with P, Q, G and T_2'\'' conditions successfully!'; \
	   echo; \
	 fi

PQGTP:
	@echo
	@echo '  +++ Building $(BINNAME) with P, Q, G, T_1 and T_2'\'' conditions'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME) DEFS="-DPQGTP"
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built with P, Q, G, T_1 and T_2'\'' conditions successfully!'; \
	   echo; \
	 fi

# -----------------------------------------------------------------------------
#   The default way to compile all source modules
# -----------------------------------------------------------------------------
%.o:	%.for Makefile
	@echo; echo "Compiling $(@:.o=.for) ..."
	$(FF) $(FFLAGS) $(SFLAGS) -c $(@:.o=.for) -o $@

%.o:	%.c Makefile
	@echo; echo "Compiling $(@:.o=.c) ..."
	$(CC) $(CFLAGS) $(SFLAGS) -c $(@:.o=.c) -o $@

%.o:	%.cpp Makefile
	@echo; echo "Compiling $(@:.o=.cpp) ..."
	$(CXX) $(CFLAGS) $(SFLAGS) $(DEFS) -c $(@:.o=.cpp) -o $@


# -----------------------------------------------------------------------------
#   Link everything together
# -----------------------------------------------------------------------------
$(BRIGHT_ROOT)/$(BINNAME):	Makefile $(OBJ) 
	@echo; echo "Linker: creating $(BRIGHT_ROOT)/$(BINNAME) ..."
	$(CXX) $(LDFLAGS) $(SFLAGS) -o $(BRIGHT_ROOT)/$(BINNAME) $(OBJ) $(LIBS)

# -----------------------------------------------------------------------------
#   Create everything newly from scratch
# -----------------------------------------------------------------------------
new:	clean all

# -----------------------------------------------------------------------------
#   Clean up all object files
# -----------------------------------------------------------------------------
clean:
	@echo -n '  +++ Cleaning all object files ... '
	@echo -n $(OBJ)
	@rm -f $(OBJ)
	@echo 'Done.'

doc:
	doxygen doc-config



# ====================== End of file 'Makefile.in' ========================== #
