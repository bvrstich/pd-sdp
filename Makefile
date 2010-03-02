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
            Matrix.cpp\
            SPM.cpp\
            TPM.cpp\
            PHM.cpp\
            DPM.cpp\
            SUP.cpp\
            EIG.cpp

OBJ	= $(CPPSRC:.cpp=.o)

# -----------------------------------------------------------------------------
#   These are the standard libraries, include paths and compiler settings
# -----------------------------------------------------------------------------

BRIGHT_ROOT= .

INCLUDE = ./include

LIBS= -llapack -lblas

CC	= gcc
CXX	= g++

# -----------------------------------------------------------------------------
#   Compiler & Linker flags
# -----------------------------------------------------------------------------
CFLAGS	= -I$(INCLUDE) -g -Wall
LDFLAGS	= -g -Wall


# =============================================================================
#   Targets & Rules
# =============================================================================
all:
	@echo
	@echo '  +++ Building $(BINNAME)...'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME)
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

# -----------------------------------------------------------------------------
#   The default way to compile all source modules
# -----------------------------------------------------------------------------
%.o:	%.for Makefile
	@echo; echo "Compiling $(@:.o=.for) ..."
	$(FF) -c $(FFLAGS) $(SFLAGS) $(@:.o=.for) -o $@

%.o:	%.c Makefile
	@echo; echo "Compiling $(@:.o=.c) ..."
	$(CC) -c $(CFLAGS) $(SFLAGS) $(@:.o=.c) -o $@

%.o:	%.cpp Makefile
	@echo; echo "Compiling $(@:.o=.cpp) ..."
	$(CXX) -c $(CFLAGS) $(SFLAGS) $(DEFS) $(@:.o=.cpp) -o $@


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


# ====================== End of file 'Makefile.in' ========================== #
