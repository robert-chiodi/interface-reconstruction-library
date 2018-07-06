include Makefile.in


default:
	@$(MAKE) library
	
all: 
	@$(MAKE)

opt:
	@$(MAKE) "FLAGS = $(OPTFLAGS)"

debug:
	@$(MAKE) "FLAGS = $(DBGFLAGS)"
	
test:
	@echo \*********************************************************
	@echo \*********************************************************
	@echo Make sure the current libs are compiled with gcov flags...
	@echo \*********************************************************
	@echo \********************************************************* 	
	@$(MAKE) "FLAGS = $(COVFLAGS)" "LDFLAGS += -lgcov -coverage" 
	cd $(TSTDIR) && $(MAKE) all

library:
	$(MAKE) -C src/

.PHONY: clean
clean:
	$(MAKE) clean -C $(SRCDIR)
	cd $(OBJDIR); rm -f *.gcov *.gcda *.gcno 
	cd $(TSTDIR) && $(MAKE) clean
