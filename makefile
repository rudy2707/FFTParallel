CC=mpiCC

RUN=mpirun
RFLAGS=-np

EXECSEQ=fftSeq
EXEC1D=fftPar1D
EXEC=fftParallel

RM=rm -f

PROC=16
SEED=19
DATA=8
OUT=out.m

$(EXECSEQ): fftSeq.o
	@$(CC) $^ -o $@

fftSeq.o: fftSeq.cc 
	@$(CC) -c $<

$(EXEC1D): fftPar1D.o
	@$(CC) $^ -o $@

fftPar1D.o: fftPar1D.cc
	@$(CC) -c $<

$(EXEC): fftParallel.o
	@$(CC) $^ -o $@

fftParallel.o: fftParallel.cc
	@$(CC) -c $<

run: $(EXEC)
	@echo "not implemented yet"
	#@$(RUN) $(RFLAGS) $(PROC) $< $(SEED) $(DATA)

seq: $(EXECSEQ)
	@$(RUN) $(RFLAGS) $(PROC) $< $(SEED) $(DATA)

one: $(EXEC1D)
	@$(RUN) $(RFLAGS) $(PROC) $< $(SEED) $(DATA)

test: $(EXEC1D)
	@$(RUN) $(RFLAGS) $(PROC) $< $(SEED) $(DATA) > $(OUT)


clean:
	@$(RM) *.o
	@$(RM) $(OUT)
	@$(RM) $(EXEC)
	@$(RM) $(EXEC1D)
	@$(RM) $(EXECSEQ)
	
