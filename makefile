CC=mpiCC

RUN=mpirun
RFLAGS=-np

EXECSEQ=fftSeq
EXEC1D=fftPar1D
EXEC=filtragePar

RM=rm -f

PROC=4
SEED=19
DATA=4
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

filtragePar.o: filtragePar.cc fftPar1D.o
	@$(CC) -c $<

run: $(EXEC)
	@$(RUN) $(RFLAGS) $(PROC) $< $(SEED) $(DATA)

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

