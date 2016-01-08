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
ARGS=evert.pgm 10
OUT=out.m

$(EXECSEQ): $(EXECSEQ).o
	@$(CC) $^ -o $@

$(EXECSEQ).o: $(EXECSEQ).cc
	@$(CC) -c $<

$(EXEC1D): $(EXEC1D).o
	@$(CC) $^ -o $@

$(EXEC1D).o: $(EXEC1D).cc
	@$(CC) -c $<

$(EXEC): $(EXEC).o fftPar1D.o
	@$(CC) $^ -o $@

filtragePar.o: filtragePar.cc
	@$(CC) -c $<

run: $(EXEC) 
	@$(RUN) $(RFLAGS) $(PROC) $^ $(ARGS)

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

