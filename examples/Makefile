random:
	@rm -f dumpdirectory/particle-*
	@mkdir -p _random
	@CC99='mpicc -std=c99' qcc -D_PTAG -D_MPI=1 -events -o a.out random.c -lm
	@cd _random/ && cp ../a.out . && mpirun -np 4 ./a.out
	@cd _random/ && gawk -f $(BASILISK)/gnuplot.awk < ../random.c > plots && sh $(BASILISK)/gnuplot.sh || rm -f plots

load:
	@rm -f dumpdirectory/particle-*
	@mkdir -p _load
	@CC99='mpicc -std=c99' qcc -D_PTAG -D_MPI=1 -events -o a.out load.c -lm
	@cd _load/ && rm -f time.dat && cp ../a.out . && for npe in 2 4 8 16; do mpirun -np $$npe ./a.out; done
	@cd _load/ && gawk -f $(BASILISK)/gnuplot.awk < ../load.c > plots && sh $(BASILISK)/gnuplot.sh || rm -f plots

ash:
	@rm -f dumpdirectory/particle-*
	@mkdir -p _ash
	@CC99='mpicc -std=c99' qcc -D_PINERTIAL -D_MPI=1 -events -o a.out ash.c -lm
	@cd _ash/ && cp ../a.out . && mpirun -np 4 ./a.out
	@cd _ash/ && gawk -f $(BASILISK)/gnuplot.awk < ../ash.c > plots && sh $(BASILISK)/gnuplot.sh || rm -f plots

confetti:
	@rm -f dumpdirectory/particle-*
	@mkdir -p _confetti
	@CC99='mpicc -std=c99' qcc -D_PTAG -D_MPI=1 -events -o a.out confetti.c -lm
	@cd _confetti/ && cp ../a.out . && mpirun -np 4 ./a.out
	@cd _confetti/ && gawk -f $(BASILISK)/gnuplot.awk < ../confetti.c > plots && sh $(BASILISK)/gnuplot.sh || rm -f plots
