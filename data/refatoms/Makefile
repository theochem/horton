all: 001__H_N01_M2.h5

refatoms.tar.bz2:
	wget http://users.ugent.be/~tovrstra/horton/refatoms.tar.bz2

001__H_N01_M2.h5: refatoms.tar.bz2
	tar -xvjf refatoms.tar.bz2
	python ./fixnames.py
