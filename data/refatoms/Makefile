all: 001__H_N01_M2.h5

refatoms.tar.bz2:
	curl -O http://users.ugent.be/~tovrstra/horton/refatoms.tar.bz2

001__H_N01_M2.h5: refatoms.tar.bz2 fixformat.py
	tar -xjf refatoms.tar.bz2
	python ./fixformat.py
