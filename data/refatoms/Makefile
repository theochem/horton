all: 001__H_N01_M2.h5

refatoms.tar.gz:
	curl -O https://users.ugent.be/~tovrstra/horton/refatoms.tar.gz

001__H_N01_M2.h5: refatoms.tar.gz fixformat.py
	tar -xzf refatoms.tar.gz
	python ./fixformat.py
