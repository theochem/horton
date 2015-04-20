default:
	@echo 'Choose one of the following targets'
	@echo '* libint: downloads and compiles libint-2.0.3-stable.tgz for mpqc'
	@echo '* libint-linux-x86_64-gcc-4.7.2: downloads and unpacks a binary libint for linux-x86_64 compiled with gcc-4.7.2'
	@echo '* libxc: downloads and compiles libxc-2.1.0'

libint: libint-2.0.3-stable/lib/libint2.a

libxc: libxc-2.1.0/src/.libs/libxc.a


# Build libint from scratch

libint-2.0.3-stable/lib/libint2.a: libint-2.0.3-stable
	cd libint-2.0.3-stable && CONFIG_SHELL=/bin/sh CFLAGS='-fPIC' CPPFLAGS='-fPIC' ./configure --with-cxx-optflags='-O1'
	cd libint-2.0.3-stable && $(MAKE)

libint-2.0.3-stable: libint-2.0.3-stable.tgz
	tar -xzvf libint-2.0.3-stable.tgz

libint-2.0.3-stable.tgz:
	curl -OL 'http://downloads.sourceforge.net/project/libint/libint-for-mpqc/libint-2.0.3-stable.tgz'


# Binary compilation (mainly for Travis CI)

libint-2.0.3-bin-horton-linux-x86_64-gcc-4.7.2.tar.bz2:
	curl -O http://users.ugent.be/~tovrstra/horton/libint-2.0.3-bin-horton-linux-x86_64-gcc-4.7.2.tar.bz2

libint-linux-x86_64-gcc-4.7.2: libint-2.0.3-bin-horton-linux-x86_64-gcc-4.7.2.tar.bz2
	tar -xjf libint-2.0.3-bin-horton-linux-x86_64-gcc-4.7.2.tar.bz2


# Build libxc from scratch

libxc-2.1.0.tar.gz:
	curl -OL 'http://www.tddft.org/programs/octopus/down.php?file=libxc/libxc-2.1.0.tar.gz'

libxc-2.1.0: libxc-2.1.0.tar.gz
	tar -xzvf libxc-2.1.0.tar.gz

libxc-2.1.0/src/.libs/libxc.a: libxc-2.1.0
	cd libxc-2.1.0 && CFLAGS='-fPIC' CPPFLAGS='-fPIC' FCCPP=-ffreestanding ./configure
	cd libxc-2.1.0 && $(MAKE)


