default:
	@echo 'Choose one of the following targets'
	@echo '* libint: downloads and compiles libint-2.0.0-stable.tgz for mpqc'
	@echo '* libint-linux-x86_64-gcc-4.7.2: downloads and unpacks a binary libint for linux-x86_64 compiled with gcc-4.7.2'
	@echo '* libxc: downloads and compiles libxc-2.0.1'

libint: libint-2.0.0-stable/lib/libint2.a

libxc: libxc-2.0.1/src/.libs/libxc.a


# Build libint from scratch

libint-2.0.0-stable/lib/libint2.a: libint-2.0.0-stable
	cd libint-2.0.0-stable && CFLAGS='-fPIC' CPPFLAGS='-fPIC' ./configure --with-cxx-optflags='-O0'
	cd libint-2.0.0-stable && $(MAKE)

libint-2.0.0-stable: libint-2.0.0-stable.tgz
	tar -xzvf libint-2.0.0-stable.tgz

libint-2.0.0-stable.tgz:
	wget 'http://downloads.sourceforge.net/project/libint/libint-for-mpqc/libint-2.0.0-stable.tgz' -O libint-2.0.0-stable.tgz


# Binary compilation (mainly for Travis CI)

libint-2.0.0-bin-horton-linux-x86_64-gcc-4.7.2.tar.bz2:
	wget http://users.ugent.be/~tovrstra/horton/libint-2.0.0-bin-horton-linux-x86_64-gcc-4.7.2.tar.bz2

libint-linux-x86_64-gcc-4.7.2: libint-2.0.0-bin-horton-linux-x86_64-gcc-4.7.2.tar.bz2
	tar -xjf libint-2.0.0-bin-horton-linux-x86_64-gcc-4.7.2.tar.bz2


# Build libxc from scratch

libxc-2.0.1.tar.gz:
	wget 'http://www.tddft.org/programs/octopus/down.php?file=libxc/libxc-2.0.1.tar.gz' -O libxc-2.0.1.tar.gz

libxc-2.0.1: libxc-2.0.1.tar.gz
	tar -xzvf libxc-2.0.1.tar.gz
	patch -d libxc-2.0.1/src < libxc-2.0.1-fix-initialization.patch

libxc-2.0.1/src/.libs/libxc.a: libxc-2.0.1
	cd libxc-2.0.1 && CFLAGS='-fPIC' CPPFLAGS='-fPIC' ./configure
	cd libxc-2.0.1 && $(MAKE)


