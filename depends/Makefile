default:
	@echo 'Choose one of the following targets'
	@echo '* libint: downloads and compiles libint-2.0.0-stable.tgz for mpqc'
	@echo '* libxc: downloads and compiles libxc-1.2.0'

libint: libint-2.0.0-stable/lib/libint2.a

libxc: libxc-1.2.0/src/.libs/libxc.a


libint-2.0.0-stable/lib/libint2.a: libint-2.0.0-stable
	cd libint-2.0.0-stable && CFLAGS='-fPIC' CPPFLAGS='-fPIC' ./configure --with-cxx-optflags='-O0' && $(MAKE)

libint-2.0.0-stable: libint-2.0.0-stable.tgz
	tar -xvzf libint-2.0.0-stable.tgz

libint-2.0.0-stable.tgz:
	wget 'http://downloads.sourceforge.net/project/libint/libint-for-mpqc/libint-2.0.0-stable.tgz'


libxc-1.2.0.tar.gz:
	wget 'http://www.tddft.org/programs/octopus/down.php?file=libxc/libxc-1.2.0.tar.gz' -O libxc-1.2.0.tar.gz

libxc-1.2.0: libxc-1.2.0.tar.gz
	tar -xvzf libxc-1.2.0.tar.gz

libxc-1.2.0/src/.libs/libxc.a: libxc-1.2.0
	cd libxc-1.2.0 && CFLAGS='-fPIC' CPPFLAGS='-fPIC' ./configure && $(MAKE)


