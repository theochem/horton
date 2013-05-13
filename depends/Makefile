default:
	@echo 'Choose one of the following targets'
	@echo '* libint: downloads and compiles libint-2.0.0-stable.tgz for mpqc'
	@echo '* libxc: downloads and compiles libxc-2.0.1'

libint: libint-2.0.0-stable/lib/libint2.a

libxc: libxc-2.0.1/src/.libs/libxc.a


libint-2.0.0-stable/lib/libint2.a: libint-2.0.0-stable
	cd libint-2.0.0-stable && CFLAGS='-fPIC' CPPFLAGS='-fPIC' ./configure --with-cxx-optflags='-O0' &> configure.out && echo 'Completed: configure libint'
	cd libint-2.0.0-stable && $(MAKE) &> make.out && echo 'Completed: make libint'

libint-2.0.0-stable: libint-2.0.0-stable.tgz
	tar -xvzf libint-2.0.0-stable.tgz &> untar-libint.out && echo 'Completed: untar libint'

libint-2.0.0-stable.tgz:
	wget 'http://downloads.sourceforge.net/project/libint/libint-for-mpqc/libint-2.0.0-stable.tgz'


libxc-2.0.1.tar.gz:
	wget 'http://www.tddft.org/programs/octopus/down.php?file=libxc/libxc-2.0.1.tar.gz' -O libxc-2.0.1.tar.gz

libxc-2.0.1: libxc-2.0.1.tar.gz
	tar -xvzf libxc-2.0.1.tar.gz &> libxc-untar.out && echo 'Completed: untar libxc'
	patch -d libxc-2.0.1/src < libxc-2.0.1-fix-initialization.patch

libxc-2.0.1/src/.libs/libxc.a: libxc-2.0.1
	cd libxc-2.0.1 && CFLAGS='-fPIC' CPPFLAGS='-fPIC' ./configure &> configure.out && echo 'Completed: configure libxc'
	cd libxc-2.0.1 && $(MAKE) &> make.out && echo 'Completed: make libxc'


