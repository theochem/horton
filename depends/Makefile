default:
	@echo 'Choose one of the following targets'
	@echo '* libint: downloads and compiles libint-2.0.0-stable.tgz for mpqc'

libint: libint-2.0.0-stable/lib/libint2.a

libint-2.0.0-stable/lib/libint2.a: libint-2.0.0-stable
	cd libint-2.0.0-stable && CFLAGS='-fPIC' CPPFLAGS='-fPIC' ./configure && $(MAKE)

libint-2.0.0-stable: libint-2.0.0-stable.tgz
	tar -xvzf libint-2.0.0-stable.tgz

libint-2.0.0-stable.tgz:
	wget http://downloads.sourceforge.net/project/libint/libint-for-mpqc/libint-2.0.0-stable.tgz


