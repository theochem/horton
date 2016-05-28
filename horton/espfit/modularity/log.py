try:
    from horton.log import log
except ImportError:
    print("Warning, using built-in logger")
    class Log(object):
        def cite(self, *args):
            print(" ".join(args))

    log = Log()