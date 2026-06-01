import sys
import time
import contextlib

__all__ = ("autolog", "log", "logblock")

class Logger:
    def __init__(self):
        self._write = sys.stderr.write
        self._flush = sys.stderr.flush
        self._time = time.time
        now = self._time()
        self.first_log_time = now
        self.last_log_time = now
        self.stack = []
        self.empty_block = False
    def _pr(self, line: str, now: float):
        t = now - self.first_log_time
        dt = now - self.last_log_time
        self._write(f"{t:.3f} +{dt:.3f} {'│'*len(self.stack)}")
        self._write(line)
        self._write("\n")
        self._flush()
        self.last_log_time = now
    def log(self, line: str):
        """Print a line to the log."""
        self._pr(line, self._time())
        self.empty_block = False
    def push_block(self, name: str):
        now = self._time()
        self._pr(f"╭{name}", now)
        self.stack.append((now, name))
        self.empty_block = True
    def pop_block(self):
        now = self._time()
        then, name = self.stack.pop()
        if self.empty_block:
            self._write("\033[F\033[K")
            self._pr(f"-{name}: {now - then:.2e}s", now)
        else:
            self._pr(f"╰{name}: {now - then:.2e}s", now)
        self.empty_block = False

LOG = Logger()

def autolog(fn):
    def fn_wrapper(*args, **kwargs):
        LOG.push_block(fn.__name__)
        try:
            return fn(*args, **kwargs)
        finally:
            LOG.pop_block()
    fn_wrapper.__name__ = fn.__name__
    fn_wrapper.__doc__ = fn.__doc__
    return fn_wrapper

def log(line):
    LOG.log(line)

@contextlib.contextmanager
def logblock(name):
    LOG.push_block(name)
    try:
        yield
    finally:
        LOG.pop_block()
