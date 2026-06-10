from functools import wraps
import contextlib
import gcad_ext

__all__ = ("autolog", "log", "logblock")

def autolog(fn):
    @wraps(fn)
    def fn_wrapper(*args, **kwargs):
        t = gcad_ext.logline_block_start(fn.__name__)
        try:
            return fn(*args, **kwargs)
        finally:
            gcad_ext.logline_block_end(t, fn.__name__)
    return fn_wrapper

def log(line):
    gcad_ext.logline(line)

@contextlib.contextmanager
def logblock(name: str):
    t = gcad_ext.logline_block_start(name)
    try:
        yield
    finally:
        gcad_ext.logline_block_end(t, name)
