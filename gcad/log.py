from contextlib import contextmanager
from functools import wraps
from gcad_ext import (
    log,
    log_trace_enter,
    log_trace_exit,
    trace_enter,
    trace_exit,
    trace_progress,
    trace_text,
)

def auto_trace(fn):
    @wraps(fn)
    def fn_wrapper(*args, **kwargs):
        trace_enter(fn.__name__)
        try:
            return fn(*args, **kwargs)
        finally:
            trace_exit()
    return fn_wrapper

def auto_log_trace(fn):
    @wraps(fn)
    def fn_wrapper(*args, **kwargs):
        log_trace_enter(fn.__name__)
        try:
            return fn(*args, **kwargs)
        finally:
            log_trace_exit()
    return fn_wrapper

@contextmanager
def trace(text: str):
    trace_enter(text)
    try:
        yield
    finally:
        trace_exit()

@contextmanager
def log_trace(text: str):
    log_trace_enter(text)
    try:
        yield
    finally:
        log_trace_exit()
