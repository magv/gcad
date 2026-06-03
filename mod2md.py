#!/usr/bin/env python3
"""
Format the contents of a module as Markdown.
"""

import importlib
import inspect
import sys
import types

def strip_indent(text: str) -> str:
    lines = text.splitlines()
    common_n = None
    for line in lines:
        line_lstrip = line.lstrip()
        if line_lstrip != "":
            n = len(line) - len(line_lstrip)
            common_n = n if common_n is None else min(common_n, n)
    return "\n".join(line[common_n:] for line in lines).strip()

# First pass

HREFS = {}
IDS = {}

def prewr_function(name: str, fun: types.FunctionType, mod: types.ModuleType):
    anchor = f"{mod.__name__}.{name}"
    HREFS[f"{fun.__module__}.{name}"] = f"#{anchor}"
    IDS[f"{fun.__module__}.{name}"] = anchor

def prewr_class(name: str, typ: type, mod: types.ModuleType):
    anchor = f"{mod.__name__}.{name}"
    HREFS[f"{typ.__module__}.{name}"] = f"#{anchor}"
    IDS[f"{typ.__module__}.{name}"] = anchor

def prewr_module(mod: types.ModuleType):
    for name, value in inspect.getmembers(mod):
        if name.startswith("_"): continue
        if inspect.isfunction(value):
            prewr_function(name, value, mod)
        elif inspect.isclass(value):
            prewr_class(name, value, mod)
        else:
            pass

# Second pass

def wr(text: str):
    sys.stdout.write(text)

def wr_type(typ):
    if isinstance(typ, types.GenericAlias):
        wr_type(typ.__origin__)
        wr("[")
        for i, arg in enumerate(typ.__args__):
            if i != 0: wr(", ")
            wr_type(arg)
        wr("]")
    elif isinstance(typ, types.UnionType):
        for i, arg in enumerate(typ.__args__):
            if i != 0: wr(" | ")
            wr_type(arg)
    elif isinstance(typ, type):
        href = HREFS.get(f"{typ.__module__}.{typ.__name__}")
        if href:
            wr(f"[{typ.__name__}]({href})")
        else:
            wr(f"{typ.__name__}")
    else:
        raise ValueError(f"Can't format type {typ} (of type {type(typ)})")

def wr_function(name: str, fun: types.FunctionType, mod: types.ModuleType):
    anchor = IDS.get(f"{fun.__module__}.{fun.__name__}")
    if anchor:
        wr(f"<a id=\"{anchor}\"></a>\n")
    wr(f"**{name}**(")
    ann = inspect.get_annotations(fun, eval_str=True)
    doc = fun.__doc__
    for i, (k, v) in enumerate(ann.items()):
        if k == "return": continue
        if i != 0: wr(", ")
        wr(f"{k}: *")
        wr_type(v)
        wr(f"*")
    wr(") → *")
    wr_type(ann['return'])
    wr("*\n\n")
    for line in strip_indent(doc).splitlines():
        wr(f"> {line}\n")
    wr("\n")

def wr_dataclass(name: str, typ: type, mod: types.ModuleType):
    anchor = IDS.get(f"{typ.__module__}.{typ.__name__}")
    if anchor:
        wr(f"<a id=\"{anchor}\"></a>\n")
    wr(f"dataclass **{name}**(")
    for i, (fld, fld_type) in enumerate(inspect.get_annotations(typ, eval_str=True).items()):
        if i != 0: wr(", ")
        wr(f"{fld}: *")
        wr_type(fld_type)
        wr("*")
    wr(f")\n\n")

def wr_class(name: str, typ: type, mod: types.ModuleType):
    anchor = IDS.get(f"{typ.__module__}.{typ.__name__}")
    if anchor:
        wr(f"<a id=\"{anchor}\"></a>\n")
    wr(f"class **{name}**\n\n")

def wr_module(mod: types.ModuleType):
    for name, value in inspect.getmembers(mod):
        if name.startswith("_"): continue
        if inspect.isfunction(value):
            wr_function(name, value, mod)
    for name, value in inspect.getmembers(mod):
        if name.startswith("_"): continue
        if inspect.isfunction(value):
            pass
        elif inspect.isclass(value):
            if hasattr(value, "__dataclass_fields__"):
                wr_dataclass(name, value, mod)
            else:
                wr_class(name, value, mod)
        elif inspect.ismodule(value):
            pass
        else:
            raise ValueError(f"Can't document {mod.__name__}.{name} (a {type(value)})")

if __name__ == "__main__":

    mod = importlib.import_module("gcad")
    prewr_module(mod)
    wr_module(mod)
