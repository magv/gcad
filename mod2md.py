#!/usr/bin/env python3
"""
Format the contents of a module as Markdown.
"""

import ast
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

def class_member_docstrings(typ: type) -> dict[str, str]:
    docstrings = {}
    tree = ast.parse(inspect.getsource(typ))
    prev_el = None
    for el in tree.body[0].body:
        match (prev_el, el):
            case (
                ast.AnnAssign(target=ast.Name(id=name)),
                ast.Expr(value=ast.Constant(comment)),
            ):
                docstrings[name] = comment
        prev_el = el
    return docstrings

# First pass

HREFS = {}
IDS = {}

def is_cython_function(obj) -> bool:
    return type(obj).__name__ == "cython_function_or_method"

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
        if inspect.isfunction(value) or is_cython_function(value):
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
    elif typ == types.NoneType:
        wr(f"None")
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
    wr(f"**{name}** (")
    ann = inspect.get_annotations(fun, eval_str=True)
    for i, (k, v) in enumerate(ann.items()):
        if k == "return": continue
        if i != 0: wr(", ")
        wr(f"{k}: *")
        wr_type(v)
        wr(f"*")
    wr(") → *")
    wr_type(ann['return'])
    wr("*\n\n")
    for line in strip_indent(fun.__doc__).splitlines():
        wr(f"{line}\n")
    wr("\n")

def wr_dataclass(name: str, typ: type, mod: types.ModuleType):
    docstrings = class_member_docstrings(typ)
    anchor = IDS.get(f"{typ.__module__}.{typ.__name__}")
    if anchor:
        wr(f"<a id=\"{anchor}\"></a>\n")
    doc = strip_indent(typ.__doc__)
    wr(f"**{name}**. {doc}\n\n")
    for i, (fld, fld_type) in enumerate(inspect.get_annotations(typ, eval_str=True).items()):
        wr(f"* **{fld}**: *")
        wr_type(fld_type)
        wr("*")
        doc = docstrings.get(fld)
        if doc:
            wr(f". {doc.strip()}")
        wr("\n")
    wr(f"\n")

def wr_class(name: str, typ: type, mod: types.ModuleType):
    anchor = IDS.get(f"{typ.__module__}.{typ.__name__}")
    if anchor:
        wr(f"<a id=\"{anchor}\"></a>\n")
    wr(f"class **{name}**\n\n")

def wr_module(mod: types.ModuleType):
    wr("### Functions\n\n")
    for name, value in inspect.getmembers(mod):
        if name.startswith("_"): continue
        if inspect.isfunction(value) or is_cython_function(value):
            wr_function(name, value, mod)
    wr("### Classes\n\n")
    for name, value in inspect.getmembers(mod):
        if name.startswith("_"): continue
        if inspect.isclass(value):
            if hasattr(value, "__dataclass_fields__"):
                wr_dataclass(name, value, mod)
            else:
                wr_class(name, value, mod)
    # Double-check there's nothing else.
    for name, value in inspect.getmembers(mod):
        if name.startswith("_"): continue
        if inspect.isfunction(value) or is_cython_function(value): continue
        if inspect.isclass(value): continue
        if inspect.ismodule(value): continue
        raise ValueError(f"Can't document {mod.__name__}.{name} (a {type(value)})")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} module_name", file=sys.stderr)
        exit(1)
    mod = importlib.import_module(sys.argv[1])
    prewr_module(mod)
    wr_module(mod)
