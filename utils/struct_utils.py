from functools import wraps
from typing import Any, Callable, Union, get_type_hints

def format_fn(in_f):
    @wraps(in_f)
    def str_fn(*args, **kwargs):
        out_float = in_f(*args, **kwargs)
        return f"{out_float:.4f}"

    str_fn.__name__ = f"{in_f.__name__}_str"
    return str_fn


def format_property(in_property: property) -> property:
    old_fn = in_property.fget
    new_fn = format_fn(old_fn)
    return property(new_fn)


def does_return_float(in_attr, in_fn_attr=None) -> bool:
    if in_fn_attr is None and callable(in_attr):
        test_fn = in_attr
    elif in_fn_attr is not None and hasattr(in_attr, in_fn_attr):
        test_fn = getattr(in_attr, in_fn_attr)
    else:
        return False
    return get_type_hints(test_fn).get("return", None) is float


def get_format_factory(in_attr: Any) -> Union[Callable, None]:
    if isinstance(in_attr, property) and does_return_float(in_attr, "fget"):
        return format_property
    elif callable(in_attr) and does_return_float(in_attr):
        return format_fn
    return None

_meta_formatter_cache: dict[str, dict[str, Callable]] = {}
class MetaFormatter(type):
    def __new__(cls, name, bases, attrs):
        
        new_fns: dict[str, Callable] = _meta_formatter_cache.get(name, dict())
        if not new_fns:
            for attr_key, attr in attrs.items():
                if attr_key.startswith("__"):
                    continue

                format_factory = get_format_factory(attr)
                if format_factory is not None:
                    new_fns[f"{attr_key}_str"] = format_factory(attr)
            _meta_formatter_cache[name] = new_fns
        
        attrs.update(new_fns)
        return super(MetaFormatter, cls).__new__(cls, name, bases, attrs)
