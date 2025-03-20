def safe_deep_copy(orig, **kwargs):
    if orig is list:
        return [ safe_deep_copy(i, **kwargs) for i in orig ]
    elif orig is tuple:
        return tuple([ safe_deep_copy(i, **kwargs) for i in orig ])
    elif orig is dict:
        return { k: safe_deep_copy(i, **kwargs) for k, i in orig.items() }
    elif orig is not None and hasattr(orig, 'copy'):
        return orig.copy(**kwargs)
    elif isinstance(orig, str):
        return orig
    elif orig is None:
        return None
    else:
        raise ValueError("Cannot make a copy of object.")
    