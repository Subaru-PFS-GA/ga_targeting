class ReadOnlyDict():
    """
    Implements a read-only wrapper around a python dict.
    """

    def __init__(self, d: dict):
        self.__d = d

    def __getitem__(self, k):
        return self.__d.__getitem__(k)

    def __setitem__(self, key, item):
        raise RuntimeError("The dictionary is read-only.")
	
    def __delitem__(self, key):
        raise RuntimeError("The dictionary is read-only.")

    def __contains__(self, o):
        return self.__d.__contains__(o)

    def __format__(self, format_spec):
        return self.__d.__format__(format_spec)

    def __len__(self):
        return self.__d.__len__()

    def __iter__(self):
        return self.__d.__iter__()

    def __str__(self):
        return self.__d.__str__()
	
    def clear(self):
        raise RuntimeError("The dictionary is read-only.")

    def copy(self):
        return self.__d.copy()
    
    def items(self):
        return self.__d.items()

    def keys(self):
        return self.__d.keys()
	
    def pop(self, key, *args):
        raise RuntimeError("The dictionary is read-only.")

    def popitem(self):
        raise RuntimeError("The dictionary is read-only.")

    def update(self, dict=None):
        raise RuntimeError("The dictionary is read-only.")

    def values(self):
        return self.__d.values()
