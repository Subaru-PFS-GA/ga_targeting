class ReadOnlyList():
    """
    Implements a read-only wrapper around a python list.
    """

    def __init__(self, l: list):
        self.__l = l

    def __getitem__(self, i):
        return self.__l.__getitem__(i)

    def __setitem__(self, i, item):
        raise RuntimeError("The list is read-only.")
	
    def __delitem__(self, i):
        raise RuntimeError("The list is read-only.")

    def __contains__(self, o):
        return self.__l.__contains__(o)

    def __format__(self, format_spec):
        return self.__l.__format__(format_spec)

    def __len__(self):
        return self.__l.__len__()

    def __iter__(self):
        return self.__l.__iter__()

    def __str__(self):
        return self.__l.__str__()
	
    def clear(self):
        raise RuntimeError("The list is read-only.")

    def copy(self):
        return self.__l.copy()
	
    def pop(self, key, *args):
        raise RuntimeError("The list is read-only.")

    def append(self, item):
        raise RuntimeError("The list is read-only.")

    def extend(self, list):
        raise RuntimeError("The list is read-only.")
