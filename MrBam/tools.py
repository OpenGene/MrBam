def try_append(d, k, v):
    if k in d:
        d[k].append(v)
    else:
        d[k] = [v]

def memo(f):
    class memodict(dict):
        def __missing__(self, key):
            ret = self[key] = f(key)
            return ret
    return memodict().__getitem__
