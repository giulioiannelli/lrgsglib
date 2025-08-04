class NestedDict(dict):
    def __getitem__(self, key):
        if key not in self:
            self[key] = NestedDict()
        return super().__getitem__(key)

    def __setitem__(self, key, value):
        if isinstance(value, dict) and not isinstance(value, NestedDict):
            value = NestedDict(value)
        super().__setitem__(key, value)

    def __delitem__(self, key):
        super().__delitem__(key)

    def get(self, key, default=None):
        return super().get(key, default)

    def set(self, key, value):
        self[key] = value

    def delete(self, key):
        del self[key]

    def to_dict(self):
        def recursive_convert(d):
            if isinstance(d, NestedDict):
                d = {k: recursive_convert(v) for k, v in d.items()}
            return d
        return recursive_convert(self)