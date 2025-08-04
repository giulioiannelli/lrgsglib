from ...shared import *



class ConditionalPartitioning:
    def __init__(self, condition: Union[Callable[[Any], bool], str, Any]):
        """
        Initialize with a condition in the form of a lambda, a key string, or a direct value.

        Parameters
        ----------
        condition : Union[Callable[[Any], bool], str, Any]
            A condition as a lambda (e.g., lambda x: x < threshold), a key string (e.g., "<10"), or a direct value.
        """
        self._original = condition
        if (isinstance(condition, str) and condition and condition[0] in "<=>!"):
            self.key = condition
            self.cond_func = self.key_to_cond()
        elif isinstance(condition, Number):
            # Direct value: use its string representation as key and treat condition as equality.
            self.key = f"={condition}"
            self.cond_func = lambda x, v=condition: x == v

    def key_to_cond(self) -> Callable[[Any], bool]:
        """
        Convert the stored key representation into a condition callable.

        If the key starts with an operator (e.g. "<", "<=", "==", ">=", ">"),
        returns a lambda that applies the corresponding comparison between its argument
        and the threshold parsed from the key. Otherwise, attempts to convert the key to a number,
        returning the value.

        Returns
        -------
        Union[Callable[[Any], bool], Any]
            A lambda function representing the condition if an operator is present,
            or the direct value if not.
        """
        op_map = {
            "<=": operator.le,
            ">=": operator.ge,
            "=": operator.eq,
            "<": operator.lt,
            ">": operator.gt,
            "!=": operator.ne
        }
        for op_str in sorted(op_map.keys(), key=len, reverse=True):
            if self.key.startswith(op_str):
                value_str = self.key[len(op_str):]
                try:
                    val = float(value_str)
                    if val.is_integer():
                        val = int(val)
                except ValueError:
                    val = value_str
                return lambda x, op=op_map[op_str], v=val: op(x, v)
        try:
            val = float(self.key)
            if val.is_integer():
                val = int(val)
            return val
        except ValueError:
            return self.key

    def __repr__(self) -> str:
        return f"ConditionalPartitioning(key='{self.key}', original={self._original})"
