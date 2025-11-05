from typing import Any, Callable, Union, Optional
import operator
from numbers import Number

__all__ = ["ConditionalPartitioning", "_normalize_conditional_partitioning"]

class ConditionalPartitioning:
    def __init__(self, condition: Union[Callable[[Any], bool], str, Any]):
        """
        Initialize with a condition in the form of a lambda, a key string, or a direct value.

        Parameters
        ----------
        condition : Union[Callable[[Any], bool], str, Any]
            A condition as a lambda (e.g., lambda x: x < threshold), a key string (e.g., "<10"), or a direct value.
            
        Raises
        ------
        ValueError
            If condition is None or an unsupported type.
        TypeError
            If condition is a callable but doesn't appear to be a valid condition function.
        """
        if condition is None:
            raise ValueError(
                "ConditionalPartitioning condition cannot be None. "
                "Please provide a number, string (e.g., '>0', '=1'), or callable."
            )
        
        self._original = condition
        
        # Handle string conditions (operators like "<10", ">=5", etc.)
        if isinstance(condition, str):
            if not condition:
                raise ValueError("ConditionalPartitioning condition string cannot be empty.")
            if condition[0] in "<=>!":
                self.key = condition
                self.cond_func = self.key_to_cond()
            else:
                raise ValueError(
                    f"Invalid condition string '{condition}'. "
                    f"String conditions must start with an operator: <, >, =, <=, >=, !=. "
                    f"Examples: '>0', '=1', '<=10'"
                )
        # Handle numeric conditions
        elif isinstance(condition, Number):
            # Direct value: use its string representation as key and treat condition as equality.
            self.key = f"={condition}"
            self.cond_func = lambda x, v=condition: x == v
        # Handle callable conditions
        elif callable(condition):
            # Try to create a reasonable key for callables
            if hasattr(condition, '__name__'):
                self.key = f"callable:{condition.__name__}"
            else:
                self.key = f"callable:{id(condition)}"
            self.cond_func = condition
            # Validate that the callable can accept at least one argument
            try:
                import inspect
                sig = inspect.signature(condition)
                # Check if there are any positional or keyword parameters
                params = [p for p in sig.parameters.values() 
                         if p.kind in (inspect.Parameter.POSITIONAL_OR_KEYWORD,
                                      inspect.Parameter.POSITIONAL_ONLY,
                                      inspect.Parameter.VAR_POSITIONAL)]
                if len(params) == 0:
                    raise TypeError(
                        f"Callable condition must accept at least one argument. "
                        f"Expected signature like 'lambda x: ...' or 'def func(x): ...', "
                        f"but got a callable that accepts no positional arguments."
                    )
            except (ValueError, TypeError) as e:
                # If inspection fails or raises our TypeError, re-raise
                if isinstance(e, TypeError) and "must accept at least one argument" in str(e):
                    raise
                # Otherwise, can't inspect (e.g., built-in), so try calling it
                # with a test value to see if it works
                try:
                    _ = condition(0)
                except TypeError as call_err:
                    raise TypeError(
                        f"Callable condition appears to not accept the expected arguments. "
                        f"Error when testing: {call_err}"
                    ) from call_err
        else:
            raise TypeError(
                f"Unsupported condition type: {type(condition).__name__}. "
                f"ConditionalPartitioning accepts: numbers, strings (e.g., '>0'), "
                f"or callables (e.g., lambda x: x > 0). Got: {condition}"
            )
        

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


def _normalize_conditional_partitioning(
    val: Union["ConditionalPartitioning", Number, str, Callable[[Any], bool], None]
) -> Optional["ConditionalPartitioning"]:
    """
    Normalize various input types to a ConditionalPartitioning object.
    
    This function provides a convenient way to convert different types of inputs
    into ConditionalPartitioning objects, making it easier to work with conditional
    partitioning across the codebase.
    
    Parameters
    ----------
    val : Union[ConditionalPartitioning, Number, str, Callable, None]
        Input value that can be:
        - A ConditionalPartitioning object (returned as-is)
        - A number (int/float) - converted to equality condition (e.g., 1 → "=1")
        - A string with operator (e.g., '>0', '<=10', '!=5')
        - A callable (e.g., lambda x: x > 0)
        - None (returned as None)
    
    Returns
    -------
    Optional[ConditionalPartitioning]
        Normalized ConditionalPartitioning object or None if input was None.
        
    Raises
    ------
    ValueError
        If val is an invalid string (doesn't start with an operator or is empty).
    TypeError
        If val is not a supported type (e.g., list, dict, tuple) or if a callable
        doesn't accept the expected number of arguments.
        
    Examples
    --------
    >>> # From a number
    >>> cp = normalize_conditional_partitioning(1)
    >>> cp.key
    '=1'
    >>> cp.cond_func(1)
    True
    
    >>> # From a string operator
    >>> cp = normalize_conditional_partitioning('>0')
    >>> cp.key
    '>0'
    >>> cp.cond_func(5)
    True
    
    >>> # From a lambda
    >>> cp = normalize_conditional_partitioning(lambda x: x > 0)
    >>> cp.cond_func(5)
    True
    
    >>> # Already a ConditionalPartitioning (returns same object)
    >>> original = ConditionalPartitioning('>0')
    >>> result = normalize_conditional_partitioning(original)
    >>> result is original
    True
    
    >>> # None input
    >>> normalize_conditional_partitioning(None)
    None
    """
    if val is None:
        return None
    if isinstance(val, ConditionalPartitioning):
        return val
    if isinstance(val, Number):
        return ConditionalPartitioning(val)
    # Try to convert strings and callables - let ConditionalPartitioning validate
    if isinstance(val, str) or callable(val):
        return ConditionalPartitioning(val)
    # If we get here, it's an unsupported type
    raise TypeError(
        f"Cannot normalize value of type {type(val).__name__} to ConditionalPartitioning. "
        f"Expected: ConditionalPartitioning, Number, string (e.g., '>0'), callable, or None. "
        f"Got: {val}"
    )
