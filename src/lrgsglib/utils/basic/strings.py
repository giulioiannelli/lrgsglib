from .common import *

def generate_random_id(
    length: int = 10,
    chars: str = string.ascii_letters + string.digits
) -> str:
    """
    Generate a random string of specified length from the given character set,
    using Python's built-in `random` module.

    Parameters
    ----------
    length : int
        Desired length of the output string (must be non-negative).
    chars : str
        Characters to sample from (must be non-empty).

    Returns
    -------
    str
        Randomly generated string of given length.

    Raises
    ------
    ValueError
        If `length` is negative or `chars` is empty.
    """
    if length < 0:
        raise ValueError("`length` must be non-negative.")
    if not chars:
        raise ValueError("`chars` must be a non-empty string.")

    return ''.join(random.choice(chars) for _ in range(length))
#
def get_first_int_in_str(s):
    """
    Finds the first sequence of digits in the given string and returns it as an integer.
    If no digits are found, returns None.

    Parameters:
    ---------------
    s : str
        The string to search for digits.

    Returns:
    ---------------
    int or None
        The first integer found in the string, or None if no digits are found.
    """
    match = re.search(r'\d+', s)
    return int(match.group()) if match else None
#
def join_non_empty(separator: str, *parts: str) -> str:
    """
    Join non-empty string parts with a given separator.

    Parameters
    ----------
    separator : str
        The string to insert between non-empty parts.
    parts : str
        Variable number of string parts to join.

    Returns
    -------
    str
        The joined string of all non-empty parts. Returns an empty string
        if no parts are non-empty.

    Raises
    ------
    TypeError
        If `separator` or any of the `parts` is not a string.
    """
    if not isinstance(separator, str):
        raise TypeError("`separator` must be a string")

    non_empty_parts = []
    for idx, part in enumerate(parts, start=1):
        if not isinstance(part, str):
            raise TypeError(f"Argument #{idx} is not a string: {part!r}")
        if part:
            non_empty_parts.append(part)

    return separator.join(non_empty_parts)