import argparse
from typing import List, Tuple, Union

class IntOrTriple(argparse.Action):
    """
    Store one integer as ``int`` or three integers as
    ``tuple[int, int, int]``.
    """

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: List[int],
        option_string: str | None = None,
    ) -> None:
        """
        Validate *values* and assign the parsed representation to *namespace*.

        Parameters
        ----------
        parser
            The argument parser invoking this action.
        namespace
            Namespace where the parsed value is stored.
        values
            Sequence of integers supplied on the command line.
        option_string
            The flag triggering this action (unused for positional arguments).
        """
        if len(values) == 1:
            parsed: Union[int, Tuple[int, int, int]] = values[0]
        elif len(values) == 3:
            parsed = tuple(values)  # type: ignore[assignment]
        else:
            raise argparse.ArgumentError(
                self, "expected either 1 or 3 integers"
            )

        setattr(namespace, self.dest, parsed)