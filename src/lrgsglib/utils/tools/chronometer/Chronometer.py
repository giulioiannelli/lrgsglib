import time
from .logger import logger
from ...basic import generate_random_id

class Chronometer:
    """
    High-precision timer for measuring and accumulating execution durations.

    Attributes
    ----------
    name : str
        Logical name of this chronometer.
    id : str
        Unique identifier.
    start_time : float
        Timestamp when timing began.
    end_time : float | None
        Timestamp when timing ended.
    auto_log : bool
        If True, logs elapsed time upon stop.

    Class Attributes
    ----------------
    _data : dict
        Accumulated { name: {total_time, count} } for all instances.
    enabled : bool
        If True, chronometers are active and will log timings.
    """
    _data = {}
    enabled: bool = False

    def __init__(self, name: str, auto_log: bool=False) -> None:
        self.name = name
        self.id = generate_random_id()
        self.start_time = time.perf_counter()
        self.end_time = None
        self.auto_log = auto_log

    def end(self) -> None:
        """Stop timing, accumulate data and (optionally) log the result."""
        if self.end_time is not None:
            if self.auto_log:
                logger.warning(
                    "Chronometer '%s' (ID=%s) already stopped.", self.name, self.id
                )
            return

        self.end_time = time.perf_counter()
        elapsed = self.get_elapsed_time()
        entry = self._data.setdefault(self.name, {"total_time": 0.0, "count": 0})
        entry["total_time"] += elapsed
        entry["count"] += 1

        if self.auto_log:
            logger.info(
                "Stopped '%s' (ID=%s): %.4g s", self.name, self.id, elapsed
            )

    def get_elapsed_time(self) -> float:
        """
        Return elapsed time in seconds.
        If not stopped, returns time since start.
        """
        end = self.end_time or time.perf_counter()
        return end - self.start_time

    @classmethod
    def summary(cls) -> None:
        """
        Log a table of accumulated timings, sorted by descending average duration.
        """
        logger.info("Summary of all chronometers:")
        rows = []
        for name, d in cls._data.items():
            avg = d["total_time"] / d["count"]
            rows.append((name, d["total_time"], d["count"], avg))
        rows.sort(key=lambda x: x[3], reverse=True)
        logger.info(
            "%-25s %-15s %-7s %-15s",
            "Name", "Total Time (s)", "Calls", "Average (s)"
        )
        for name, total, count, avg in rows:
            logger.info(
                "%-25s %-15.4g %-7d %-15.4g",
                name, total, count, avg
            )

    @classmethod
    def enable(cls) -> None:
        """
        Enable all chronometers, allowing them to log timings.

        Parameters:
        -----------
        None

        Returns:
        -----------
        None

        Notes:
        -----------
        Sets the class attribute `enabled` to True.
        """
        cls.enabled = True


# def _generate_id(length: int = 8) -> str:
#     """Return a random alphanumeric ID."""
#     import random, string
#     return ''.join(random.choices(string.ascii_letters + string.digits, k=length))