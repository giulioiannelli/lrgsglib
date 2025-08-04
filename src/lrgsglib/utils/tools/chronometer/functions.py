from .Chronometer import Chronometer

def print_accumulated_timings() -> None:
    """
    Log elapsed time for each live chronometer instance.
    """
    for name, chrono in Chronometer._data.items():
        elapsed = chrono.get_elapsed_time()
        logger = __import__("chronometer.logger", fromlist=["logger"]).logger
        logger.info(
            "Function '%s' (ID=%s): %.4g s",
            name, getattr(chrono, "id", "?"), elapsed
        )