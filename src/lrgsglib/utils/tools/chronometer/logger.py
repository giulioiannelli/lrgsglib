import logging

# create a dedicated logger
logger = logging.getLogger("chronometer")
logger.setLevel(logging.INFO)
logger.propagate = False  # ← do not pass messages up to root

# attach a single StreamHandler
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
handler.setFormatter(
    logging.Formatter("%(asctime)s %(name)s [%(levelname)s] %(message)s",
                      datefmt="%Y-%m-%d %H:%M:%S")
)
logger.addHandler(handler)
