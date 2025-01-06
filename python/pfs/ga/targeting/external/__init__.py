# Import and patch broken external libraries here

import logging
from contextlib import contextmanager

@contextmanager
def preserve_logging():
    original_handlers = logging.root.handlers[:]
    original_level = logging.root.level
    original_formatters = [handler.formatter for handler in logging.root.handlers]

    yield

    logging.root.handlers = original_handlers
    logging.root.level = original_level
    for handler, formatter in zip(logging.root.handlers, original_formatters):
        handler.setFormatter(formatter)

# Prevent cobraCharmer to hijack the logging
with preserve_logging():
    from ics.cobraOps.BlackDotsCalibrationProduct import BlackDotsCalibrationProduct
    from ics.cobraCharmer.pfiDesign import PFIDesign
    from ics.cobraCharmer.cobraCoach.cobraCoach import CobraCoach
    from ics.cobraCharmer.pfi import PFI