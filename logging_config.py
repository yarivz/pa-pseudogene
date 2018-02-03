import logging
from logging.handlers import QueueHandler


def listener_configurer():
    root = logging.getLogger()
    formatter = logging.Formatter('[%(levelname)s %(asctime)s %(module)s] - %(message)s')
    file_handler = logging.FileHandler('pseudogene.log', 'w')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    root.addHandler(file_handler)
    root.addHandler(console_handler)


def worker_configurer(queue):
    q_handler = QueueHandler(queue)  # Just the one handler needed
    root = logging.getLogger()
    root.addHandler(q_handler)
    # send all messages, for demo; no other level or filter logic applied.
    root.setLevel(logging.DEBUG)


def listener_process(queue, configurer):
    configurer()
    while True:
        try:
            record = queue.get()
            if record is None:  # We send this as a sentinel to tell the listener to quit.
                break
            logger = logging.getLogger(record.name)
            logger.handle(record)  # No level or filter logic applied - just do it!
        except Exception:
            import sys, traceback
            print('Whoops! Problem:', file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
