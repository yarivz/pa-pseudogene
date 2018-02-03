import logging
from logging.handlers import QueueHandler


def listener_configurer():
    root = logging.getLogger()
    formatter = logging.Formatter('[%(levelname)s %(asctime)s %(name)s] - %(message)s')
    file_handler = logging.FileHandler('pseudogene.log', 'w')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    console_handler.setFormatter(formatter)
    root.addHandler(file_handler)
    root.addHandler(console_handler)


def worker_configurer(queue):
    q_handler = QueueHandler(queue)
    root = logging.getLogger()
    root.addHandler(q_handler)
    root.setLevel(logging.DEBUG)


def listener_process(queue, configurer):
    configurer()
    while True:
        try:
            record = queue.get()
            if record is None:
                break
            logger = logging.getLogger(record.name)
            logger.handle(record)
        except Exception:
            import sys, traceback
            print('Whoops! Problem:', file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
