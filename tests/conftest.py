import os

os.environ["RAVEN_KN_BROKERS"] = "0"     # no request to the alert brokers while a transient is assessed
