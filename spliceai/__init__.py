# Original source code modified to add prediction batching support by Invitae in 2021.
# Modifications copyright (c) 2021 Invitae Corporation.

import signal
from pkg_resources import get_distribution


try:
    signal.signal(signal.SIGINT, lambda x, y: exit(0))
except ValueError:
    # Continue if we're not able to set the signal handler due to which thread is running the code
    pass

name = 'spliceai'
__version__ = get_distribution(name).version
