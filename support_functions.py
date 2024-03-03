# this script include some support functions for LeafcutterITI

import time
from functools import wraps


def timing_decorator(func):
    # the timing decorator, use to keep track the running time of functions
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time

        # Initialize an empty list to collect time parts
        time_parts = []
        
        # Calculate hours, minutes, and seconds
        hours, remainder = divmod(elapsed_time, 3600)
        minutes, seconds = divmod(remainder, 60)
        
        # Append non-zero time parts to the list
        if hours:
            time_parts.append(f"{int(hours)}h")
        if minutes:
            time_parts.append(f"{int(minutes)}m")
        # Always include seconds
        time_parts.append(f"{seconds:.2f}s")
        
        # Join the non-zero parts with a space
        time_str = " ".join(time_parts)
        
        print(f"{func.__name__!r} executed in {time_str}")
        return result
    return wrapper


def write_options_to_file(options, filename):
    with open(filename, 'w') as file:
        for attr in vars(options):
            value = getattr(options, attr)
            file.write(f"{attr}: {value}\n")



