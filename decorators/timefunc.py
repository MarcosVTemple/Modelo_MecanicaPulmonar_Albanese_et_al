import time
from datetime import timedelta


def timefunc(func):
    def inner(*args, **kwargs):
        start = time.perf_counter()
        func(*args, **kwargs)
        end = time.perf_counter()
        delta = timedelta(seconds=end - start)
        if delta.seconds >= 60:
            minutos, segundos = divmod(delta.seconds, 60)
            print(f"A função {func.__name__} tomou {minutos} minutos e {segundos} segundos")
        elif delta.seconds >= 1:
            segundos = delta.seconds
            milissegundos = delta.microseconds / 1000 if delta.microseconds > 1000 else 0
            print(f"A função {func.__name__} tomou {segundos},{int(milissegundos)} segundos")
        elif delta.microseconds >= 1000:
            milissegundos, microssegundos = divmod(delta.microseconds, 1000)
            print(f"A função {func.__name__} tomou {int(milissegundos)} milissegundos")
        else:
            microssegundos = delta.microseconds
            print(f"A função {func.__name__} tomou {int(microssegundos)} microssegundos")
    return inner
