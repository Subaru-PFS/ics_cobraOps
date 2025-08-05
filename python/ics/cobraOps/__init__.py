import importlib.metadata

try:
    __version__ = importlib.metadata.version("ics.cobraOps")
except importlib.metadata.PackageNotFoundError:
    __version__ = "unknown"
