"""Optional Model Context Protocol integration. Importing this module needs no SDK."""


def create_server():
    """Create the local server, loading optional dependencies on demand."""
    try:
        from .server import create_server as factory
    except ImportError as error:
        raise ImportError("MCP support requires: pip install 'paftacular[mcp]'") from error
    return factory()
