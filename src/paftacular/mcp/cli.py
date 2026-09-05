"""Local stdio entry point with no mandatory MCP imports."""

import argparse
import sys

from paftacular import __version__


def main() -> None:
    parser = argparse.ArgumentParser(description="Expose paftacular tools, resources, and prompts over MCP stdio.")
    parser.add_argument("--version", action="version", version=f"paftacular-mcp {__version__}")
    parser.parse_args()
    from . import create_server

    try:
        server = create_server()
    except ImportError as error:
        parser.exit(2, f"{error}\n")
    try:
        server.run(transport="stdio")
    except KeyboardInterrupt:
        sys.exit(0)
