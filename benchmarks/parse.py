"""Run with Python from the environment whose installation is being measured."""

import argparse
import statistics
import time
import tracemalloc

import paftacular as p


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count", type=int, default=10000)
    args = parser.parse_args()
    cases = {
        "simple": "y5",
        "modified": "&2@y5-H2O+i13C[M+H+Na]^2/-0.55ppm*0.85",
        "multiple": "b2,y3-H2O,p^2",
    }
    print(f"Package: {p.__file__}")
    for name, text in cases.items():
        samples = []
        for _ in range(5):
            start = time.perf_counter()
            for _ in range(args.count):
                p.parse_multi(text)
            samples.append((time.perf_counter() - start) / args.count * 1e6)
        print(f"{name}: {statistics.median(samples):.2f} microseconds per record")
    tracemalloc.start()
    for index in range(args.count * 2):
        p.parse_single(f"r[benchmark-{index}]")
    retained, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    print(f"Unique references: {len(p.ReferenceIon._cache)} cached, {retained / 1024:.0f} KiB retained, {peak / 1024:.0f} KiB peak")


if __name__ == "__main__":
    main()
