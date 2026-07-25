FROM python:3.12-slim

# fuse3 + libfuse-dev needed for the optional VFS mount feature.
# build-essential not required -- everything here is pure Python.
RUN apt-get update && apt-get install -y --no-install-recommends \
        fuse3 \
        libfuse-dev \
        pkg-config \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY pyproject.toml README.md ./
COPY dnastore ./dnastore

# Install the package itself plus optional extras. If you don't need
# the FUSE mount feature, drop the [vfs] extra to skip that dependency.
RUN pip install --no-cache-dir -e ".[vfs,dev]"

COPY tests ./tests
COPY examples ./examples
COPY benchmarks ./benchmarks
COPY smoke_test.py ./

# Install the demo plugin package too, so smoke_test.py / pytest exercise
# the plugin-discovery path inside the container as well (matches CI).
RUN pip install --no-cache-dir --no-deps -e examples/plugin_example

# FUSE mounts inside a container need the device and appropriate
# capabilities; run with:
#   docker run --rm -it --cap-add SYS_ADMIN --device /dev/fuse dnastore
RUN mkdir -p /mnt/dna

CMD ["python", "smoke_test.py"]
