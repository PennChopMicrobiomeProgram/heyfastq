# syntax=docker/dockerfile:1

FROM python:3.11-slim

# Set up the working directory for the application
WORKDIR /app

# Copy project metadata first to leverage Docker layer caching
COPY pyproject.toml README.md /app/

# Copy the source code into the image
COPY src /app/src

# Install the package and its CLI entrypoint
RUN pip install --no-cache-dir .

# Default command prints usage instructions; override to run other subcommands
ENTRYPOINT ["heyfastq"]
CMD ["--help"]
