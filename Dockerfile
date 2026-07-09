FROM python:3.11-slim

WORKDIR /app

# Scientific dependencies ship as manylinux wheels, so no system build tools are
# needed. The package itself is mounted at /app for development, or added for a
# release build.
RUN pip install --no-cache-dir \
    "duckdb>=1.0" "polars>=1.0" "pandas>=2.0" "numpy>=1.26" "scipy>=1.11" \
    "pydantic>=2.6" "cma>=3.3" "anthropic>=0.116" "matplotlib>=3.8" "pytest>=8.0"

ENV PYTHONPATH=/app
CMD ["pytest", "-q"]
