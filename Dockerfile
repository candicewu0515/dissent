# Dissent — a container that needs no configuration to run.
#
# No API keys and no database: every source is a public endpoint, and the one
# bundled data file (the ClinVar baseline) ships in the image. So this is a
# single build and a single run, with nothing to set up first.

FROM python:3.13-slim

# Bytecode files and a stdout buffer both make container logs lie about where a
# process got to.
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

# Dependencies first, so a code change does not re-resolve them. Only the two
# files the build backend reads are copied at this point.
COPY pyproject.toml README.md ./
COPY dissent/__init__.py dissent/__init__.py
RUN pip install --upgrade pip && pip install .

COPY . .
RUN pip install --no-deps -e .

# The cache is an optimisation that must never fail a request, and it already
# turns itself off if it cannot write. Giving it a writable directory it owns
# means it does not have to: a read-only root filesystem then costs speed, not
# correctness.
ENV DISSENT_CACHE_PATH=/var/cache/dissent/dissent.sqlite
RUN useradd --create-home --uid 10001 dissent \
    && mkdir -p /var/cache/dissent \
    && chown -R dissent:dissent /var/cache/dissent /app
USER dissent

# PaaS platforms assign the port at run time; 8000 is the fallback for a plain
# `docker run -p 8000:8000`.
ENV PORT=8000
EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD python -c "import os,urllib.request;urllib.request.urlopen(f\"http://127.0.0.1:{os.environ.get('PORT','8000')}/api/health\",timeout=4)"

# A shell is needed so $PORT expands, and `exec` is what makes that safe: it
# replaces the shell with uvicorn, so SIGTERM reaches the server rather than
# stopping at a /bin/sh that never forwards it. Without the exec a platform's
# stop signal is ignored until it gives up and sends SIGKILL, which means no
# graceful shutdown and connections cut mid-response on every deploy.
#
# One worker on purpose: every request is I/O bound on six public APIs, so the
# concurrency comes from asyncio, and a second worker would only split the
# response cache in two.
CMD ["sh", "-c", "exec uvicorn dissent.api.app:app --host 0.0.0.0 --port ${PORT:-8000}"]
