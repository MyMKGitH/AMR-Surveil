# syntax=docker/dockerfile:1
# AMR-Surveil: Streamlit dashboard + CLI. Hugging Face Docker Space compatible (port 7860, non-root user).
FROM python:3.12-slim

ARG VERSION=1.1.1
LABEL org.opencontainers.image.title="AMR-Surveil" \
      org.opencontainers.image.description="Reproducible genomic and epidemiological AMR surveillance workbench" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.source="https://github.com/MyMKGitH/amr-surveil" \
      org.opencontainers.image.licenses="MIT"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_SERVER_ADDRESS=0.0.0.0 \
    STREAMLIT_SERVER_PORT=7860 \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false \
    HOME=/home/appuser

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && useradd --create-home --uid 1000 --shell /usr/sbin/nologin appuser

# Dependencies first so code edits do not invalidate the pip layer.
COPY requirements.txt ./
RUN python -m pip install --upgrade pip && python -m pip install -r requirements.txt

COPY pyproject.toml README.md LICENSE ./
COPY amr_surveil ./amr_surveil
COPY app.py ./app.py
COPY .streamlit ./.streamlit
COPY config ./config
COPY data ./data
COPY schemas ./schemas
RUN python -m pip install --no-deps . \
    && amr-surveil --version

USER appuser
EXPOSE 7860

HEALTHCHECK --interval=30s --timeout=5s --start-period=30s --retries=3 \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:7860/_stcore/health', timeout=3)" || exit 1

CMD ["streamlit", "run", "app.py", "--server.port=7860", "--server.address=0.0.0.0", "--server.headless=true"]
