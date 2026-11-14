"""Safe, auditable retrieval of public data over HTTPS.

Design rules (each is tested):

* **HTTPS only, host allow-list.** A source may only reach hosts it declares; redirects to any other
  host are refused, so a compromised or mistyped URL cannot make the app fetch arbitrary addresses.
* **Bounded.** Hard limits on download size (compressed and decompressed) and on time.
* **Auditable.** Every retrieval returns the URL, UTC timestamp, byte count and SHA-256 of the body,
  and the bytes are cached on disk with a JSON sidecar so a later run can prove what it analysed.
* **Injectable transport.** Tests supply an ``opener``; no test touches the network.

Nothing here redistributes data: files land in a local cache directory on the machine running the
app, at the user's explicit request.
"""

from __future__ import annotations

import gzip
import hashlib
import http.client
import io
import json
import os
import re
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable, Iterable
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from amr_surveil._version import __version__
from amr_surveil.exceptions import AMRSurveilError

DEFAULT_MAX_BYTES = 50 * 1024 * 1024
DEFAULT_TIMEOUT_S = 60.0
_CHUNK = 64 * 1024


class FetchError(AMRSurveilError):
    """A public-data retrieval was refused or failed (message says why and what to try)."""


@dataclass(frozen=True)
class FetchResult:
    url: str
    path: Path
    sha256: str
    n_bytes: int
    retrieved_at: str
    from_cache: bool

    def provenance(self) -> dict[str, Any]:
        record = asdict(self)
        record["path"] = str(self.path)
        return record


Opener = Callable[[urllib.request.Request, float], Any]


def _default_opener(allowed_hosts: frozenset[str]) -> Opener:
    class _Guard(urllib.request.HTTPRedirectHandler):
        def redirect_request(self, req, fp, code, msg, headers, newurl):  # type: ignore[override]
            _check_url(newurl, allowed_hosts)
            return super().redirect_request(req, fp, code, msg, headers, newurl)

    director = urllib.request.build_opener(_Guard)

    def open_(request: urllib.request.Request, timeout: float) -> Any:
        return director.open(request, timeout=timeout)

    return open_


def _check_url(url: str, allowed_hosts: Iterable[str]) -> urllib.parse.ParseResult:
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme != "https":
        raise FetchError(f"Only https:// URLs are allowed (got {parsed.scheme or 'no scheme'!r}).")
    host = (parsed.hostname or "").lower()
    if host not in {h.lower() for h in allowed_hosts}:
        raise FetchError(
            f"Host {host!r} is not on this source's allow-list {sorted(allowed_hosts)}. "
            "Nothing was downloaded."
        )
    if parsed.username or parsed.password:
        raise FetchError("URLs with embedded credentials are not allowed.")
    return parsed


def _read_limited(response: Any, limit: int) -> bytes:
    declared = None
    headers = getattr(response, "headers", None)
    if headers is not None:
        declared = headers.get("Content-Length")
    if declared is not None and str(declared).isdigit() and int(declared) > limit:
        raise FetchError(
            f"The server announces {int(declared):,} bytes, above the {limit:,}-byte limit. "
            "Narrow the query (organism, drug, years, country) or download the file yourself and upload it."
        )
    buffer = bytearray()
    while True:
        chunk = response.read(_CHUNK)
        if not chunk:
            break
        buffer.extend(chunk)
        if len(buffer) > limit:
            raise FetchError(
                f"Download exceeded the {limit:,}-byte limit and was stopped. "
                "Narrow the query or upload a file you downloaded yourself."
            )
    return bytes(buffer)


def _gunzip_limited(data: bytes, limit: int) -> bytes:
    out = bytearray()
    with gzip.GzipFile(fileobj=io.BytesIO(data)) as handle:
        while True:
            chunk = handle.read(_CHUNK)
            if not chunk:
                break
            out.extend(chunk)
            if len(out) > limit:
                raise FetchError("Decompressed data exceeds the size limit; refusing to continue.")
    return bytes(out)


def fetch(
    url: str,
    *,
    allowed_hosts: Iterable[str],
    cache_dir: str | Path,
    source_id: str,
    headers: dict[str, str] | None = None,
    max_bytes: int = DEFAULT_MAX_BYTES,
    timeout: float = DEFAULT_TIMEOUT_S,
    refresh: bool = False,
    opener: Opener | None = None,
    now: Callable[[], datetime] | None = None,
) -> FetchResult:
    """Download ``url`` (GET) to the cache and return its provenance. Cached copies are reused
    unless ``refresh`` is true, so a repeated query is reproducible and polite to the provider."""
    hosts = frozenset(h.lower() for h in allowed_hosts)
    _check_url(url, hosts)
    if not re.fullmatch(r"[A-Za-z0-9_-]{1,40}", source_id):
        raise FetchError(f"Invalid source_id {source_id!r}.")
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256(url.encode("utf-8")).hexdigest()[:24]
    body_path = cache / f"{source_id}-{key}.bin"
    meta_path = cache / f"{source_id}-{key}.json"
    if not refresh and body_path.is_file() and meta_path.is_file():
        try:  # a damaged or foreign cache entry is simply re-downloaded
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            data = body_path.read_bytes()
            if (
                isinstance(meta, dict)
                and meta.get("url") == url
                and len(data) <= max_bytes
                and hashlib.sha256(data).hexdigest() == meta.get("sha256")
            ):
                return FetchResult(
                    url,
                    body_path,
                    meta["sha256"],
                    len(data),
                    str(meta["retrieved_at"]),
                    from_cache=True,
                )
        except (ValueError, KeyError, OSError):
            pass
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": f"AMR-Surveil/{__version__} (research; public-data client)",
            "Accept-Encoding": "gzip",
            **(headers or {}),
        },
        method="GET",
    )
    transport = opener or _default_opener(hosts)
    try:
        response = transport(request, timeout)
        try:
            raw = _read_limited(response, max_bytes)
            encoding = (getattr(response, "headers", None) or {}).get("Content-Encoding", "")
        finally:
            close = getattr(response, "close", None)
            if callable(close):
                close()
    except FetchError:
        raise
    except urllib.error.HTTPError as exc:
        raise FetchError(f"{source_id}: the server answered HTTP {exc.code} for {url}.") from exc
    except (urllib.error.URLError, TimeoutError, OSError, http.client.HTTPException) as exc:
        raise FetchError(
            f"{source_id}: could not reach {urllib.parse.urlparse(url).hostname} ({exc}). "
            "Check your internet connection, proxy or firewall, or upload a file you downloaded yourself."
        ) from exc
    if encoding.lower() == "gzip" or raw[:2] == b"\x1f\x8b":
        try:
            raw = _gunzip_limited(raw, max_bytes)
        except (OSError, EOFError) as exc:
            raise FetchError(f"{source_id}: response was not valid gzip data ({exc}).") from exc
    stamp = (now() if now else datetime.now(timezone.utc)).strftime("%Y-%m-%dT%H:%M:%SZ")
    digest = hashlib.sha256(raw).hexdigest()
    body_path.write_bytes(raw)
    meta_path.write_text(
        json.dumps(
            {
                "url": url,
                "sha256": digest,
                "n_bytes": len(raw),
                "retrieved_at": stamp,
                "source_id": source_id,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return FetchResult(url, body_path, digest, len(raw), stamp, from_cache=False)


def default_cache_dir() -> Path:
    """Per-user, owner-only cache directory (``AMR_SURVEIL_CACHE`` overrides)."""
    configured = os.environ.get("AMR_SURVEIL_CACHE")
    base = (
        Path(configured) if configured else Path.home() / ".cache" / "amr_surveil" / "public_data"
    )
    try:
        base.mkdir(parents=True, exist_ok=True, mode=0o700)
        if os.name == "posix":
            base.chmod(0o700)
    except OSError:  # read-only home (some containers): fall back to a private temp directory
        import tempfile

        base = Path(tempfile.mkdtemp(prefix="amr_surveil_cache_"))
    return base
