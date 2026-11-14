"""Domain-specific exceptions used by AMR-Surveil."""


class AMRSurveilError(Exception):
    """Base exception for expected application failures."""


class InputFileError(AMRSurveilError):
    """Raised when an input file cannot be read or has an invalid schema."""


class QualityGateError(AMRSurveilError):
    """Raised when strict quality controls detect blocking problems."""


class ConfigurationError(AMRSurveilError):
    """Raised when configuration values are unsupported or inconsistent."""
