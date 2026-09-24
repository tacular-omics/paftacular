import os

try:
    from hypothesis import HealthCheck, settings
except ImportError:  # the wheel-installation CI jobs install only pytest
    settings = None

if settings is not None:
    # The default profile keeps CI fast. Run `HYPOTHESIS_PROFILE=thorough uv run pytest` for a deep search.
    settings.register_profile("default", max_examples=150, deadline=None, suppress_health_check=[HealthCheck.too_slow])
    settings.register_profile("thorough", max_examples=5000, deadline=None, suppress_health_check=[HealthCheck.too_slow])
    settings.load_profile(os.environ.get("HYPOTHESIS_PROFILE", "default"))
